#!/usr/bin/env python3
"""Build and test the production Fortran eigensolvers in a temporary directory.

Requires a GNU-compatible Fortran compiler, NumPy and SciPy. FC and LAPACK_LIBS
may override the compiler and link options. Does not build/install LAPW1 or
modify calculation directories. The default link options are Accelerate on
macOS and -llapack -lblas elsewhere.
"""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[1]
WORK = Path(tempfile.mkdtemp(prefix="dftks_solver_tests_"))
ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1")
for key in list(ENV):
    if key.startswith("DFTKS_EIG") or key.startswith("DFTKS_DAVIDSON"):
        ENV.pop(key)
RESULTS = []


def build(real):
    directory = WORK / ("real" if real else "complex")
    directory.mkdir()
    compiler = shlex.split(os.environ.get("FC", "gfortran"))
    libs = shlex.split(os.environ.get("LAPACK_LIBS",
                      "-framework Accelerate" if sys.platform == "darwin"
                      else "-llapack -lblas"))
    command = compiler + ["-cpp", "-ffree-form", "-fopenmp", "-O2", "-g", "-fcheck=all",
                          "-ffree-line-length-none"]
    if real:
        command += ["-D_REAL"]
    command += [str(ROOT / "tests/eigensolver_driver.F90"),
                str(ROOT / "diagonalize.FPP"), "-o", "solver"] + libs
    proc = subprocess.run(command, cwd=directory, env=ENV,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        raise RuntimeError(proc.stdout + proc.stderr)
    return directory / "solver"


def fixture(real, n=80, condition=30.0):
    rng = np.random.default_rng(3751)
    def orthogonal():
        x = rng.normal(size=(n, n))
        if not real:
            x = x + 1j * rng.normal(size=(n, n))
        return np.linalg.qr(x)[0]
    a, u = orthogonal(), orthogonal()
    root = (a * np.sqrt(np.geomspace(1., condition, n))) @ a.conj().T
    spectrum = np.linspace(-4., 8., n)
    h = root @ (u * spectrum) @ u.conj().T @ root
    o = root @ root
    return (h + h.conj().T) / 2, (o + o.conj().T) / 2


def solve(exe, h, o, name, method, cache, *, fixed=False, bands=8,
          lo=-2.7, hi=-2., capacity=None, labels_shift=0, kshift=0.,
          extra=None, fails=False):
    n = len(h)
    real = np.isrealobj(h)
    dtype = np.dtype("<f8" if real else "<c16")
    raw = np.zeros_like(h)
    lower = np.tril_indices(n, -1)
    upper = np.triu_indices(n)
    raw[lower] = h[lower].conj()
    raw[upper] = o[upper]
    labels = np.vstack((np.arange(n) + labels_shift, np.zeros(n), np.ones(n))).astype("<i4")
    inp, out = WORK / (name + ".input"), WORK / (name + ".output")
    with inp.open("wb") as stream:
        stream.write(np.array([n, capacity or n, int(fixed), bands, 0], dtype="<i4").tobytes())
        stream.write(np.array([lo, hi, .125 + kshift, .25, .375], dtype="<f8").tobytes())
        stream.write(labels.tobytes(order="F"))
        stream.write(raw.astype(dtype).tobytes(order="F"))
        stream.write(np.diag(h).astype(dtype).tobytes())
    env = dict(ENV, DFTKS_EIGENSOLVER=method, **(extra or {}))
    proc = subprocess.run([str(exe), str(inp), str(out), str(cache)],
                          cwd=WORK, env=env, text=True, capture_output=True)
    if fails:
        assert proc.returncode != 0, proc.stdout
        RESULTS.append(dict(name=name, expected_error=True))
        return proc.stdout
    assert proc.returncode == 0, proc.stdout + proc.stderr
    with out.open("rb") as stream:
        nn, m = np.fromfile(stream, dtype="<i4", count=2)
        assert nn == n
        e = np.fromfile(stream, dtype="<f8", count=m)
        z = np.fromfile(stream, dtype=dtype, count=n*m).reshape((n, m), order="F")
    ref = eigh(h, o, eigvals_only=True)
    ref = ref[:bands] if fixed else ref[(ref > lo) & (ref <= hi)]
    assert len(ref) == len(e), (name, len(ref), len(e), proc.stdout)
    error = np.max(np.abs(ref-e), initial=0.)
    residuals = h @ z - (o @ z)*e
    scale = np.maximum(1., np.linalg.norm(h@z, axis=0) + np.abs(e)*np.linalg.norm(o@z, axis=0))
    residual = np.max(np.linalg.norm(residuals, axis=0)/scale, initial=0.)
    metric = np.max(np.abs(z.conj().T @ o @ z - np.eye(m)), initial=0.)
    assert error < 2e-8, (name, error, proc.stdout)
    assert residual < 2e-9, (name, residual, proc.stdout)
    assert metric < 2e-8, (name, metric, proc.stdout)
    result = dict(name=name, method=method, bands=int(m), eigenvalue_error=float(error),
                  relative_residual=float(residual), metric_error=float(metric),
                  log=proc.stdout)
    RESULTS.append(result)
    return proc.stdout


def suite(real):
    tag = "real" if real else "complex"
    exe = build(real)
    h, o = fixture(real)
    cache = WORK / (tag + "cache")
    def run(name, method="GVD", hh=h, oo=o, **kwargs):
        return solve(exe, hh, oo, tag+"_"+name, method, cache, **kwargs)
    run("legacy", "GVX")
    run("dense_window")
    run("dense_fixed", fixed=True)
    run("empty_window", lo=20., hi=21.)
    log = run("seed", "DAVIDSON", fixed=True)
    assert "dense seed" in log
    log = run("identical", "DAVIDSON", fixed=True)
    assert ":DFTKS_SOLVER: DAVIDSON" in log, log
    rng = np.random.default_rng(41)
    perturbation = rng.normal(size=h.shape)
    if not real:
        perturbation = perturbation+1j*rng.normal(size=h.shape)
    perturbation = (perturbation+perturbation.conj().T)*.0001
    log = run("warm", "DAVIDSON", hh=h+perturbation, fixed=True)
    assert ":DFTKS_SOLVER: DAVIDSON" in log, log
    log = run("moving_overlap", "DAVIDSON", hh=h+perturbation,
              oo=o+perturbation*.1, fixed=True)
    assert ":DFTKS_SOLVER: DAVIDSON" in log, log
    # Reuse the same cache but deliberately force a failed iteration.
    log = run("fallback", "DAVIDSON", hh=h-perturbation, fixed=True,
              extra={"DFTKS_DAVIDSON_MAXITER":"1"})
    assert "dense fallback" in log, log
    log = run("label_change", "DAVIDSON", fixed=True, labels_shift=1)
    assert "dense seed" in log
    log = run("k_change", "DAVIDSON", fixed=True, kshift=.01)
    assert "dense seed" in log
    # A truncated cache must be ignored rather than producing an I/O crash.
    cachefile = Path(str(cache) + f".dftks-{'r' if real else 'c'}.k00000001")
    cachefile.write_bytes(b"DFTKSE01")
    log = run("truncated", "DAVIDSON", fixed=True)
    assert "dense seed" in log
    # Large requested subspace is routed to the dense solver in AUTO.
    run("many_seed", "DAVIDSON", fixed=True, bands=50)
    log = run("many_auto", "AUTO", fixed=True, bands=50)
    assert "band fraction" in log
    run("capacity", capacity=1, fails=True)
    run("bad_solver", "typo", fails=True)
    run("bad_setting", "DAVIDSON", extra={"DFTKS_DAVIDSON_TOL":"nan"}, fails=True)
    bad_o=o.copy(); bad_o[0,0]=-1000.
    run("indefinite_metric", oo=bad_o, fails=True)
    log = run("tiny_space", "DAVIDSON", hh=h+perturbation, fixed=True,
              extra={"DFTKS_DAVIDSON_MAXSPACE":"8"})
    assert "dense fallback" in log
    # Diagonal problems make hidden band crossings particularly dangerous:
    # old states can have identically zero residual after an omitted band drops.
    hd = np.diag(np.arange(60, dtype=float))
    od = np.eye(60)
    if not real:
        hd=hd.astype(complex); od=od.astype(complex)
    run("cross_seed", "DAVIDSON", hh=hd, oo=od, fixed=True)
    moved = hd.copy(); moved[-1,-1]=-10.
    log = run("hidden_cross", "DAVIDSON", hh=moved, oo=od, fixed=True)
    assert "inertia" in log and "dense fallback" in log, log
    degenerate=hd.copy();degenerate[4,4]=3.;degenerate[5,5]=3.
    run("degenerate_seed", "DAVIDSON", hh=degenerate, oo=od, fixed=True)
    run("degenerate_warm", "DAVIDSON", hh=degenerate, oo=od, fixed=True)
    run("window_seed", "DAVIDSON", fixed=False)
    log = run("window_warm", "DAVIDSON", hh=h+perturbation)
    assert ":DFTKS_SOLVER: DAVIDSON" in log, log
    # Force the shift onto an exact eigenvalue; factorization fails safely.
    log = run("singular_shift", "DAVIDSON", hh=hd, oo=od, lo=-1., hi=8.)
    log = run("singular_shift_warm", "DAVIDSON", hh=hd, oo=od, lo=-1., hi=8.)
    assert "dense fallback" in log
    hi, oi = fixture(real, condition=1e5)
    run("conditioned_dense", hh=hi, oo=oi, fixed=True)
    run("conditioned_seed", "DAVIDSON", hh=hi, oo=oi, fixed=True)
    run("conditioned_warm", "DAVIDSON", hh=hi, oo=oi, fixed=True)


if __name__ == "__main__":
    try:
        suite(True)
        suite(False)
    finally:
        (WORK / "results.json").write_text(json.dumps(RESULTS, indent=2))
        print(f"Results: {WORK / 'results.json'}")
    print(f"PASS: {len(RESULTS)} cases (real and complex)")
