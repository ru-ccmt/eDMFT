#!/usr/bin/env python3
"""Execute prepared dmft2, then an optional frozen-Sigma response.

The callers own input preparation and sigma splitting. This module never
regenerates those inputs. run_dmft defers response until a single final pass;
a standalone x_dmft execution may respond immediately. The present response
analytic backend saves a complete tensor for publication during SCF assembly.
The optional legacy native backend exports density/potential tangents only.
"""
from __future__ import annotations

from contextlib import contextmanager
import csv
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import uuid

ALIASES = {"STRT": "STR", "STRFT": "STRF"}
SOLVER_KEYS = {
    "backend",
    "physical_policy",
    "lapw0_worker",
    "lapw1_worker",
    "lcore_worker",
    "library",
    "charge_integration",
    "strain",
    "probe_step",
    "tolerance",
    "max_iterations",
    "mpi",
    "native_mpi",
    "cache_mb",
    "basis_cache_mb",
    "dmft2",
    "wien_root",
}
SNAPSHOT_SUFFIXES = {
    ".struct",
    ".in0",
    ".in1",
    ".in1c",
    ".in2",
    ".in2c",
    ".inc",
    ".inm",
    ".klist",
    ".kgen",
    ".rotlm",
    ".indmfl",
    ".indmfi",
    ".indmf2",
    ".indmfldn",
    ".indmf2dn",
    ".clmsum",
    ".vsp",
    ".vns",
}


def entries(text):
    """Read DEF filenames structurally; never replace arbitrary '_x' text."""
    result = {}
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = next(csv.reader([line], quotechar="'", skipinitialspace=True))
        if len(fields) < 5:
            raise ValueError(f"invalid DEF record: {line}")
        unit = int(fields[0])
        if unit in result:
            raise ValueError(f"duplicate DEF unit {unit}")
        result[unit] = fields[1].strip()
    return result


def replace_filenames(text, changes):
    lines = []
    for line in text.splitlines(True):
        match = re.match(r"(\s*(\d+)\s*,\s*)'([^']+)'", line)
        if match and int(match[2]) in changes:
            name = str(changes[int(match[2])])
            if "'" in name or "\n" in name:
                raise ValueError("DEF filenames cannot contain quotes or newlines")
            line = line[: match.start(3)] + name + line[match.end(3) :]
        lines.append(line)
    return "".join(lines)


def local_definition(text):
    """Isolated response runs use local files and unsplit WIEN vectors.

    Preserve case names containing '_x'. Only the output-file suffix '_x'
    is removed. All output paths are rebased into the response workspace;
    this prevents inherited absolute DEF paths overwriting the parent run.
    """
    names = {}
    for unit, filename in entries(text).items():
        name = Path(filename).name
        if name in {".", "..", ""}:
            raise ValueError(f"invalid DEF filename {filename!r}")
        names[unit] = name[:-2] if name.endswith("_x") else name
    return replace_filenames(text, names)


def input_mode(path):
    lines = Path(path).read_text().splitlines()
    if not lines or not lines[0].split():
        raise ValueError(f"empty IN2 mode: {path}")
    return lines[0].split()[0].upper()


def charge_mode(path):
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            return line.split()[0].lower()
    raise ValueError(f"empty dmft2 input: {path}")


def request(case, definition="dmft2.def", config_path=None, disabled=False):
    """Read live controls; an explicit disabled setting overrides IN2 aliases."""
    defs = entries(Path(definition).read_text())
    mode = input_mode(defs[5])
    path = Path(config_path) if config_path else Path(case + ".stress_response.json")
    if config_path and not path.is_file():
        raise ValueError(f"response configuration does not exist: {path}")
    config = json.loads(path.read_text()) if path.is_file() else {}
    if not isinstance(config, dict):
        raise ValueError("response configuration must be a JSON object")
    unknown = set(config) - {"enabled", "solver", "workspace_root"}
    if unknown:
        raise ValueError(f"unknown response settings: {sorted(unknown)}")
    enabled = config.get("enabled", mode in ALIASES)
    if not isinstance(enabled, bool):
        raise ValueError("response enabled must be true or false")
    return {
        "enabled": enabled and not disabled,
        "mode": mode,
        "config": config,
        "config_path": str(path.resolve()) if path.is_file() else None,
        "definition": str(Path(definition).resolve()),
        "files": defs,
    }


def requested(case, definition="dmft2.def"):
    return request(case, definition)["enabled"]


def atomic_json(path, record):
    path = Path(path)
    with tempfile.NamedTemporaryFile(
        mode="w", dir=path.parent, prefix=path.name + ".", delete=False
    ) as stream:
        temporary = Path(stream.name)
        try:
            json.dump(record, stream, indent=2, allow_nan=False)
            stream.write("\n")
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(path)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@contextmanager
def execution_definition(definition):
    """Translate response aliases for native executables without editing IN2.

    Native WIEN/dmft2 readers have four-character modes. STRFT is therefore
    a wrapper request, not a new Fortran input format. Temporary definitions
    are in the original working directory so other relative paths retain
    their meaning, and are removed even on subprocess failure.
    """
    text = Path(definition).read_text()
    source = Path(entries(text)[5])
    mode = input_mode(source)
    if mode not in ALIASES:
        yield str(definition)
        return
    with tempfile.TemporaryDirectory(prefix=".dmft2-response-input-", dir=".") as temp:
        folder = Path(temp)
        native = folder / source.name
        lines = source.read_text().splitlines(True)
        lines[0] = re.sub(r"^\s*\S+", f"{ALIASES[mode]:4s}", lines[0], count=1)
        native.write_text("".join(lines))
        df = folder / "dmft2.def"
        df.write_text(replace_filenames(text, {5: native}))
        yield str(df)


def run_core(root, mpi, definition, output):
    """Fail on nonzero exits AND Fortran STOP errors that return zero.

    Keep arbitrary site MPI launch prefixes, as the historical launchers do;
    quote executable and DEF paths. Stream merged output so users see progress.
    """
    command = (
        (mpi + " " if mpi else "")
        + shlex.quote(str(Path(root) / "dmft2"))
        + " "
        + shlex.quote(str(definition))
    )
    print("DMFT2 ..... running:", command, file=output, flush=True)
    finished = False
    with subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    ) as child:
        for line in child.stdout:
            print(line, end="", file=output, flush=True)
            if "DMFT2 END" in line:
                finished = True
        code = child.wait()
    if code:
        raise RuntimeError(f"dmft2 exited with status {code}")
    errors = [
        str(p)
        for p in Path(".").glob("dmft2*.error*")
        if p.is_file() and p.stat().st_size
    ]
    if errors or not finished:
        raise RuntimeError(
            "dmft2 did not finish successfully: "
            + ", ".join(errors or ["missing DMFT2 END"])
        )


def preflight(case, req, root, mpi):
    """Reject unsupported physics before starting an expensive response."""
    config = req["config"]
    solver = config.get("solver", {})
    if not isinstance(solver, dict) or set(solver) - SOLVER_KEYS:
        raise ValueError(
            "unsupported solver settings; see dmft2_response_integration.txt"
        )
    backend = solver.get("backend", "analytic")
    if backend not in {"native", "analytic"}:
        raise ValueError("response backend must be native or analytic")
    if solver.get("charge_integration", "native") not in {"native", "legacy_spline"}:
        raise ValueError("charge_integration must be native or legacy_spline")
    defs = req["files"]
    for suffix in [".inso", ".clmup", ".clmdn", ".incup", ".incdn"]:
        p = Path(case + suffix)
        if p.is_file() and p.stat().st_size:
            raise ValueError("response requires no SOC and an unpolarized KS/core potential")
    control = [
        line.split("#", 1)[0].split()
        for line in Path(defs[2]).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ][0]
    if len(control) < 3 or control[0] != "c" or int(control[2]) not in (0, 1):
        raise ValueError(
            "response requires charge mode c and recomputeEF=0 or 1"
        )
    if req["mode"] not in {"STR", "STRF", "STRT", "STRFT"}:
        raise ValueError(
            "response requires STR, STRF, STRT or STRFT in the actual IN2 input"
        )
    if not Path("lapw1.def").is_file():
        raise ValueError("response requires a prepared lapw1.def")
    from wstruct import Struct
    import numpy as np

    s = Struct()
    s.ReadStruct(case + ".struct")
    from lapw_lattice import NativeLattice
    from lapw_response_policy import detect_policy
    NativeLattice.from_structure(s)
    policy = detect_policy(".", case, control=defs[2], override=solver.get("physical_policy"))
    if backend == "native" and (policy.magnetic or not policy.fixed_number or s.lattyp.strip() != "P" or not np.allclose([s.alpha,s.beta,s.gamma],90)):
        raise ValueError("legacy native backend requires nonmagnetic fixed-N primitive orthogonal inputs")
    if policy.magnetic:
        if Path(defs[7]).name != case+".indmfl":
            raise ValueError("launch the primary dmft2 definition; response handles both spins")
        companion = Path("dmft2dn.def")
        if not companion.is_file():
            raise ValueError("magnetic response requires prepared dmft2dn.def")
        down = entries(companion.read_text())
        for unit in [5,9,18,20,30,9919]:
            if unit in defs and (unit not in down or Path(defs[unit]).resolve()!=Path(down[unit]).resolve()):
                raise ValueError("spin branches must share native KS inputs at unit "+str(unit))
        dh = Path(down[7]).read_text().splitlines()[0].split()
        if len(dh)<4 or int(dh[2])!=1 or abs(int(dh[3]))!=5:
            raise ValueError("second spin requires normalized projector 5 and fixed bands")
        for name in ["combineud"]+(["fermif"] if policy.ferro and policy.fixed_number else []):
            helper = Path(root)/name
            if not helper.is_file() or not os.access(helper,os.X_OK):
                raise ValueError("missing magnetic helper: "+str(helper))
    header = Path(defs[7]).read_text().splitlines()[0].split()
    if len(header) < 4 or int(header[2]) != 1 or abs(int(header[3])) != 5:
        raise ValueError(
            "response requires normalized projector 5 and fixed band indices"
        )
    wi = (
        Path(solver.get("wien_root", os.environ.get("WIENROOT", "")))
        .expanduser()
        .resolve()
    )
    exe = (
        Path(solver.get("dmft2", str(Path(root) / "dmft2_response")))
        .expanduser()
        .resolve()
    )
    # Inherit the ordinary dmft2 MPI prefix automatically. Python's mpi4py
    # and the Fortran workers may link different MPI implementations, so the
    # shared resolver verifies/adapts the Python launcher without changing
    # the ordinary dmft2 command. Explicit JSON argv lists still override it.
    from response_mpi import resolve_launches
    launch, native_launch = resolve_launches(mpi, solver, analytic=backend == "analytic")
    for p in [wi / "x", wi / "lapw1", wi / "lapw1c", exe]:
        if not p.is_file() or not os.access(p, os.X_OK):
            raise ValueError(f"missing executable: {p}")
    for name in [case + ".clmsum", "EF.dat", "projectorw.dat"]:
        if not Path(name).is_file():
            raise ValueError(f"missing response input: {name}")
    if backend == "analytic":
        cache_mb = solver.get("cache_mb", 256.0)
        if (
            not isinstance(cache_mb, (float, int))
            or not np.isfinite(cache_mb)
            or cache_mb < 0
        ):
            raise ValueError("solver.cache_mb must be a nonnegative per-rank budget")
        basis_cache_mb = solver.get("basis_cache_mb", 128.0)
        if (
            not isinstance(basis_cache_mb, (float, int))
            or not np.isfinite(basis_cache_mb)
            or basis_cache_mb < 0
        ):
            raise ValueError(
                "solver.basis_cache_mb must be a nonnegative per-rank budget"
            )
        legacy = set(solver) & {"strain", "probe_step", "tolerance", "max_iterations"}
        if legacy:
            raise ValueError(
                "native-reference solver controls do not apply to analytic stress: "
                + ", ".join(sorted(legacy))
                + "; remove them or select backend=native"
            )
        if not Path(case + ".in0").is_file() or input_mode(case + ".in0") != "STR":
            raise ValueError("analytic stress publication requires STR in case.in0")
        if solver.get("charge_integration", "native") != "native":
            raise ValueError("analytic backend requires native charge integration")
        from lapw_total_stress import invariant_strains

        invariant_strains(s)
        names = {
            "lapw0_worker": "lapw0_response",
            "lcore_worker": "lcore_response",
            "lapw1_worker": (
                "lapw1c_response" if req["files"][5].endswith("c") else "lapw1_response"
            ),
            "library": "liblapw_response.so",
        }
        solver = dict(solver)
        for key, name in names.items():
            path = Path(solver.get(key, str(Path(root) / name))).expanduser().resolve()
            if not path.is_file() or (
                key != "library" and not os.access(path, os.X_OK)
            ):
                raise ValueError(f"missing analytic response worker: {path}")
            solver[key] = str(path)
    return solver | {
        "backend": backend,
        "physical_policy": policy.record(),
        "wien_root": str(wi),
        "dmft2": str(exe),
        "mpi": launch,
        "native_mpi": native_launch,
        "charge_integration": solver.get("charge_integration", "native"),
    }


def compute_response(case, req, root, mpi, run_id, output):
    """Snapshot the final physical inputs, solve once, and publish a sidecar.

    No cached correction is applied to another charge iterate. The status
    file names this exact snapshot. Solver children invoke executables
    directly, avoiding recursive x_dmft/run_dmft response execution.
    """
    solver = preflight(case, req, root, mpi)
    if solver.get("backend") == "analytic":
        for label, key in [("Python", "mpi"), ("native", "native_mpi")]:
            command = solver.get(key, [])
            print("RESPONSE ..... " + label + " MPI: " +
                  (shlex.join(command) if command else "serial"),
                  file=output, flush=True)
    base = (
        Path(req["config"].get("workspace_root", f"{case}.response_runs"))
        .expanduser()
        .resolve()
    )
    base.mkdir(parents=True, exist_ok=True)
    work = base / run_id
    work.mkdir()
    snapshot = work / "reference"
    snapshot.mkdir()
    fingerprints = {}
    for p in sorted(Path(".").iterdir()):
        if p.is_file() and (
            (p.name.startswith(case + ".") and p.suffix in SNAPSHOT_SUFFIXES)
            or p.name in {"EF.dat", "projectorw.dat", "Eorb_imp.dat", "params.dat"}
            or re.fullmatch(r"sig\.inp\d*(?:dn)?", p.name)
        ):
            shutil.copy2(p, snapshot / p.name)
            fingerprints[p.name] = sha256(snapshot / p.name)
    # Actual DEF input names must match the case-local snapshot expected by
    # the backend. Nonstandard aliases should be explicit rather than silently
    # replaced by similarly named files from another calculation.
    for unit, suffix in [
        (2, ".indmf2"),
        (5, ".in2c" if req["files"][5].endswith("c") else ".in2"),
        (7, ".indmfl"),
        (20, ".struct"),
    ]:
        if Path(req["files"][unit]).resolve() != Path(case + suffix).resolve():
            raise ValueError(f"unsupported response DEF alias at unit {unit}")
    for name, source in [
        ("dmft2.def", Path(req["definition"])),
        ("lapw1.def", Path("lapw1.def")),
    ] + ([("dmft2dn.def", Path("dmft2dn.def"))] if solver["physical_policy"]["magnetic"] else []):
        (snapshot / name).write_text(local_definition(source.read_text()))
        fingerprints[name] = sha256(snapshot / name)
    # Preserve identities of split and unsplit vector/energy inputs without
    # copying those large files. The response builds its own central vectors.
    for unit in [9, 30]:
        name = req["files"].get(unit)
        if name:
            for p in sorted(Path(name).parent.glob(Path(name).name + "*")):
                if p.is_file():
                    fingerprints["parent/" + str(p)] = sha256(p)
    identities = {
        str(Path(root) / "dmft2"): sha256(Path(root) / "dmft2"),
        solver["dmft2"]: sha256(solver["dmft2"]),
    }
    for name in ["lapw0", "lapw1", "lapw1c", "lcore"]:
        p = Path(solver["wien_root"]) / name
        if p.is_file():
            identities[str(p)] = sha256(p)
    if solver.get("backend", "native") == "analytic":
        for key in ["lapw0_worker", "lapw1_worker", "lcore_worker", "library"]:
            identities[solver[key]] = sha256(Path(solver[key]))
    scripts = {}
    for module in [
        "response_parallel",
        "response_mpi",
        "lapw_response_policy",
        "dmft2_spin",
        "lapw_magnetic_response",
        "lapw_lattice",
        "frozen_sigma_tangent",
        "frozen_sigma_response",
        "lapw_valence_state",
        "lapw_screened_state",
        "lapw_geometry_response",
        "lapw_strain_charge",
        "lapw_strain_functional",
        "lapw_response_kpool",
        "lapw_response_prepare",
        "lapw_total_stress_driver",
        "lapw_total_stress",
        "lapw_klein_functional",
        "lapw_density_functional",
        "lapw_energy_response",
        "lapw_radial_energy_policy",
        "lapw_native_potential_export",
        "stress_response_publish",
        "dmft2_step",
        "wien_frozen_sigma_response",
        "frozen_sigma_screening",
        "wien_response_density",
        "lapw_charge_response",
        "lapw_hartree_periodic",
        "lapw_hartree_radial",
        "lapw_density_response",
        "stress_response_adjoint",
    ]:
        spec = importlib.util.find_spec(module)
        if spec is not None and spec.origin:
            scripts[spec.origin] = sha256(spec.origin)
    identity = {
        "inputs": fingerprints,
        "executables": identities,
        "scripts": scripts,
        "python": sys.version,
        "solver": solver,
    }
    ref = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    atomic_json(work / "reference.json", {"reference_id": ref, **identity})
    if solver.get("backend", "native") == "analytic":
        cfg = {
            "case": case,
            "physical_policy": solver["physical_policy"],
            "fermif": str(Path(root)/"fermif"),
            "reference": str(snapshot),
            "workspace": str(work / "central"),
            "wien_root": solver["wien_root"],
            "dmft2": solver["dmft2"],
            "native_mpi": solver.get("native_mpi", []),
            "lapw1_mpi": solver.get("mpi", []),
            "lapw1_definition": str(snapshot / "lapw1.def"),
            "dmft2_definition": str(snapshot / "dmft2.def"),
            **{k: solver[k] for k in ["lapw0_worker", "lapw1_worker", "lcore_worker"]},
        }
        if solver["physical_policy"]["magnetic"]:
            cfg["dmft2dn_definition"] = str(snapshot/"dmft2dn.def")
        atomic_json(work / "prepare.json", cfg)
        commands = [
            [
                sys.executable,
                str(Path(__file__).with_name("lapw_response_prepare.py")),
                str(work / "prepare.json"),
            ],
            [
                sys.executable,
                str(Path(__file__).with_name("lapw_total_stress_driver.py")),
                "--snapshot",
                str(work / "central"),
                "--case",
                case,
                "--potential-core",
                str(work / "central/potential_core.npz"),
                "--library",
                solver["library"],
                "--output",
                str(work / "total_stress"),
            ],
        ]
        commands[1] += ["--cache-mb", str(solver.get("cache_mb", 256.0))]
        commands[1] += ["--basis-cache-mb", str(solver.get("basis_cache_mb", 128.0))]
        if solver.get("mpi"):
            commands[1] = solver["mpi"] + commands[1] + ["--mpi"]
        for i, command in enumerate(commands):
            print(
                "RESPONSE ..... running:", shlex.join(command), file=output, flush=True
            )
            log = work / ("prepare.log" if i == 0 else "response.log")
            with log.open("w") as stream:
                # Each response rank owns k points. Avoid oversubscribed BLAS
                # inside thousands of small angular/band products. Production
                # dmft2 retains the caller's MPI/OpenMP launch unchanged.
                child_env = os.environ | {
                    "OMP_NUM_THREADS": "1",
                    "OPENBLAS_NUM_THREADS": "1",
                    "VECLIB_MAXIMUM_THREADS": "1",
                    "MKL_NUM_THREADS": "1",
                }
                result = subprocess.run(
                    command,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    env=child_env,
                )
            if result.returncode:
                raise RuntimeError(f"analytic response failed; see {log}")
        path = work / "total_stress/report.json"
        report = json.loads(path.read_text())
        if (
            report.get("converged") is not True
            or report.get("full_tensor_within_symmetry") is not True
            or report.get("diagnostic_only", False)
        ):
            raise RuntimeError(
                "analytic response did not supply a complete supported tensor"
            )
        # Only physical parent inputs are checked at SCF assembly. The isolated
        # preparation creates its own DEF files and updates its own EF.dat.
        parent = {
            k: v
            for k, v in fingerprints.items()
            if not k.startswith("parent/") and k not in {"lapw1.def", "dmft2.def", "dmft2dn.def"}
        }
        return {
            "backend": "analytic",
            "reference_id": ref,
            "workspace": str(work),
            "report": str(path),
            "report_sha256": sha256(path),
            "parent_fingerprints": parent,
            "converged": True,
            "stress_ready": True,
            "stress_applied": False,
            "scope": "complete stress within imposed native charge symmetry; pending SCF assembly",
        }
    cfg = solver | {
        "case": case,
        "reference": str(snapshot),
        "workspace": str(work / "screening"),
        "lapw1_definition": str(snapshot / "lapw1.def"),
        "dmft2_definition": str(snapshot / "dmft2.def"),
    }
    atomic_json(work / "config.json", cfg)
    command = [
        sys.executable,
        str(Path(__file__).with_name("wien_frozen_sigma_response.py")),
        str(work / "config.json"),
    ]
    print("RESPONSE ..... running:", shlex.join(command), file=output, flush=True)
    with (work / "response.log").open("w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    report_path = work / "screening/report.json"
    report = json.loads(report_path.read_text()) if report_path.exists() else {}
    if result.returncode or report.get("converged") is not True:
        raise RuntimeError(
            f'response failed or did not converge; see {work / "response.log"}'
        )
    outputs = {"density": work / "screening/response.density"}
    for side in ["plus", "minus"]:
        for suffix in ["struct", "vsp", "vns"]:
            outputs[side + "_" + suffix] = (
                work / f"screening/potential_{side}" / (case + "." + suffix)
            )
    for p in outputs.values():
        if not p.is_file() or not p.stat().st_size:
            raise RuntimeError(f"missing response output: {p}")
    return {
        "reference_id": ref,
        "workspace": str(work),
        "report": str(report_path),
        "outputs": {k: str(v) for k, v in outputs.items()},
        "output_sha256": {k: sha256(v) for k, v in outputs.items()},
        "converged": True,
        "stress_applied": False,
        "scope": "screened_density_and_potential; no stress tensor modification",
    }


def execute(
    case,
    root,
    mpi="",
    definition="dmft2.def",
    *,
    phase="standalone",
    config_path=None,
    disabled=False,
    output=None,
):
    """Shared execution path. phase='defer' NEVER runs response.

    The final run_dmft measurement calls phase='final' exactly once, after
    termination is known. This avoids evaluating response on every inner
    charge or impurity iteration, including when finish changes on the fly.
    """
    if phase not in {"standalone", "defer", "final"}:
        raise ValueError("invalid response phase")
    output = output or sys.stdout
    req = request(case, definition, config_path, disabled)
    is_charge = charge_mode(req["files"][2]) == "c"
    enabled = req["enabled"] and is_charge
    due = enabled and phase != "defer"
    run_id = time.strftime("%Y%m%dT%H%M%S") + "-" + uuid.uuid4().hex[:12]
    status_path = Path(case + ".stress_response.status.json")
    status = {
        "schema": 1,
        "run_id": run_id,
        "phase": phase,
        "requested_mode": req["mode"],
        "config": req["config_path"],
        "stress_applied": False,
        "status": "running_dmft2",
    }
    # Even disabled/deferred calls invalidate a previous published result.
    tracking = req["enabled"] or status_path.exists() or req["mode"] in ALIASES
    lock = None
    locked = False
    try:
        if tracking:
            lock = open("." + case + ".dmft2-response.lock", "a+")
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            locked = True
        if due:
            preflight(case, req, root, mpi)
        if tracking:
            atomic_json(status_path, status)
        magnetic = due and Path(case+".indmfldn").is_file() and Path(case+".indmfldn").stat().st_size > 0
        if magnetic:
            from lapw_response_policy import detect_policy
            from dmft2_spin import run_channels, run_fermif, average_charge
            policy = detect_policy(".",case,control=req["files"][2],
                                   override=req["config"].get("solver",{}).get("physical_policy"))
            run_channels([definition,"dmft2dn.def"],policy,
                lambda label,df: run_core(root,mpi,df,output),
                lambda bands: run_fermif(Path(root)/"fermif",bands,output))
            average_charge(case,root,output)
            status["physical_policy"] = policy.record()
        else:
            with execution_definition(definition) as native:
                run_core(root, mpi, native, output)
        if due:
            status["status"] = "running_response"
            atomic_json(status_path, status)
            status.update(compute_response(case, req, root, mpi, run_id, output))
            status["status"] = "complete"
        else:
            status["status"] = (
                "deferred_until_final"
                if enabled
                else "disabled" if is_charge else "skipped_noncharge"
            )
        if tracking:
            atomic_json(status_path, status)
        return status
    except BaseException as exc:
        # Failure must stop the caller before lcore/mixer; do not leave an old
        # successful sidecar looking like a correction to this failed step.
        if tracking and locked:
            status.update(status="failed", error=str(exc))
            atomic_json(status_path, status)
        raise
    finally:
        if lock is not None:
            lock.close()
