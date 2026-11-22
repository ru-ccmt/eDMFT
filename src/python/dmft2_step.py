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

ALIASES = {"STRT": "STR", "STRFT": "STRF", "FORT": "FOR"}
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
    "disk_matching",
    "dmft2",
    "wien_root",
    "force_atoms",
    "force_steps_Bohr",
    "force_source",
    "soc_operator_force",
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
def execution_definition(definition, final_force=False):
    """Translate response aliases for native executables without editing IN2.

    Native WIEN/dmft2 readers have four-character modes. STRFT and FORT are
    wrapper requests, not new Fortran input formats. Temporary definitions
    are in the original working directory so other relative paths retain
    their meaning, and are removed even on subprocess failure.
    """
    text = Path(definition).read_text()
    source = Path(entries(text)[5])
    mode = input_mode(source)
    native_mode = ALIASES.get(mode, mode)
    if final_force and mode == "TOT":
        native_mode = "FOR"
    if native_mode == mode:
        yield str(definition)
        return
    with tempfile.TemporaryDirectory(prefix=".dmft2-response-input-", dir=".") as temp:
        folder = Path(temp)
        native = folder / source.name
        lines = source.read_text().splitlines(True)
        lines[0] = re.sub(r"^\s*\S+", f"{native_mode:4s}", lines[0], count=1)
        native.write_text("".join(lines))
        df = folder / "dmft2.def"
        df.write_text(replace_filenames(text, {5: native}))
        yield str(df)


def run_core(root, mpi, definition, output, *, extra_env=None):
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
        env=None if extra_env is None else os.environ | extra_env,
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


def soc_operator_scope(case, files):
    """Eligibility for the explicit second-variational operator in a native FOR.

    Read the actual DEF and scalar output. A leftover case.inso alone does
    not identify a spinor calculation. Unsupported workflows retain their
    native FOR unless the operator is explicitly required in the config.
    """
    if any(Path(files.get(unit, "")).name.removesuffix("_x") != case + suffix
           for unit, suffix in ((9, ".vectorso"), (10, ".vectorsodn"))):
        return "the native definition does not use both spinor components"
    for suffix in (".clmup", ".clmdn", ".incup", ".incdn", ".indmfldn"):
        path = Path(case + suffix)
        if path.is_file() and path.stat().st_size:
            return "spin-polarized KS or separate DMFT spin channels"
    from lapw_response_policy import control_fields

    if int(control_fields(files[2])[2][2]) not in (0, 1):
        return "the native chemical-potential mode is outside 0 or 1"
    with Path(files[7]).open() as inp:
        first, second = inp.readline().split(), inp.readline().split()
    if len(first) < 4 or int(first[2]) != 1 or abs(int(first[3])) != 5 or not int(second[0]):
        return "the spinor derivative requires normalized projector 5 on Matsubara frequencies"
    for name in ("lapw1.def", "lapwso.def"):
        if not Path(name).is_file():
            return "missing central " + name
    so = entries(Path("lapwso.def").read_text())
    if 12 in so:
        return "LAPWSO orbital potential"
    # The normalized spinor projector currently represents an unshifted
    # atomic shell in each impurity block, with fixed local axes.
    from indmffile import Indmfl
    inl = Indmfl(case)
    inl.read()
    if any(len(shells) != 1 for shells in inl.cix.values()):
        return "several correlated shells in a single impurity block"
    if any(shift is not None or len(axes) not in (0, 3)
           for _, axes, shift, _ in inl.atoms.values()):
        return "shifted shell or incomplete local axes"
    projector = Path("projectorw.dat")
    if not projector.is_file():
        return "missing projectorw.dat"
    with projector.open() as fh:
        nsets = int(fh.readline().split()[1])
        for _ in range(nsets):
            n, inside = map(int, fh.readline().split()[1:3])
            if n != inside:
                return "radial projector extends outside its muffin tin"
            for _ in range(n):
                next(fh)
    inso = Path(so[5]).read_text().splitlines()
    if len(inso) > 4 and inso[4].split() and int(inso[4].split()[0]):
        return "relativistic p1/2 local orbitals"
    scalar = entries(Path("lapw1.def").read_text())
    output = Path(scalar[6])
    candidates = [output]
    if output.name.endswith("_x"):
        candidates += [output.with_name(output.name[:-1] + "1")]
    for name in candidates:
        if name.is_file():
            with name.open(errors="replace") as stream:
                if stream.readline().strip().lower() in {"complex lapw1c version", "real lapw1 version"}:
                    return None
    return "the central scalar output does not identify dftKS"


def preflight(case, req, root, mpi):
    """Reject unsupported physics before starting an expensive response."""
    import numpy as np
    config = req["config"]
    solver = config.get("solver", {})
    if not isinstance(solver, dict) or set(solver) - SOLVER_KEYS:
        raise ValueError(
            "unsupported solver settings; see dmft2_response_integration.txt"
        )
    backend = solver.get("backend", "analytic")
    soc_forward = backend == "native_soc_forward"
    if soc_forward:
        if req["mode"] != "FORT" or not Path(case + ".inso").is_file() or not Path("lapwso.def").is_file():
            raise ValueError("native_soc_forward requires FORT, case.inso and lapwso.def")
        supported = {"backend", "force_atoms", "force_steps_Bohr", "tolerance", "max_iterations",
                     "wien_root", "native_mpi", "physical_policy"}
        if set(solver) - supported:
            raise ValueError("unsupported native SOC solver setting: "+str(sorted(set(solver)-supported)))
    if backend not in {"native", "analytic", "native_soc_forward"}:
        raise ValueError("response backend must be native, analytic or native_soc_forward")
    if solver.get("charge_integration", "native") not in {"native", "legacy_spline"}:
        raise ValueError("charge_integration must be native or legacy_spline")
    defs = req["files"]
    so_input = Path(case + ".inso")
    soc_analytic = (backend == "analytic" and req["mode"] in {"FOR", "FORT"}
                    and so_input.is_file() and so_input.stat().st_size > 0)
    spinor = soc_forward or soc_analytic
    if spinor and not Path("lapwso.def").is_file():
        raise ValueError("spinor FORT needs the LAPWSO definition used for the central spinors")
    if soc_analytic and solver.get("force_source", "analytic") != "analytic":
        raise ValueError("select native_soc_forward explicitly for displaced spinor diagnostics")
    for suffix in ([] if spinor else [".inso"]) + [".clmup", ".clmdn", ".incup", ".incdn"]:
        p = Path(case + suffix)
        if p.is_file() and p.stat().st_size:
            if suffix == ".inso":
                raise ValueError(f"SOC is unsupported for {req['mode']} with backend {backend}")
            raise ValueError(f"response requires an unpolarized KS/core potential; found {p.name}")
    control = [
        line.split("#", 1)[0].split()
        for line in Path(defs[2]).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ][0]
    if len(control) < 3 or control[0] != "c" or int(control[2]) not in (0, 1):
        raise ValueError(
            "response requires charge mode c and recomputeEF=0 or 1"
        )
    if req["mode"] == "FOR" and not soc_analytic:
        raise ValueError("a direct SOC force requires the analytic spinor backend")
    if req["mode"] not in {"FOR", "STR", "STRF", "STRT", "STRFT", "FORT"}:
        raise ValueError(
            "response requires STR, STRF, STRT, STRFT or FORT in the actual IN2 input"
        )
    force_requested = req["mode"] in {"FOR", "FORT", "STRFT"}
    if force_requested and backend != "analytic" and not soc_forward:
        raise ValueError("screened atomic forces require the analytic charge-adjoint backend")
    force_source = solver.get("force_source", "analytic")
    if force_source not in {"analytic", "displaced"}:
        raise ValueError("solver.force_source must be analytic or displaced")
    force_atoms = solver.get("force_atoms")
    if force_atoms is not None and (
        not isinstance(force_atoms, list) or not force_atoms
        or any(type(atom) is not int or atom < 1 for atom in force_atoms)
    ):
        raise ValueError("solver.force_atoms must be a nonempty list of 1-based atom indices")
    force_steps = solver.get("force_steps_Bohr", [0.01, 0.005])
    if force_requested and (
        not isinstance(force_steps, list) or len(force_steps) not in (1, 2)
        or any(type(h) not in (int, float) or not np.isfinite(h) or h <= 0 for h in force_steps)
        or len(set(force_steps)) != len(force_steps)
    ):
        raise ValueError("solver.force_steps_Bohr needs one or two distinct positive steps")
    if soc_forward and len(force_steps) != 2:
        raise ValueError("native_soc_forward needs two step sizes for extrapolation")
    if not Path("lapw1.def").is_file():
        raise ValueError("response requires a prepared lapw1.def")
    from wstruct import Struct

    s = Struct()
    s.ReadStruct(case + ".struct")
    from lapw_lattice import NativeLattice
    from lapw_response_policy import detect_policy
    NativeLattice.from_structure(s)
    policy = detect_policy(".", case, control=defs[2], override=solver.get("physical_policy"), native_soc=spinor)
    if req["mode"] in {"FOR", "FORT"} and policy.magnetic:
        raise NotImplementedError("FORT magnetic atomic sources have not been validated")
    all_forces_fixed = False
    if force_requested and not policy.magnetic:
        from force_response_modes import internal_modes, check_force_frame
        modes = internal_modes(s, force_atoms)
        all_forces_fixed = not modes
        check_force_frame(s, force_atoms or range(1, s.nat + 1))
    if soc_analytic and req["mode"] == "FORT" and not all_forces_fixed:
        reason = soc_operator_scope(case, defs)
        if reason:
            raise ValueError("unsupported analytic spinor FORT input: " + reason)
    if soc_forward:
        if len(force_steps) != 2 or not 0 < solver.get("tolerance", 1e-3) < 1:
            raise ValueError("SOC forward response requires two steps and a relative tolerance between zero and one")
        if type(solver.get("max_iterations", 50)) is not int or solver.get("max_iterations",50)<1:
            raise ValueError("SOC forward response requires positive integer max_iterations")
    if spinor:
        if Path(defs[9]).name.removesuffix("_x") != case+".vectorso" or Path(defs[10]).name.removesuffix("_x") != case+".vectorsodn":
            raise ValueError("SOC forward response needs both spinor components")
        scalar_control = entries(Path("lapw1.def").read_text())[5]
        scalar_worker = "lapw1c" if scalar_control.endswith(".in1c") else "lapw1"
        for name in [scalar_worker, "lapwso", "dmft2"]:
            if not os.access(Path(root)/name,os.X_OK):
                raise ValueError("missing native SOC worker: "+name)
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
        Path(solver.get("dmft2", str(Path(root) / ("dmft2" if spinor else "dmft2_response"))))
        .expanduser()
        .resolve()
    )
    # Inherit the ordinary dmft2 MPI prefix automatically. Python's mpi4py
    # and the Fortran workers may link different MPI implementations, so the
    # shared resolver verifies/adapts the Python launcher without changing
    # the ordinary dmft2 command. Explicit JSON argv lists still override it.
    from response_mpi import resolve_launches
    launch, native_launch = resolve_launches(mpi, solver, analytic=backend == "analytic")
    direct_soc = req["mode"] == "FOR" and soc_analytic
    for p in ([exe] if direct_soc else [wi / "x", wi / "lapw1", wi / "lapw1c", exe]):
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
        if type(solver.get("disk_matching", False)) is not bool:
            raise ValueError("solver.disk_matching must be true or false")
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
        if req["mode"] not in {"FOR", "FORT"} and (
            not Path(case + ".in0").is_file() or input_mode(case + ".in0") != "STR"
        ):
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
        if direct_soc:
            # The native force already supplied its potential and electronic
            # moments. Neither potential nor core response is used here.
            del names["lapw0_worker"], names["lcore_worker"]
        if soc_analytic:
            names["lapw1_worker"] = scalar_worker
            names["lapwso_worker"] = "lapwso"
        solver = dict(solver)
        if force_requested:
            solver["force_source"] = force_source
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
        "all_forces_fixed": all_forces_fixed,
        "spinor": spinor,
        "wien_root": str(wi),
        "dmft2": str(exe),
        "mpi": launch,
        "native_mpi": native_launch,
        "charge_integration": solver.get("charge_integration", "native"),
    }


def compute_response(case, req, root, mpi, run_id, output, *, captured_soc=None):
    """Snapshot the final physical inputs, solve once, and publish a sidecar.

    No cached correction is applied to another charge iterate. The status
    file names this exact snapshot. Solver children invoke executables
    directly, avoiding recursive x_dmft/run_dmft response execution.
    """
    solver = preflight(case, req, root, mpi)
    from stress_response_summary import progress, summary_path
    journal = summary_path(case).resolve()
    if req["mode"] == "STRFT" and solver["physical_policy"]["magnetic"]:
        warning = "Magnetic STRFT: screened stress is available; screened atomic force is not yet validated"
        print("RESPONSE ..... WARNING:", warning, file=output, flush=True)
        progress(journal, "WARNING: " + warning)
    if solver.get("backend") == "analytic":
        for label, key in [("Python", "mpi"), ("native", "native_mpi")]:
            command = solver.get(key, [])
            print("RESPONSE ..... " + label + " MPI: " +
                  (shlex.join(command) if command else "serial"),
                  file=output, flush=True)
            progress(journal, label + " MPI: " + (shlex.join(command) if command else "serial"))
    base = (
        Path(req["config"].get("workspace_root", f"{case}.response_runs"))
        .expanduser()
        .resolve()
    )
    base.mkdir(parents=True, exist_ok=True)
    work = base / run_id
    work.mkdir()
    progress(journal, "Saving isolated central inputs in " + str(work))
    snapshot = work / "reference"
    snapshot.mkdir()
    fingerprints = {}
    direct_soc = req["mode"] == "FOR"
    soc_forward = solver["backend"] == "native_soc_forward"
    spinor = solver.get("spinor", soc_forward)
    suffixes = SNAPSHOT_SUFFIXES | ({".inso", ".in0abp", ".scfso"} if spinor else set())
    for p in sorted(Path(".").iterdir()):
        if p.is_file() and (
            (p.name.startswith(case + ".") and p.suffix in suffixes)
            or p.name in {"EF.dat", "projectorw.dat", "Eorb_imp.dat", "params.dat"}
            or (spinor and p.name == "SOlapm.dat")
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
    ] + ([("lapwso.def", Path("lapwso.def"))] if spinor else []) + ([("dmft2dn.def", Path("dmft2dn.def"))] if solver["physical_policy"]["magnetic"] else []):
        (snapshot / name).write_text(local_definition(source.read_text()))
        fingerprints[name] = sha256(snapshot / name)
    if spinor and not solver["all_forces_fixed"]:
        from lapw_soc_rank import parent_ranks
        from wstruct import Struct
        struct = Struct()
        struct.ReadStruct(case + ".struct")
        scalar = entries(Path("lapw1.def").read_text())[10]
        manifest = parent_ranks(scalar, req["files"][9], req["files"][10], struct.nat)
        atomic_json(snapshot / "spinor_ranks.json", manifest)
        fingerprints["spinor_ranks.json"] = sha256(snapshot / "spinor_ranks.json")
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
    for name in ([] if direct_soc else ["lapw0", "lapw1", "lapw1c", "lcore"]):
        p = Path(solver["wien_root"]) / name
        if p.is_file():
            identities[str(p)] = sha256(p)
    if solver.get("backend", "native") == "analytic":
        for key in ("lapw0_worker", "lapw1_worker", "lcore_worker", "lapwso_worker", "library"):
            if key in solver:
                identities[solver[key]] = sha256(Path(solver[key]))
    if spinor:
        for name in ("lapwso","lapw1c" if req["files"][3].endswith("c") else "lapw1"):
            executable = Path(root)/name
            identities[str(executable)] = sha256(executable)
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
        "stress_response_summary",
        "dmft2_step",
        "wien_frozen_sigma_response",
        "frozen_sigma_screening",
        "wien_response_density",
        "lapw_charge_response",
        "lapw_hartree_periodic",
        "lapw_hartree_radial",
        "lapw_density_response",
        "stress_response_adjoint",
        "force_response_modes",
        "lapw_atomic_force",
        "lapw_response_scratch",
        "lapw_augmentation_response",
        "lapw_compact_response",
        "lapw_complement_response",
        "lapw_energy_readout",
        "lapw_nonspherical_response",
        "lapw_physical_potential",
        "lapw_potential_response",
        "lapw_spherical_response",
        "lapw_strain_field",
        "lapw_strain_potential",
        "force_screened_diagnostic",
        "lapw_soc_force_reference",
        "lapw_soc_subspace", "lapw_soc_projector", "lapw_soc_response",
        "lapw_soc_valence", "lapw_soc_electronic", "lapw_soc_screened", "lapw_soc_rank",
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
    if soc_forward:
        from lapw_soc_force_reference import run_reference
        progress(journal,"Solving the native spinor charge response for each requested atomic displacement")
        force = run_reference(snapshot,case,root,solver["wien_root"],work/"soc_forward",
            atoms=solver.get("force_atoms"),steps=solver.get("force_steps_Bohr",[0.01,0.005]),
            tolerance=solver.get("tolerance",1e-3),max_iterations=solver.get("max_iterations",50),
            mpi=solver["native_mpi"],progress=lambda line: progress(journal,line))
        force_path = work/"total_force"/"report.json"
        force_path.parent.mkdir()
        atomic_json(force_path,force)
        digest=sha256(force_path)
        parent={k:v for k,v in fingerprints.items() if not k.startswith("parent/") and not k.endswith(".def")}
        return dict(backend=solver["backend"],reference_id=ref,workspace=str(work),
                    converged=True,stress_ready=False,force_ready=True,force_applied=False,
                    report=str(force_path),report_sha256=digest,force_report=str(force_path),
                    force_report_sha256=digest,parent_fingerprints=parent,
                    scope="screened frozen-Sigma spinor forces from the native forward charge response")
    if solver.get("backend", "native") == "analytic":
        # A new ordinary FOR is possible on every electronic/ionic step.
        # Use the parent measurement's temporary area for its large k files.
        # Retain the small immutable reference, timings and force reports.
        central_dir = Path(captured_soc).resolve() / "central" if captured_soc else work / "central"
        cfg = {
            "case": case,
            "physical_policy": solver["physical_policy"],
            "fermif": str(Path(root)/"fermif"),
            "reference": str(snapshot),
            "workspace": str(central_dir),
            "wien_root": solver["wien_root"],
            "dmft2": solver["dmft2"],
            "native_mpi": solver.get("native_mpi", []),
            # Each Python/mpi4py k-pool rank starts a serial LAPW1 worker.
            # Its launcher must match mpi4py, including on mixed-MPI systems.
            "lapw1_mpi": solver.get("mpi", []),
            "lapw1_definition": str(snapshot / "lapw1.def"),
            "dmft2_definition": str(snapshot / "dmft2.def"),
            **{k: solver[k] for k in ("lapw0_worker", "lapw1_worker", "lcore_worker") if k in solver},
        }
        if solver["physical_policy"]["magnetic"]:
            cfg["dmft2dn_definition"] = str(snapshot/"dmft2dn.def")
        if spinor:
            if captured_soc is not None:
                cfg["captured_soc"] = str(Path(captured_soc).resolve())
            cfg.update(lapwso_worker=solver["lapwso_worker"], lapwso_definition=str(snapshot / "lapwso.def"),
                       preserve_soc_normalization=direct_soc)
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
                str(central_dir),
                "--case",
                case,
                "--potential-core",
                str(central_dir / "potential_core.npz"),
                "--library",
                solver["library"],
                "--output",
                str(work / "total_stress"),
                "--progress-file",
                str(journal),
            ],
        ]
        commands[1] += ["--cache-mb", str(solver.get("cache_mb", 256.0))]
        commands[1] += ["--basis-cache-mb", str(solver.get("basis_cache_mb", 128.0))]
        if solver.get("disk_matching", False):
            commands[1].append("--disk-matching")
        if req["mode"] == "FORT":
            commands[1].append("--adjoint-only")
        force_ready = req["mode"] in {"FOR", "FORT", "STRFT"} and not solver["physical_policy"]["magnetic"]
        fixed_by_symmetry = force_ready and solver.get("all_forces_fixed", False)
        analytic_forces = force_ready and not fixed_by_symmetry and solver["force_source"] == "analytic"
        if analytic_forces:
            commands[1].append("--atomic-forces")
            if solver.get("force_atoms"):
                commands[1] += ["--force-atoms", *map(str, solver["force_atoms"])]
        if spinor:
            commands[1] = [sys.executable, str(Path(__file__).with_name("lapw_soc_screened.py")),
                           "--snapshot", str(central_dir), "--case", case,
                           "--library", solver["library"], "--output", str(work / "total_stress"),
                           "--progress-file", str(journal),
                           "--cache-mb", str(solver.get("cache_mb", 256.0))]
            if solver.get("force_atoms"):
                commands[1] += ["--force-atoms", *map(str, solver["force_atoms"])]
            if direct_soc:
                commands[1].append("--soc-operator-only")
        if solver.get("mpi"):
            if spinor:
                # Abort peers when one rank discovers an incompatible local
                # record before the others reach the next collective.
                commands[1][1:1] = ["-m", "mpi4py"]
            commands[1] = solver["mpi"] + commands[1] + ["--mpi"]
        if fixed_by_symmetry and req["mode"] in {"FORT", "FOR"}:
            # No invariant first derivative exists for any selected atom.
            # Preserve the same fingerprinted publication path as ordinary forces.
            commands = []
            path = work / "total_stress/report.json"
            path.parent.mkdir()
            atomic_json(path, {"symmetry_fixed": True, "converged": True})
            progress(journal, "All selected forces vanish by site symmetry. No adjoint is needed.")
        for i, command in enumerate(commands):
            if i == 0:
                message = ("Preparing the central scalar matrices and SOC radial integrals" if direct_soc
                           else "Preparing central LAPW0/LAPW1/DMFT2/LCORE exports")
            elif direct_soc:
                message = "Computing the explicit frozen-radial SOC force"
            elif req["mode"] == "FORT":
                message = "Starting shared atomic-force adjoint"
            else:
                message = "Starting analytic screened response and stress sources"
            progress(journal, message)
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
            progress(journal, ("Central preparation" if i == 0 else "Analytic response") + " finished")
        path = work / "total_stress/report.json"
        report = json.loads(path.read_text())
        if captured_soc is not None and commands:
            # Mixer needs the site rotations after the temporary matrices are gone.
            report["structure"] = str(snapshot / (case + ".struct"))
            atomic_json(path, report)
            shutil.copy2(central_dir / "preparation.json", work / "preparation.json")
        if (
            report.get("converged") is not True
            or (commands and not direct_soc and report.get("adjoint", {}).get("converged") is not True)
            or (req["mode"] not in {"FOR", "FORT"} and report.get("full_tensor_within_symmetry") is not True)
            or report.get("diagnostic_only", False)
        ):
            raise RuntimeError(
                "analytic response did not supply the requested converged adjoint/tensor"
            )
        force_path = None
        if force_ready:
            if fixed_by_symmetry:
                from wstruct import Struct
                s = Struct(inAngs=False)
                structure = snapshot / (case + ".struct")
                s.ReadStruct(str(structure))
                result = dict(case=case, structure=str(structure), components=[],
                              selected_atoms=solver.get("force_atoms") or list(range(1, s.nat+1)),
                              force_mRy_Bohr=[[0., 0., 0.] for _ in range(s.nat)],
                              source_method="site_symmetry", scope="all selected sites fixed by symmetry")
                if direct_soc:
                    result["application"] = "add_soc_operator"
            elif analytic_forces:
                result = report if spinor else report["atomic_forces"]
            else:
                from force_response_modes import run_force_modes
                progress(journal, "Preparing symmetry-preserving atomic force sources with the shared adjoint")
                result = run_force_modes(
                    work, work / "force_sources",
                    steps=solver.get("force_steps_Bohr", [0.01, 0.005]),
                    atoms=solver.get("force_atoms"),
                    use_mpi=True,
                    progress=lambda message: progress(journal, message),
                )
            force_path = work / "total_force" / "report.json"
            force_path.parent.mkdir()
            atomic_json(force_path, result)
            progress(journal, ("Explicit SOC force additions: " if direct_soc else "Screened atomic forces: ") + json.dumps(result["force_mRy_Bohr"]))
        # Only physical parent inputs are checked at SCF assembly. The isolated
        # preparation creates its own DEF files and updates its own EF.dat.
        parent = {
            k: v
            for k, v in fingerprints.items()
            if not k.startswith("parent/") and k not in {"spinor_ranks.json", "lapw1.def", "lapwso.def", "dmft2.def", "dmft2dn.def"}
        }
        return {
            "backend": "analytic",
            "reference_id": ref,
            "workspace": str(work),
            "report": str(path),
            "report_sha256": sha256(path),
            "parent_fingerprints": parent,
            "converged": True,
            "stress_ready": req["mode"] not in {"FOR", "FORT"},
            "force_ready": force_ready,
            "force_report": None if force_path is None else str(force_path),
            "force_report_sha256": None if force_path is None else sha256(force_path),
            "stress_applied": False,
            "scope": "requested symmetry-preserving frozen-Sigma response; pending SCF assembly",
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
        "force_applied": False,
        "scope": "screened_density_and_potential; no stress tensor modification",
    }


def execute(
    case,
    root,
    mpi="",
    definition="dmft2.def",
    *,
    phase="standalone",
    final_force=False,
    config_path=None,
    disabled=False,
    output=None,
):
    """Run the native measurement and any enabled force or stress terms.

    The final run_dmft measurement calls phase='final' exactly once, after
    termination is known. This avoids evaluating response on every inner
    charge or impurity iteration, including when finish changes on the fly.
    """
    if phase not in {"standalone", "defer", "final"}:
        raise ValueError("invalid response phase")
    output = output or sys.stdout
    req = request(case, definition, config_path, disabled)
    is_charge = charge_mode(req["files"][2]) == "c"
    # FOR includes the explicit SOC operator on its actual charge pass.
    # The screened charge adjoint remains final-only. A private TOT -> FOR
    # measurement for phonons follows exactly the same operator selection.
    solver_control = req["config"].get("solver", {})
    option = solver_control.get("soc_operator_force")
    if "soc_operator_force" in solver_control and type(option) is not bool:
        raise ValueError("solver.soc_operator_force must be true or false")
    native_force = req["mode"] == "FOR" or (req["mode"] == "TOT" and final_force and phase == "final")
    direct_soc = False
    if (is_charge and native_force and option is not False and not disabled
            and req["config"].get("enabled", True)
            and Path(req["files"].get(9, "")).name.removesuffix("_x") == case + ".vectorso"
            and solver_control.get("backend", "analytic") == "analytic"):
        reason = soc_operator_scope(case, req["files"])
        if reason is None:
            direct_soc = True
        elif option is True:
            raise ValueError("explicit SOC FOR requested: " + reason)
        else:
            print("FOR ..... native SOC force only: " + reason, file=output, flush=True)
    if direct_soc:
        req = dict(req, mode="FOR", enabled=True)
    enabled = req["enabled"] and is_charge
    due = enabled and (direct_soc or phase != "defer")
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
    capture = None
    from stress_response_summary import begin, finish, progress, summary_path
    journal = summary_path(case)
    try:
        if tracking:
            lock = open("." + case + ".dmft2-response.lock", "a+")
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            locked = True
            begin(case, status)
        if due:
            settings = preflight(case, req, root, mpi)
            if direct_soc and not settings["all_forces_fixed"]:
                capture = tempfile.TemporaryDirectory(dir=Path.cwd(), prefix=".dmft2-soc-force-")
        if tracking:
            atomic_json(status_path, status)
            progress(journal, "Running ordinary dmft2" +
                     ("; response follows this step" if due else "; response not due on this step"))
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
            # A phonon measurement uses a private FOR input on the final
            # charge pass. The next run can start from the same TOT input.
            measure_force = final_force and is_charge and phase == "final"
            with execution_definition(definition, final_force=measure_force) as native:
                if capture is None:
                    run_core(root, mpi, native, output)
                else:
                    run_core(root, mpi, native, output, extra_env={
                        "EDMFT_RESPONSE_EXPORT": str(Path(capture.name)/"response"),
                        "EDMFT_SOC_OPERATOR_READOUT": "1"})
        if due:
            status["status"] = "running_response"
            atomic_json(status_path, status)
            capture_args = {"captured_soc": capture.name} if capture is not None else {}
            status.update(compute_response(case, req, root, mpi, run_id, output, **capture_args))
            status["status"] = "complete"
        else:
            status["status"] = (
                "deferred_until_final"
                if enabled
                else "disabled" if is_charge else "skipped_noncharge"
            )
        if tracking:
            atomic_json(status_path, status)
            finish(case, status)
        return status
    except BaseException as exc:
        # Failure must stop the caller before lcore/mixer; do not leave an old
        # successful sidecar looking like a correction to this failed step.
        if tracking and locked:
            status.update(status="failed", error=str(exc))
            atomic_json(status_path, status)
            finish(case, status)
        raise
    finally:
        if capture is not None:
            capture.cleanup()
        if lock is not None:
            lock.close()
