#!/usr/bin/env python3
"""Execute prepared dmft2, then an optional frozen-Sigma response.

The callers own input preparation and sigma splitting. This module never
regenerates those inputs. run_dmft defers response until a single final pass;
a standalone x_dmft execution may respond immediately. The present response
backend exports density/potential tangents, NOT a corrected mixer tensor.
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
    "strain",
    "probe_step",
    "tolerance",
    "max_iterations",
    "mpi",
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
    defs = req["files"]
    # Backend initially has one unpolarized density and a fixed-N solver.
    for suffix in [".inso", ".indmfldn", ".clmup", ".clmdn"]:
        p = Path(case + suffix)
        if p.is_file() and p.stat().st_size:
            raise ValueError(
                "response currently supports nonmagnetic calculations without SOC only"
            )
    control = [
        line.split("#", 1)[0].split()
        for line in Path(defs[2]).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ][0]
    if len(control) < 3 or control[0] != "c" or int(control[2]) != 1:
        raise ValueError(
            "response currently requires charge mode c and recomputeEF=1 (fixed N)"
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
    if s.lattyp.strip() != "P" or not np.allclose([s.alpha, s.beta, s.gamma], 90.0):
        raise ValueError("response currently supports primitive orthogonal cells only")
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
    launch = solver.get("mpi", shlex.split(mpi))
    if not isinstance(launch, list) or any(not isinstance(x, str) for x in launch):
        raise ValueError("solver.mpi must be an argument list")
    if any("$" in x or any(t in x for t in [";", "|", ">", "<"]) for x in launch):
        raise ValueError(
            "set solver.mpi to explicit arguments, without shell expansions"
        )
    for p in [wi / "x", wi / "lapw1", wi / "lapw1c", exe]:
        if not p.is_file() or not os.access(p, os.X_OK):
            raise ValueError(f"missing executable: {p}")
    for name in [case + ".clmsum", "EF.dat", "projectorw.dat"]:
        if not Path(name).is_file():
            raise ValueError(f"missing response input: {name}")
    return solver | {"wien_root": str(wi), "dmft2": str(exe), "mpi": launch}


def compute_response(case, req, root, mpi, run_id, output):
    """Snapshot the final physical inputs, solve once, and publish a sidecar.

    No cached correction is applied to another charge iterate. The status
    file names this exact snapshot. Solver children invoke executables
    directly, avoiding recursive x_dmft/run_dmft response execution.
    """
    solver = preflight(case, req, root, mpi)
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
            or p.name in {"EF.dat", "projectorw.dat", "Eorb_imp.dat"}
            or re.fullmatch(r"sig\.inp\d*", p.name)
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
    ]:
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
    scripts = {}
    for module in [
        "dmft2_step",
        "wien_frozen_sigma_response",
        "frozen_sigma_screening",
        "wien_response_density",
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
