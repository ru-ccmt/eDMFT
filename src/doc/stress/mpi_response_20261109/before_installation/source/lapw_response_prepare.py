#!/usr/bin/env python3
"""Export one isolated central state for the analytic stress driver.

This is a single central LAPW0 -> LAPW1 -> DMFT2 -> LCORE pass, not a strained
SCF loop. DMFT2 first determines the central chemical potential, then its
fixed-mu export pass records that same state. The analytic derivative later
enforces fixed N. Parent files are copied, never linked or modified.

Workers are explicit configuration paths: lapw0_worker/lapw1_worker/lcore_worker
must support the native response exports; dmft2 is the response-enabled worker.
Stock executables without these exports fail validation. The native workers
are separate from ordinary production executables, and all their scratch
outputs live under the new workspace. MPI preparation is not supported yet.
"""
from pathlib import Path
import argparse, json, os, re, shutil, subprocess, time
from dmft2_step import local_definition, entries
from lapw_native_potential_export import prepare_native_potential_core
from lapw_potential_core_state import fingerprint_files


def prepare(config):
    required = {
        "case",
        "reference",
        "workspace",
        "wien_root",
        "lapw0_worker",
        "lapw1_worker",
        "lcore_worker",
        "dmft2",
        "lapw1_definition",
        "dmft2_definition",
    }
    if set(config) != required:
        raise ValueError("preparation config keys must be " + str(sorted(required)))
    case = config["case"]
    reference = Path(config["reference"]).resolve()
    work = Path(config["workspace"]).resolve()
    if not re.fullmatch(r"[A-Za-z0-9_.+-]+", case):
        raise ValueError("invalid case name")
    work.mkdir(parents=True, exist_ok=False)
    env = {
        k: v
        for k, v in os.environ.items()
        if not k.startswith(("EDMFT_", "DFTKS_", "WIEN_"))
    }
    env.update(
        WIENROOT=str(Path(config["wien_root"]).resolve()),
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
    )
    suffixes = [
        "struct",
        "in0",
        "in1",
        "in1c",
        "in2",
        "in2c",
        "inc",
        "inm",
        "klist",
        "kgen",
        "rotlm",
        "indmfl",
        "indmfi",
        "indmf2",
        "clmsum",
    ]
    names = (
        [case + "." + s for s in suffixes]
        + ["EF.dat", "projectorw.dat", "Eorb_imp.dat"]
        + [p.name for p in reference.glob("sig.inp*")]
    )
    for suffix in ["inso", "indmfldn", "clmup", "clmdn"]:
        p = reference / (case + "." + suffix)
        if p.is_file() and p.stat().st_size:
            raise NotImplementedError("nonmagnetic, no-SOC preparation required")
    files = {name: reference / name for name in names if (reference / name).is_file()}
    before = fingerprint_files(files)
    for name, path in files.items():
        shutil.copy2(path, work / name)
    for label in ["lapw1", "dmft2"]:
        (work / (label + ".def")).write_text(
            local_definition(Path(config[label + "_definition"]).read_text())
        )
    for suffix in ["in0", "in2", "in2c"]:
        p = work / (case + "." + suffix)
        if p.is_file() and p.stat().st_size:
            lines = p.read_text().splitlines()
            lines[0] = re.sub(r"^\s*\S+", "TOT", lines[0], count=1)
            p.write_text("\n".join(lines) + "\n")
    control = work / (case + ".indmf2")
    lines = control.read_text().splitlines()
    j = next(
        i for i, l in enumerate(lines) if l.strip() and not l.lstrip().startswith("#")
    )
    words = lines[j].split()
    if words[0] != "c" or int(words[2]) != 1:
        raise NotImplementedError("charge mode with recomputeEF=1 required")
    times = []

    def run(label, exe, args, extra=None, marker=None):
        start = time.perf_counter()
        log = work / (label + ".log")
        with log.open("w") as stream:
            result = subprocess.run(
                [str(exe), *args],
                cwd=work,
                env=env | (extra or {}),
                stdout=stream,
                stderr=subprocess.STDOUT,
            )
        text = log.read_text(errors="replace")
        if result.returncode or any(
            x in text
            for x in ["Fortran runtime error", "ERROR STOP", "Error termination"]
        ):
            raise RuntimeError(str(log))
        if marker is not None and marker + " END" not in text:
            raise RuntimeError("missing completion in " + str(log))
        times.append(dict(program=label, seconds=time.perf_counter() - start))

    x = Path(config["wien_root"]) / "x"
    for program in ["lapw0", "lcore"]:
        run(program + "_definition", x, [program, "-f", case, "-d"])
    run("lapw0", config["lapw0_worker"], ["lapw0.def"])
    run(
        "lapw1",
        config["lapw1_worker"],
        ["lapw1.def"],
        {"WIEN_LAPW1_MATRIX_DUMP": str(work / "native")},
        "LAPW1",
    )
    run("dmft2_number", config["dmft2"], ["dmft2.def"], marker="DMFT2")
    words[1] = (work / "EF.dat").read_text().split()[0]
    words[2] = "0"
    lines[j] = " ".join(words)
    control.write_text("\n".join(lines) + "\n")
    run(
        "dmft2",
        config["dmft2"],
        ["dmft2.def"],
        {"EDMFT_RESPONSE_EXPORT": str(work / "response")},
        "DMFT2",
    )
    run("lcore", config["lcore_worker"], ["lcore.def"])
    prepare_native_potential_core(work, case, work / "potential_core.npz")
    if before != fingerprint_files(files):
        raise RuntimeError("parent inputs changed during preparation")
    if not (work / "response.meta").is_file() or not list(
        work.glob("native.ham.k*.bin")
    ):
        raise RuntimeError("workers did not export native response data")
    report = dict(
        case=case,
        parent_fingerprints=before,
        config=config,
        timings=times,
        scope="one central snapshot, no strained SCF loop",
        stress_applied=False,
    )
    (work / "preparation.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(prepare(json.loads(args.config.read_text())), indent=2), flush=True
    )
