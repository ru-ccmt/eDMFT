#!/usr/bin/env python3
"""Export one isolated central state for the analytic stress driver.

This is a single central LAPW0 -> LAPW1 -> DMFT2 -> LCORE pass, not a strained
SCF loop. The configured physical ensemble determines whether DMFT2/fermif first
updates the common chemical potential. Fixed-mu exports then record both spin
branches. Preparation metadata retains the original response ensemble. Parent files are copied, never linked or modified.

Workers are explicit configuration paths: lapw0_worker/lapw1_worker/lcore_worker
must support the native response exports; dmft2 is the response-enabled worker.
Stock executables without these exports fail validation. The native workers
are separate from ordinary production executables, and all their scratch
outputs live under the new workspace. native_mpi distributes the two DMFT2
passes. lapw1_mpi runs the unchanged serial WIEN worker in independent k-point
pools and merges validated streams. Only LAPW0/LCORE remain replicated setup.
"""
from pathlib import Path
from contextlib import contextmanager
from lapw_response_policy import detect_policy
from dmft2_spin import run_channels, run_fermif
import argparse, json, os, re, shutil, subprocess, time, sys
from dmft2_step import local_definition, entries
from lapw_native_potential_export import prepare_native_potential_core
from lapw_potential_core_state import fingerprint_files


@contextmanager
def working_directory(path):
    """Scoped native-file context, also available on Python 3.10."""
    previous = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(previous)


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
    if set(config) - {"native_mpi", "lapw1_mpi", "physical_policy", "dmft2dn_definition", "fermif"} != required:
        raise ValueError("preparation config keys must be " + str(sorted(required)))
    native_mpi = config.get("native_mpi", [])
    if not isinstance(native_mpi, list) or any(
        not isinstance(x, str) or any(c in x for c in "$;|><") for x in native_mpi
    ):
        raise ValueError("native_mpi must be an explicit argument list")
    lapw1_mpi = config.get("lapw1_mpi", [])
    if not isinstance(lapw1_mpi, list) or any(
        not isinstance(x, str) or any(c in x for c in "$;|><") for x in lapw1_mpi
    ):
        raise ValueError("lapw1_mpi must be an explicit argument list")
    case = config["case"]
    reference = Path(config["reference"]).resolve()
    work = Path(config["workspace"]).resolve()
    if not re.fullmatch(r"[A-Za-z0-9_.+-]+", case):
        raise ValueError("invalid case name")
    policy = detect_policy(reference, case, override=config.get("physical_policy"))
    if policy.magnetic and "dmft2dn_definition" not in config:
        raise ValueError("two-spin preparation requires dmft2dn_definition")
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
        "indmfldn",
        "indmf2dn",
        "clmsum",
    ]
    names = (
        [case + "." + s for s in suffixes]
        + ["EF.dat", "projectorw.dat", "Eorb_imp.dat", "params.dat"]
        + [p.name for p in reference.glob("sig.inp*")]
    )
    files = {name: reference / name for name in names if (reference / name).is_file()}
    before = fingerprint_files(files)
    for name, path in files.items():
        shutil.copy2(path, work / name)
    for label in ["lapw1", "dmft2"] + (["dmft2dn"] if policy.magnetic else []):
        (work / (label + ".def")).write_text(
            local_definition(Path(config[label + "_definition"]).read_text())
        )
    for suffix in ["in0", "in2", "in2c"]:
        p = work / (case + "." + suffix)
        if p.is_file() and p.stat().st_size:
            lines = p.read_text().splitlines()
            lines[0] = re.sub(r"^\s*\S+", "TOT", lines[0], count=1)
            p.write_text("\n".join(lines) + "\n")
    times = []

    def run(label, exe, args, extra=None, marker=None):
        start = time.perf_counter()
        log = work / (label + ".log")
        with log.open("w") as stream:
            result = subprocess.run(
                (native_mpi if label.startswith("dmft2") else []) + [str(exe), *args],
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
    if lapw1_mpi:
        run(
            "lapw1",
            lapw1_mpi[0],
            lapw1_mpi[1:]
            + [
                sys.executable,
                str(Path(__file__).with_name("lapw_response_kpool.py")),
                "--workspace",
                str(work),
                "--case",
                case,
                "--worker",
                config["lapw1_worker"],
            ],
            marker="LAPW1",
        )
    else:
        run(
            "lapw1",
            config["lapw1_worker"],
            ["lapw1.def"],
            {"WIEN_LAPW1_MATRIX_DUMP": str(work / "native")},
            "LAPW1",
        )
    definitions = [str(work/('dmft2'+channel+'.def')) for channel in policy.channels]
    fermif = Path(config.get("fermif", Path(config["dmft2"]).parent/"fermif")).resolve()
    def common_fermi(bands):
        with (work/"fermif.log").open("w") as stream:
            run_fermif(fermif, bands, stream)
    with working_directory(work):
        run_channels(definitions, policy,
            lambda label, definition: run(label+"_number", config["dmft2"], [definition], marker="DMFT2"),
            common_fermi, search_only=True)
        def export(label, definition):
            suffix = "_dn" if label.endswith("dn") else ""
            run(label, config["dmft2"], [definition],
                {"EDMFT_RESPONSE_EXPORT": str(work/("response"+suffix))}, "DMFT2")
        run_channels(definitions, policy, export, common_fermi, export=True)
    run("lcore", config["lcore_worker"], ["lcore.def"])
    prepare_native_potential_core(work, case, work / "potential_core.npz")
    if before != fingerprint_files(files):
        raise RuntimeError("parent inputs changed during preparation")
    if any(not (work / ("response"+("_dn" if ch else "")+".meta")).is_file() for ch in policy.channels) or not list(
        work.glob("native.ham.k*.bin")
    ):
        raise RuntimeError("workers did not export native response data")
    report = dict(
        case=case,
        physical_policy=policy.record(),
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
