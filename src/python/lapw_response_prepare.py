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
For a spinor force, the replay regenerates the parent's per-k ranks and LO
matching convention. An ordinary FOR reuses the primary DMFT2 moments and
potential, so it needs just dftKS and a radial-only LAPWSO export.
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
    electronic_model = config.get("electronic_model", "dmft")
    if electronic_model not in {"dmft", "dft"}:
        raise ValueError("electronic_model must be dmft or dft")
    required = {
        "case",
        "reference",
        "workspace",
        "wien_root",
        "lapw0_worker",
        "lapw1_worker",
        "lcore_worker",
        "lapw1_definition",
    }
    if electronic_model == "dmft":
        required |= {"dmft2", "dmft2_definition"}
    soc = "lapwso_worker" in config or "lapwso_definition" in config
    if soc:
        if electronic_model != "dmft":
            raise ValueError("analytic SOC force requires an eDMFT snapshot")
        required |= {"lapwso_worker", "lapwso_definition"}
    allowed = {"native_mpi", "lapw1_mpi", "physical_policy", "dmft2dn_definition", "fermif", "electronic_model", "preserve_soc_normalization", "captured_soc"}
    if "captured_soc" in config:
        if not soc:
            raise ValueError("captured SOC moments require a spinor preparation")
        required -= {"lapw0_worker", "lcore_worker"}
        allowed |= {"lapw0_worker", "lcore_worker"}
    if set(config) - allowed != required:
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
    policy = detect_policy(
        reference, case, override=config.get("physical_policy"),
        electronic_model=electronic_model, native_soc=soc,
    )
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
    if soc:
        suffixes += ["inso","in0abp"]
        if policy.magnetic:
            raise NotImplementedError("spin-polarized SOC response")
    captured = Path(config["captured_soc"]) if "captured_soc" in config else None
    if captured is not None:
        suffixes += ["vsp", "vns"]
    names = (
        [case + "." + s for s in suffixes]
        + ["EF.dat", "projectorw.dat", "Eorb_imp.dat", "params.dat"]
        + [p.name for p in reference.glob("sig.inp*")]
    )
    if soc and config.get("preserve_soc_normalization", False):
        names.append("SOlapm.dat")
    files = {name: reference / name for name in names if (reference / name).is_file()}
    before = fingerprint_files(files)
    for name, path in files.items():
        shutil.copy2(path, work / name)
    definitions = ["lapw1"] + (["lapwso"] if soc else [])
    if electronic_model == "dmft":
        definitions += ["dmft2"] + (["dmft2dn"] if policy.magnetic else [])
    for label in definitions:
        (work / (label + ".def")).write_text(
            local_definition(Path(config[label + "_definition"]).read_text())
        )
    for suffix in ["in0", "in2", "in2c"]:
        p = work / (case + "." + suffix)
        if p.is_file() and p.stat().st_size:
            lines = p.read_text().splitlines()
            # WIEN reads MODUS as five characters: TOT needs two spaces.
            # Merely replacing STRFT with TOT can leave "TOT (", which
            # silently skips the muffin-tin valence density in LAPW2.
            fields = lines[0].split(maxsplit=1)
            lines[0] = "TOT  " + (fields[1] if len(fields) > 1 else "")
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
    for program in ([] if captured is not None else ["lapw0", "lcore"]):
        run(program + "_definition", x, [program, "-f", case, "-d"])
    if captured is None:
        run("lapw0", config["lapw0_worker"], ["lapw0.def"])
    if soc:
        # Keep both scalar and spinor eigenvector streams concatenated in this
        # private snapshot. The response distributes their k points itself.
        from lapw_soc_rank import enforce_ranks, parent_lo_filter
        from wstruct import Struct
        struct=Struct()
        struct.ReadStruct(str(work / (case+".struct")))
        layout_file=reference / "spinor_ranks.json"
        layout=json.loads(layout_file.read_text()) if layout_file.is_file() else None
        run("lapw1", config["lapw1_worker"], ["lapw1.def"],
            {"DFTKS_OVERLAP_DUMP": str(work/"native")}, "LAPW1")
        if layout is not None:
            defs=entries((work/"lapw1.def").read_text())
            enforce_ranks(work/defs[10], work/defs[11], struct.nat, layout["scalar"])
        legacy_matching = parent_lo_filter(reference / (case + ".scfso"))
        soc_env = {"EDMFT_SOC_RESPONSE_GEOMETRY": str(work/"soc.geometry"),
                   "EDMFT_SOC_LEGACY_LO_FILTER": str(int(legacy_matching))}
        if captured is not None:
            soc_env["EDMFT_SOC_GEOMETRY_ONLY"] = "1"
        run("lapwso", config["lapwso_worker"], ["lapwso.def"], soc_env, "LAPWSO")
        if layout is not None and captured is None:
            defs=entries((work/"lapwso.def").read_text())
            for vector,energy,norm in ((41,51,45),(42,52,46)):
                enforce_ranks(work/defs[vector],work/defs[energy],struct.nat,layout["spinor"],work/defs[norm])
    elif lapw1_mpi:
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
    if captured is not None:
        if not (captured / "response.spin").is_file() or not list(captured.glob("response.operator.k????????.bin")):
            raise ValueError("ordinary SOC FOR did not export its spinor moments")
        for exported in captured.glob("response.*"):
            shutil.move(str(exported),work / exported.name)
    elif electronic_model == "dft":
        # Native LAPW2 supplies its own Mermin occupations and valence charge.
        # Its DEF and all outputs are isolated in this response workspace.
        # Follow the actual last LAPW1 definition; an old in1c file can
        # remain in a real calculation and must not silently switch LAPW2.
        complex_vectors = entries((work / "lapw1.def").read_text())[5].endswith(".in1c")
        run("lapw2_definition", x, ["lapw2", "-f", case, "-d"] + (["-c"] if complex_vectors else []))
        run("lapw2", Path(config["wien_root"]) / ("lapw2c" if complex_vectors else "lapw2"), ["lapw2.def"])
    else:
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
                export_env = {"EDMFT_RESPONSE_EXPORT": str(work/("response"+suffix))}
                if soc and config.get("preserve_soc_normalization", False):
                    export_env["EDMFT_SOC_OPERATOR_READOUT"] = "1"
                run(label, config["dmft2"], [definition], export_env, "DMFT2")
            run_channels(definitions, policy, export, common_fermi, export=True)
    if captured is None:
        run("lcore", config["lcore_worker"], ["lcore.def"])
        prepare_native_potential_core(work, case, work / "potential_core.npz")
    if before != fingerprint_files(files):
        raise RuntimeError("parent inputs changed during preparation")
    missing_export = (
        any(not (work / ("response"+("_dn" if ch else "")+".meta")).is_file() for ch in policy.channels)
        if electronic_model == "dmft" else
        not (work / (case + ".weight")).is_file()
    )
    if missing_export or not list(work.glob("native.ham.k*.bin")):
        raise RuntimeError("workers did not export native response data")
    report = dict(
        case=case,
        electronic_model=electronic_model,
        physical_policy=policy.record(),
        parent_fingerprints=before,
        config=config,
        timings=times,
        spinor=soc,
        scope="one central snapshot, no strained SCF loop",
        stress_applied=False,
    )
    if soc:
        report["soc_lo_matching"] = "LEGACY_COMPONENT_CUTOFF" if legacy_matching else "LINEAR"
    if soc and (reference / "SOlapm.dat").is_file():
        from frozen_sigma_response import Metadata
        from lapw_soc_electronic import saved_normalization_error
        metadata = Metadata.read(str(work / "response.meta"))
        report["maximum_parent_normalization_change"] = saved_normalization_error(reference / "SOlapm.dat", metadata)
    (work / "preparation.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(prepare(json.loads(args.config.read_text())), indent=2), flush=True
    )
