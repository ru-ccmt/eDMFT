#!/usr/bin/env python3
"""Build and check the native discrete radial tangent/adjoint in isolation.

No installed executable or material file is modified. Without --work, build
files are temporary. Optional --case-directory entries read actual spherical
potentials and test l=0,2,4 at a stated test energy (0.17 Ha); they do NOT
compute material stresses or assume these are the case's linearization energies.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile


def run(compiler, root, work, case_directories):
    tests = Path(__file__).resolve().parent
    atpar = root / "src/dmft2/atpar"
    flags = ["-O2", "-fcheck=bounds,pointer,mem,do", "-ffpe-trap=invalid,zero,overflow"]
    results = {}
    programs = {
        "test_outwin_response": [],
        "test_atpar_response": [
            "atpar",
            "rint13",
            "abc",
            "diracout",
            "dergl",
            "inouh",
            "inth",
        ],
        "test_lo_response": ["abc"],
    }
    if case_directories:
        programs["radial_material_driver"] = []
    for name, dependencies in programs.items():
        command = [
            compiler,
            *flags,
            str(atpar / "outwin.f90"),
            str(atpar / "outwin_response.f90"),
            *[str(atpar / (dep + ".f90")) for dep in dependencies],
            str(tests / (name + ".f90")),
            "-o",
            str(work / name),
        ]
        built = subprocess.run(command, cwd=work, capture_output=True, text=True)
        (work / (name + ".build.log")).write_text(built.stdout + built.stderr)
        built.check_returncode()
        if name != "radial_material_driver":
            tested = subprocess.run(
                [str(work / name)], cwd=work, capture_output=True, text=True
            )
            (work / (name + ".log")).write_text(tested.stdout + tested.stderr)
            tested.check_returncode()
            if "PASS" not in tested.stdout:
                raise RuntimeError("missing success marker: " + name)
            results[name] = tested.stdout.strip().splitlines()[-1]
    if case_directories:
        sys.path.insert(0, str(root / "src/python"))
        from wstruct import Struct

        material_rows = []
        for directory in case_directories:
            structures = list(directory.glob("*.struct"))
            if len(structures) != 1:
                raise ValueError(
                    "expected exactly one case.struct in " + str(directory)
                )
            structure = structures[0]
            case = structure.stem
            s = Struct()
            s.ReadStruct(str(structure))
            mode = int("RELA" in structure.read_text().splitlines()[2])
            blocks = re.split(
                r"ATOMNUMBER\s*=\s*\d+", (directory / (case + ".vsp")).read_text()
            )[1:]
            if len(blocks) != s.nat:
                raise ValueError("VSP atom count differs from structure")
            for atom, n in enumerate(s.jrj):
                data = re.split(r"VLM\(R\)\s+FOR L\s+0\s+M=\s*0", blocks[atom])[1]
                values = re.findall(r"[-+]?\d+\.\d+[DE][-+]\d+", data)[:n]
                if len(values) != n:
                    raise ValueError("truncated spherical potential")
                potential = [float(v.replace("D", "E")) / 2 for v in values]
                dx = math.log(s.rmt[atom] / s.r0[atom]) / (n - 1)
                for ell in [0, 2, 4]:
                    path = work / (case + f"_atom{atom+1}_l{ell}.input")
                    path.write_text(
                        f"{n} {mode} {s.r0[atom]:.17g} {dx:.17g} 0.17 {ell} {s.Znuc[atom]}\n"
                        + "\n".join(f"{v:.17g}" for v in potential)
                        + "\n"
                    )
                    result = subprocess.run(
                        [str(work / "radial_material_driver"), str(path)],
                        capture_output=True,
                        text=True,
                    )
                    result.check_returncode()
                    fd, adj, forward, reverse = map(float, result.stdout.split())
                    material_rows.append(
                        dict(
                            case=case,
                            atom=atom + 1,
                            ell=ell,
                            rela=bool(mode),
                            energy_Ha=0.17,
                            finite_difference_relative=fd,
                            transpose_relative=adj,
                            forward_seconds=forward,
                            adjoint_seconds=reverse,
                        )
                    )
        results["material_radial_controls"] = material_rows
    results["source_sha256"] = {
        name: hashlib.sha256((atpar / name).read_bytes()).hexdigest()
        for name in [
            "outwin.f90",
            "outwin_response.f90",
            "atpar.f90",
            "rint13.f90",
            "abc.f90",
        ]
    }
    (work / "report.json").write_text(json.dumps(results, indent=2) + "\n")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default=shutil.which("gfortran"))
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[3]
    )
    parser.add_argument("--work", type=Path)
    parser.add_argument("--case-directory", type=Path, action="append", default=[])
    args = parser.parse_args()
    if args.compiler is None:
        parser.error("gfortran not found; supply --compiler")
    if args.work:
        args.work = args.work.resolve()
        args.work.mkdir(parents=True, exist_ok=False)
        result = run(args.compiler, args.root.resolve(), args.work, args.case_directory)
    else:
        with tempfile.TemporaryDirectory(prefix="edmft-radial-response-") as folder:
            result = run(
                args.compiler, args.root.resolve(), Path(folder), args.case_directory
            )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
