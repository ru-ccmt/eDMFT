#!/usr/bin/env python3
"""Independent ATPARN endpoint and adjoint checks of radial basis assembly.

Only scratch build files are written. Material tests use spherical potentials,
linearization energies/APW flags/ordinary LO slots from the energy-file header,
and the on-disk projector where present. Header energies have file precision
(usually five decimal places); vector files are not needed. These tests do not
compute stresses. Spherical Bessel matching data are sampled at three qR values.
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

import numpy as np
from scipy.special import spherical_jn


def run(root, module, tests, work, compiler, directories):
    atpar = root / "src/dmft2/atpar"
    sources = [
        atpar / (s + ".f90")
        for s in [
            "outwin",
            "outwin_response",
            "atpar",
            "rint13",
            "abc",
            "diracout",
            "dergl",
            "inouh",
            "inth",
        ]
    ]
    sources += [module, tests / "test_lapw_radial_response.f90"]
    executable = work / "test_lapw_radial_response"
    build = subprocess.run(
        [
            compiler,
            "-O2",
            "-fcheck=bounds,pointer,mem,do",
            "-ffpe-trap=invalid,zero,overflow",
            *map(str, sources),
            "-o",
            str(executable),
        ],
        cwd=work,
        capture_output=True,
        text=True,
    )
    (work / "build.log").write_text(build.stdout + build.stderr)
    build.check_returncode()
    cases = []
    n = 401
    r0, rmt, z = 1e-4, 2.1, 26.0
    dx = math.log(rmt / r0) / (n - 1)
    r = r0 * np.exp(dx * np.arange(n))
    # LAPW with two ordinary LOs and APW+lo with/without a higher LO.
    for mode in (0, 1):
        for lapw, energies in [
            (True, []),
            (True, [-1.2, 0.8]),
            (False, [0.34]),
            (False, [0.34, -1.2]),
        ]:
            cases.append(
                dict(
                    name=f"synthetic_{mode}_{int(lapw)}_{len(energies)}",
                    mode=mode,
                    ell=2,
                    lapw=lapw,
                    energy=0.34,
                    elo=energies,
                    z=z,
                    r0=r0,
                    rmt=rmt,
                    dx=dx,
                    potential=-2 * z * np.exp(-r) + 0.3 * r * r,
                    projector=np.column_stack((r**3 * np.exp(-r), r * r * np.exp(-r))),
                    projector_source="synthetic; finite at RMT",
                )
            )
    provenance = {}
    sys.path.insert(0, str(root / "src/python"))
    from wstruct import Struct

    for directory in directories:
        structures = list(directory.glob("*.struct"))
        if len(structures) != 1:
            raise ValueError("need exactly one structure in " + str(directory))
        structure = structures[0]
        case = structure.stem
        s = Struct()
        s.ReadStruct(str(structure))
        mode = int("RELA" in structure.read_text().splitlines()[2])
        energy_file = directory / (case + ".energy")
        if not energy_file.is_file():
            energy_file = directory / (case + ".energy_0")
        lines = energy_file.read_text().splitlines()
        vsp = directory / (case + ".vsp")
        blocks = re.split(r"ATOMNUMBER\s*=\s*\d+", vsp.read_text())[1:]
        projectors = {}
        projector_file = directory / "projectorw.dat"
        if projector_file.is_file():
            with projector_file.open() as stream:
                count = int(stream.readline().split()[1])
                for _ in range(count):
                    nr, nr0, atom, ell = map(int, stream.readline().split()[1:])
                    rows = np.array(
                        [list(map(float, stream.readline().split())) for _ in range(nr)]
                    )
                    projectors[(atom - 1, ell)] = rows[:nr0, 1:3]
        for path in [structure, energy_file, vsp, projector_file]:
            if path.is_file():
                provenance[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        for atom, n in enumerate(s.jrj):
            header = [
                float(lines[2 * atom][i : i + 9])
                for i in range(0, len(lines[2 * atom]), 9)
            ]
            lo = [
                float(lines[2 * atom + 1][i : i + 9])
                for i in range(0, len(lines[2 * atom + 1]), 9)
            ]
            # Standard energy header: lomax=3, three LO slots; the final
            # slot denotes DIRACOUT. Refuse active special LOs explicitly.
            if len(lo) != 12:
                raise ValueError("unsupported energy-header LO shape")
            lo = np.array(lo).reshape(3, 4)
            data = re.split(r"VLM\(R\)\s+FOR L\s+0\s+M=\s*0", blocks[atom])[1]
            values = re.findall(r"[-+]?\d+\.\d+[DE][-+]\d+", data)[:n]
            potential = np.array([float(v.replace("D", "E")) for v in values])
            r0, rmt, z = s.r0[atom], s.rmt[atom], s.Znuc[atom]
            dx = math.log(rmt / r0) / (n - 1)
            r = r0 * np.exp(dx * np.arange(n))
            for ell in range(4):
                if lo[-1, ell] < 995:
                    raise ValueError("DIRACOUT response not supported")
                selected = lo[:, ell][lo[:, ell] < 995].tolist()
                if any(lo[j, ell] >= 995 for j in range(len(selected))):
                    raise ValueError("noncompact LO slots")
                lapw = header[ell] <= 150
                actual = (atom, ell) in projectors
                projector = projectors.get(
                    (atom, ell),
                    np.column_stack((r ** (ell + 1) * np.exp(-r), r**ell * np.exp(-r))),
                )
                if projector.shape != (n, 2):
                    raise ValueError("projector radial dimensions differ")
                cases.append(
                    dict(
                        name=f"{case}_atom{atom+1}_l{ell}",
                        mode=mode,
                        ell=ell,
                        lapw=lapw,
                        energy=header[ell] - (0 if lapw else 200),
                        elo=selected,
                        z=z,
                        r0=r0,
                        rmt=rmt,
                        dx=dx,
                        potential=potential,
                        projector=projector,
                        projector_source=(
                            str(projector_file)
                            if actual
                            else "synthetic; finite at RMT"
                        ),
                    )
                )
    rows = []
    keys = [
        "central_relative",
        "basis_fd_relative",
        "forms_fd_relative",
        "basis_adjoint_relative",
        "forms_adjoint_relative",
        "wronskian",
        "wronskian_directional_derivative",
        "forward_seconds",
        "adjoint_seconds",
    ]
    for case in cases:
        n = len(case["potential"])
        x = np.array([0.5, 1.4, 2.8])
        boundary = np.array(
            [
                spherical_jn(case["ell"], x),
                spherical_jn(case["ell"], x, derivative=True) * x / case["rmt"],
            ]
        )
        path = work / (case["name"] + ".input")
        with path.open("w") as f:
            f.write(
                f"{n} {case['mode']} {case['ell']} {int(case['lapw'])} {len(case['elo'])}\n"
            )
            f.write(
                " ".join(
                    format(case[k], ".17g") for k in ["r0", "dx", "rmt", "energy", "z"]
                )
                + "\n"
            )
            if case["elo"]:
                f.write(" ".join(format(e, ".17g") for e in case["elo"]) + "\n")
            for a in [boundary, case["potential"], case["projector"]]:
                np.savetxt(f, np.asarray(a).ravel(order="F"), fmt="%.17g")
        result = subprocess.run(
            [str(executable), str(path)], cwd=work, capture_output=True, text=True
        )
        (work / (case["name"] + ".log")).write_text(result.stdout + result.stderr)
        result.check_returncode()
        values = list(map(float, result.stdout.split()))
        if len(values) != len(keys):
            raise RuntimeError("unexpected native test output")
        rows.append(
            {k: v for k, v in case.items() if k not in ["potential", "projector"]}
            | dict(zip(keys, values))
        )
    report = dict(
        rows=rows,
        source_sha256={
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources
        },
        material_sha256=provenance,
    )
    (work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return {
        "cases": len(rows),
        "maximum_errors": {k: max(r[k] for r in rows) for k in keys[:5]},
        "maximum_action_seconds": max(r["forward_seconds"] for r in rows),
        "report": str(work / "report.json"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[3]
    )
    parser.add_argument("--module", type=Path)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--compiler", default=shutil.which("gfortran"))
    parser.add_argument("--case-directory", type=Path, action="append", default=[])
    args = parser.parse_args()
    args.work = args.work.resolve()
    args.work.mkdir(parents=True, exist_ok=False)
    module = args.module or args.root / "src/dmft2/atpar/lapw_radial_response.f90"
    print(
        json.dumps(
            run(
                args.root.resolve(),
                module.resolve(),
                Path(__file__).resolve().parent,
                args.work,
                args.compiler,
                args.case_directory,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
