#!/usr/bin/env python3
"""Test the native radial H/O tangent against independent ATPARN endpoints.

Use inputs from check_lapw_radial_response.py --work FIXTURES. This builds
only in a fresh scratch directory, and never changes an installed executable.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[3]
    )
    parser.add_argument("--module-directory", type=Path)
    parser.add_argument("--fixtures", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--compiler", default=shutil.which("gfortran"))
    args = parser.parse_args()
    root = args.root.resolve()
    atpar = root / "src/dmft2/atpar"
    module = args.module_directory or atpar
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    names = [
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
    sources = [atpar / (name + ".f90") for name in names]
    sources[2:2] = [
        module / "lapw_radial_response.f90",
        module / "lapw_radial_hamiltonian.f90",
    ]
    sources.append(Path(__file__).with_name("test_radial_hamiltonian.f90"))
    build = subprocess.run(
        [
            args.compiler,
            "-O2",
            "-fcheck=bounds,pointer,mem,do",
            "-ffpe-trap=invalid,zero,overflow",
            *map(str, sources),
            "-o",
            str(work / "check"),
        ],
        cwd=work,
        capture_output=True,
        text=True,
    )
    (work / "build.log").write_text(build.stdout + build.stderr)
    build.check_returncode()
    rows = []
    keys = [
        "central_relative",
        "H_fd_relative",
        "O_fd_relative",
        "adjoint_relative",
        "forward_seconds",
        "adjoint_seconds",
    ]
    inputs = sorted(args.fixtures.glob("*.input"))
    if not inputs:
        raise ValueError("no native fixtures found")
    for p in inputs:
        result = subprocess.run(
            [str(work / "check"), str(p.resolve())],
            cwd=work,
            capture_output=True,
            text=True,
        )
        (work / (p.stem + ".log")).write_text(result.stdout + result.stderr)
        result.check_returncode()
        values = list(map(float, result.stdout.split()))
        if len(values) != len(keys):
            raise RuntimeError("unexpected test output")
        rows.append(dict(case=p.stem, **dict(zip(keys, values))))
    report = dict(
        rows=rows,
        source_sha256={
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources
        },
        input_sha256={
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs
        },
        maximum_errors={key: max(row[key] for row in rows) for key in keys[:4]},
    )
    (work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            dict(
                cases=len(rows),
                maximum_errors=report["maximum_errors"],
                report=str(work / "report.json"),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
