#!/usr/bin/env python3
"""Fixed-G-vector finite-difference regression test for WIEN2k stress.

The matching optional LAPW1 patch is activated only by marker files:

  lapw1.fixgvec.write  write the ordinary central-cell G-vector lists
  lapw1.fixgvec.read   reuse those lists in a strained cell

This driver prepares the markers/reference and compares converged total-energy
derivatives with the analytic STR_DENS - STR_CORE - STR_VAL tensor.  It never
edits a WIEN2k source file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


RY_BOHR3_TO_GPA = 14710.5164057


def detect_case(directory: Path, explicit: str | None) -> str:
    if explicit:
        struct = directory / f"{explicit}.struct"
        if not struct.is_file():
            raise RuntimeError(f"missing structure file: {struct}")
        return explicit
    structs = sorted(directory.glob("*.struct"))
    if len(structs) != 1:
        names = ", ".join(path.name for path in structs) or "none"
        raise RuntimeError(
            f"cannot detect case in {directory}: found {names}; use --case"
        )
    return structs[0].stem


def last_float(path: Path, tag: str) -> float:
    if not path.is_file():
        raise RuntimeError(f"missing file: {path}")
    pattern = re.compile(rf"^{re.escape(tag)}.*?([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s*(?:\([^)]*\))?\s*$")
    value: float | None = None
    for line in path.read_text(errors="replace").splitlines():
        match = pattern.search(line)
        if match:
            value = float(match.group(1))
    if value is None:
        raise RuntimeError(f"tag {tag!r} not found in {path}")
    return value


def total_energy(directory: Path, case: str) -> float:
    return last_float(directory / f"{case}.scf", ":ENE  :")


def volume(directory: Path, case: str) -> float:
    return last_float(directory / f"{case}.scf", ":VOL  :")


def raw_energy_parts(directory: Path, case: str) -> dict[str, float]:
    scf = directory / f"{case}.scf"
    ene = last_float(scf, ":ENE  :")
    val = last_float(scf, ":SUM  :")
    density = last_float(scf, ":DEN  :")
    return {
        "valence_eigenvalue": val,
        "core_eigenvalue": ene - val - density,
        "density": density,
        "total": ene,
    }


def tensor_from_file(path: Path, tag: str) -> list[list[float]]:
    tensor = [[0.0] * 3 for _ in range(3)]
    rows: dict[int, list[float]] = {}
    pattern = re.compile(
        rf"^{re.escape(tag)}00([123]):\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)"
    )
    for line in path.read_text(errors="replace").splitlines():
        match = pattern.match(line)
        if match:
            rows[int(match.group(1)) - 1] = [float(match.group(i)) for i in range(2, 5)]
    if set(rows) != {0, 1, 2}:
        raise RuntimeError(f"incomplete {tag} tensor in {path}")
    for row, values in rows.items():
        tensor[row] = values
    return tensor


def analytic_tensor(directory: Path, case: str) -> tuple[list[list[float]], dict[str, list[list[float]]]]:
    parts = {
        "density": tensor_from_file(directory / f"{case}.scf0", ":STR_DENS"),
        "core": tensor_from_file(directory / f"{case}.scfc", ":STR_CORE"),
        "valence": tensor_from_file(directory / f"{case}.scf2", ":STR_VAL"),
    }
    total = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            total[i][j] = (
                parts["density"][i][j]
                - parts["core"][i][j]
                - parts["valence"][i][j]
            ) * RY_BOHR3_TO_GPA
    return total, parts


def matrix_sizes(path: Path) -> list[int]:
    if not path.is_file():
        raise RuntimeError(f"missing LAPW1 output: {path}")
    pattern = re.compile(r"Matrix size\s+(\d+)", re.IGNORECASE)
    return [int(match.group(1)) for match in pattern.finditer(path.read_text(errors="replace"))]


def reference_sizes(path: Path) -> list[int]:
    lines = path.read_text(errors="strict").splitlines()
    if not lines or lines[0].strip() != "LAPW1_FIXGVEC_V1":
        raise RuntimeError(f"invalid fixed-G reference: {path}")
    sizes: list[int] = []
    index = 1
    while index < len(lines):
        fields = lines[index].split()
        if len(fields) != 4:
            raise RuntimeError(f"invalid block header at line {index + 1} in {path}")
        count = int(fields[0])
        sizes.append(count)
        index += 1 + count
    if index != len(lines):
        raise RuntimeError(f"truncated fixed-G reference: {path}")
    return sizes


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cmd_prepare(args: argparse.Namespace) -> int:
    center = args.center.resolve()
    minus = args.minus.resolve()
    plus = args.plus.resolve()
    case = detect_case(center, args.case)
    for directory in (minus, plus):
        detect_case(directory, case)

    write_marker = center / "lapw1.fixgvec.write"
    read_marker = center / "lapw1.fixgvec.read"
    reference = center / "lapw1.fixgvec"
    if read_marker.exists() or write_marker.exists():
        raise RuntimeError(f"remove stale fixed-G marker(s) in {center}")

    write_marker.touch()
    try:
        command = [args.x_command, "lapw1"]
        if (center / f"{case}.in1c").is_file():
            command.append("-c")
        command.extend(["-f", case])
        subprocess.run(command, cwd=center, check=True)
    finally:
        if write_marker.exists():
            write_marker.unlink()

    sizes = reference_sizes(reference)
    for directory in (minus, plus):
        shutil.copy2(reference, directory / "lapw1.fixgvec")
        marker = directory / "lapw1.fixgvec.read"
        marker.touch(exist_ok=True)

    print(f"case={case}")
    print(f"reference={reference}")
    print(f"sha256={sha256(reference)}")
    print(f"k-points={len(sizes)}")
    print(f"basis range={min(sizes)}..{max(sizes)} plane waves")
    print("Staged lapw1.fixgvec.read in minus and plus directories.")
    return 0


def mode_components(mode: str) -> tuple[list[int], int]:
    if mode == "iso":
        return [0, 1, 2], 3
    if mode == "ab":
        return [0, 1], 2
    return [{"x": 0, "y": 1, "z": 2}[mode]], 1


def cmd_compare(args: argparse.Namespace) -> int:
    center = args.center.resolve()
    minus = args.minus.resolve()
    plus = args.plus.resolve()
    case = detect_case(center, args.case)
    for directory in (minus, plus):
        detect_case(directory, case)

    reference = center / "lapw1.fixgvec"
    expected = reference_sizes(reference)
    reference_hash = sha256(reference)
    staged_ok = all(
        (directory / "lapw1.fixgvec.read").is_file()
        and (directory / "lapw1.fixgvec").is_file()
        and sha256(directory / "lapw1.fixgvec") == reference_hash
        for directory in (minus, plus)
    )
    observed = {
        "center": matrix_sizes(center / f"{case}.output1"),
        "minus": matrix_sizes(minus / f"{case}.output1"),
        "plus": matrix_sizes(plus / f"{case}.output1"),
    }
    # ``Matrix size`` includes local-orbital columns, whereas the reference
    # stores only the interstitial plane waves.  Its absolute value therefore
    # need not equal the reference block length.  The invariant required by
    # the finite difference is that LAPW1 uses the same full matrix dimensions
    # in the center and both strained calculations.  COORS has a sizing and a
    # construction pass, so current output1 files normally contain two matrix
    # size records per k point.
    record_count_ok = all(
        len(values) in (len(expected), 2 * len(expected))
        for values in observed.values()
    )
    basis_ok = (
        staged_ok
        and record_count_ok
        and observed["center"] == observed["minus"] == observed["plus"]
    )

    omega = volume(center, case)
    e_minus = total_energy(minus, case)
    e_plus = total_energy(plus, case)
    indices, multiplicity = mode_components(args.mode)
    derivative = (e_plus - e_minus) / (2.0 * args.strain)
    numerical = derivative / (multiplicity * omega) * RY_BOHR3_TO_GPA

    analytic, analytic_parts = analytic_tensor(center, case)
    analytic_selected = sum(analytic[index][index] for index in indices) / multiplicity
    difference = analytic_selected - numerical

    minus_parts = raw_energy_parts(minus, case)
    plus_parts = raw_energy_parts(plus, case)
    raw_derivatives = {
        name: (plus_parts[name] - minus_parts[name])
        / (2.0 * args.strain * multiplicity * omega)
        * RY_BOHR3_TO_GPA
        for name in minus_parts
    }

    result = {
        "case": case,
        "mode": args.mode,
        "strain": args.strain,
        "volume_bohr3": omega,
        "basis_identical": basis_ok,
        "k_points": len(expected),
        "energy_minus_ry": e_minus,
        "energy_plus_ry": e_plus,
        "numerical_gpa": numerical,
        "analytic_gpa": analytic_selected,
        "analytic_minus_numerical_gpa": difference,
        "raw_numerical_parts_gpa": raw_derivatives,
        "analytic_tensor_gpa": analytic,
        "analytic_raw_groups_ry_bohr3": analytic_parts,
        "pass": basis_ok and abs(difference) <= args.tolerance,
    }

    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(f"case                         {case}")
        print(f"mode                         {args.mode}")
        print(f"strain                       {args.strain:.12e}")
        print(f"reference k-points           {len(expected)}")
        print(f"reference staged exactly     {'YES' if staged_ok else 'NO'}")
        print(f"fixed basis identical        {'YES' if basis_ok else 'NO'}")
        print(f"E(-)                         {e_minus:.12f} Ry")
        print(f"E(+)                         {e_plus:.12f} Ry")
        print(f"numerical stress             {numerical:.9f} GPa")
        print(f"analytic stress              {analytic_selected:.9f} GPa")
        print(f"analytic - numerical         {difference:+.9f} GPa")
        print("raw numerical energy derivatives (not analytic term counterparts):")
        for name, value in raw_derivatives.items():
            print(f"  {name:22s} {value:+.9f} GPa")
        print(f"result                       {'PASS' if result['pass'] else 'FAIL'}")
    return 0 if result["pass"] else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare", help="write and stage a central fixed-G reference")
    prepare.add_argument("--center", type=Path, required=True)
    prepare.add_argument("--minus", type=Path, required=True)
    prepare.add_argument("--plus", type=Path, required=True)
    prepare.add_argument("--case")
    prepare.add_argument("--x-command", default="x")
    prepare.set_defaults(func=cmd_prepare)

    compare = subparsers.add_parser("compare", help="compare converged numerical and analytic stress")
    compare.add_argument("--center", type=Path, required=True)
    compare.add_argument("--minus", type=Path, required=True)
    compare.add_argument("--plus", type=Path, required=True)
    compare.add_argument("--case")
    compare.add_argument("--strain", type=float, required=True)
    compare.add_argument("--mode", choices=("iso", "ab", "x", "y", "z"), required=True)
    compare.add_argument("--tolerance", type=float, default=0.1, help="PASS tolerance in GPa")
    compare.add_argument("--json", action="store_true")
    compare.set_defaults(func=cmd_compare)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        return args.func(args)
    except (OSError, RuntimeError, subprocess.CalledProcessError, ValueError) as exc:
        parser.error(str(exc))
        return 2


if __name__ == "__main__":
    sys.exit(main())
