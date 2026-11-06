#!/usr/bin/env python3
"""Analyze finite-strain LAPW radial-function dumps.

The dumps are produced only when LAPW2 is run with
``WIEN_STRESS_RADIAL_DUMP=FILE``.  This utility forms the central finite
difference of the radial functions from converged ``-eps`` and ``+eps``
calculations.  It is an incomplete-basis-set (IBS) diagnostic precursor, not
a stress correction: no result printed here is added to the stress tensor.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys


COLUMNS = (
    "potential",
    "u_large",
    "u_small",
    "udot_large",
    "udot_small",
    "du_dr_large",
    "du_dr_small",
    "dudot_dr_large",
    "dudot_dr_small",
)


class RadialResponseError(RuntimeError):
    pass


def read_dump(path: Path) -> dict[tuple[int, int], tuple[list[float], list[list[float]]]]:
    records: dict[tuple[int, int], list[tuple[int, float, list[float]]]] = {}
    if not path.is_file():
        raise RadialResponseError(f"missing dump: {path}")
    for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        if not line.startswith(":RADIAL_FD"):
            continue
        fields = line.split()
        if len(fields) != 14:
            raise RadialResponseError(
                f"{path}:{line_number}: expected 14 fields, found {len(fields)}"
            )
        atom, ell, radial_index = (int(fields[i]) for i in range(1, 4))
        radius = float(fields[4].replace("D", "E"))
        values = [float(value.replace("D", "E")) for value in fields[5:]]
        records.setdefault((atom, ell), []).append((radial_index, radius, values))
    if not records:
        raise RadialResponseError(f"no :RADIAL_FD records in {path}")

    result = {}
    for key, rows in records.items():
        rows.sort()
        expected = list(range(1, len(rows) + 1))
        if [row[0] for row in rows] != expected:
            raise RadialResponseError(f"nonconsecutive radial mesh for atom/l={key}")
        result[key] = ([row[1] for row in rows], [row[2] for row in rows])
    return result


def trapz_square(radius: list[float], first: list[float], second: list[float]) -> float:
    total = 0.0
    for i in range(len(radius) - 1):
        left = first[i] * first[i] + second[i] * second[i]
        right = first[i + 1] * first[i + 1] + second[i + 1] * second[i + 1]
        total += 0.5 * (radius[i + 1] - radius[i]) * (left + right)
    return total


def compatible_mesh(
    reference: dict[tuple[int, int], tuple[list[float], list[list[float]]]],
    candidate: dict[tuple[int, int], tuple[list[float], list[list[float]]]],
    label: str,
) -> None:
    if reference.keys() != candidate.keys():
        raise RadialResponseError(f"{label}: atom/l channels do not match")
    for key in reference:
        r0 = reference[key][0]
        r1 = candidate[key][0]
        if len(r0) != len(r1):
            raise RadialResponseError(f"{label}: radial size differs for atom/l={key}")
        for i, (left, right) in enumerate(zip(r0, r1), 1):
            if not math.isclose(left, right, rel_tol=2.0e-12, abs_tol=1.0e-14):
                raise RadialResponseError(
                    f"{label}: radial mesh differs for atom/l={key}, point {i}"
                )


def analyze(
    minus_path: Path,
    plus_path: Path,
    strain: float,
    center_path: Path | None,
) -> dict[str, object]:
    minus = read_dump(minus_path)
    plus = read_dump(plus_path)
    compatible_mesh(minus, plus, "plus")
    center = read_dump(center_path) if center_path is not None else None
    if center is not None:
        compatible_mesh(minus, center, "center")

    channels = []
    for atom, ell in sorted(minus):
        radius, minus_values = minus[(atom, ell)]
        plus_values = plus[(atom, ell)][1]
        response = [
            [(p - m) / (2.0 * strain) for m, p in zip(mrow, prow)]
            for mrow, prow in zip(minus_values, plus_values)
        ]
        if center is None:
            base = [
                [0.5 * (m + p) for m, p in zip(mrow, prow)]
                for mrow, prow in zip(minus_values, plus_values)
            ]
            nonlinearity = None
        else:
            base = center[(atom, ell)][1]
            numerator = sum(
                (p + m - 2.0 * c) ** 2
                for mrow, prow, crow in zip(minus_values, plus_values, base)
                for m, p, c in zip(mrow, prow, crow)
            )
            denominator = sum(
                (p - m) ** 2
                for mrow, prow in zip(minus_values, plus_values)
                for m, p in zip(mrow, prow)
            )
            nonlinearity = math.sqrt(numerator / denominator) if denominator else 0.0

        u_norm = math.sqrt(
            max(0.0, trapz_square(radius, [row[1] for row in base], [row[2] for row in base]))
        )
        du_norm = math.sqrt(
            max(
                0.0,
                trapz_square(
                    radius,
                    [row[1] for row in response],
                    [row[2] for row in response],
                ),
            )
        )
        udot_norm = math.sqrt(
            max(0.0, trapz_square(radius, [row[3] for row in base], [row[4] for row in base]))
        )
        dudot_norm = math.sqrt(
            max(
                0.0,
                trapz_square(
                    radius,
                    [row[3] for row in response],
                    [row[4] for row in response],
                ),
            )
        )
        channels.append(
            {
                "atom": atom,
                "l": ell,
                "mesh_points": len(radius),
                "rmt_bohr": radius[-1],
                "u_response_norm": du_norm,
                "u_relative_response": du_norm / u_norm if u_norm else None,
                "udot_response_norm": dudot_norm,
                "udot_relative_response": dudot_norm / udot_norm if udot_norm else None,
                "u_boundary_response": response[-1][1],
                "du_dr_boundary_response": response[-1][5],
                "potential_boundary_response": response[-1][0],
                "central_nonlinearity_ratio": nonlinearity,
            }
        )

    return {
        "strain": strain,
        "minus": str(minus_path),
        "plus": str(plus_path),
        "center": None if center_path is None else str(center_path),
        "interpretation": (
            "finite-strain radial-response diagnostic only; not an IBS stress "
            "contraction and not a stress correction"
        ),
        "channels": channels,
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minus", type=Path, required=True)
    parser.add_argument("--plus", type=Path, required=True)
    parser.add_argument("--strain", type=float, required=True)
    parser.add_argument("--center", type=Path)
    parser.add_argument("--json", type=Path, help="write complete machine-readable output")
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    try:
        if not math.isfinite(options.strain) or options.strain <= 0.0:
            raise RadialResponseError("--strain must be finite and positive")
        result = analyze(options.minus, options.plus, options.strain, options.center)
        print("RADIAL RESPONSE DIAGNOSTIC ONLY -- DO NOT ADD TO STRESS")
        print(" atom   l       ||du/deps||    relative      u(R) response   nonlinearity")
        for row in result["channels"]:
            nonlinearity = row["central_nonlinearity_ratio"]
            nonlinear_text = "       n/a" if nonlinearity is None else f"{nonlinearity:10.3e}"
            relative = row["u_relative_response"]
            relative_text = "       n/a" if relative is None else f"{relative:10.3e}"
            print(
                f"{row['atom']:5d} {row['l']:3d} {row['u_response_norm']:16.8e} "
                f"{relative_text} {row['u_boundary_response']:16.8e} {nonlinear_text}"
            )
        if options.json is not None:
            options.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
            print(f"wrote {options.json}")
        return 0
    except RadialResponseError as error:
        print(f"radial_response_fd.py: error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
