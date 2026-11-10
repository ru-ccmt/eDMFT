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

LO_COLUMNS = (
    "u_large",
    "u_small",
    "du_dr_large",
    "du_dr_small",
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


def read_lo_dump(
    path: Path,
) -> tuple[
    dict[tuple[int, int, int], tuple[list[float], list[list[float]]]],
    dict[tuple[int, int, int], tuple[bool, float, float, float, float]],
]:
    """Read bare LO radial functions and their ABC construction data."""
    records: dict[tuple[int, int, int], list[tuple[int, float, list[float]]]] = {}
    coefficients: dict[tuple[int, int, int], tuple[bool, float, float, float, float]] = {}
    for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        fields = line.split()
        if line.startswith(":RADIAL_LO_FD"):
            if len(fields) != 10:
                raise RadialResponseError(
                    f"{path}:{line_number}: expected 10 LO fields, found {len(fields)}"
                )
            atom, ell, local_orbital, radial_index = (
                int(fields[i]) for i in range(1, 5)
            )
            radius = float(fields[5].replace("D", "E"))
            values = [float(value.replace("D", "E")) for value in fields[6:]]
            records.setdefault((atom, ell, local_orbital), []).append(
                (radial_index, radius, values)
            )
        elif line.startswith(":RADIAL_LO_COEFF"):
            if len(fields) != 9:
                raise RadialResponseError(
                    f"{path}:{line_number}: expected 9 LO-coefficient fields, "
                    f"found {len(fields)}"
                )
            key = tuple(int(fields[i]) for i in range(1, 4))
            lapw = bool(int(fields[4]))
            values = tuple(float(value.replace("D", "E")) for value in fields[5:])
            coefficients[key] = (lapw, *values)

    if not records and not coefficients:
        return {}, {}
    if records.keys() != coefficients.keys():
        raise RadialResponseError(
            f"{path}: LO radial and coefficient channel sets do not match"
        )

    result = {}
    for key, rows in records.items():
        rows.sort()
        expected = list(range(1, len(rows) + 1))
        if [row[0] for row in rows] != expected:
            raise RadialResponseError(f"nonconsecutive LO radial mesh for atom/l/j={key}")
        result[key] = ([row[1] for row in rows], [row[2] for row in rows])
    return result, coefficients


def complete_lo(
    radial: list[list[float]],
    lo_radial: list[list[float]],
    coefficients: tuple[bool, float, float, float, float],
) -> list[list[float]]:
    """Form A*u+B*udot+C*u_LO and its stored radial derivative."""
    _lapw, a, b, c, _energy = coefficients
    return [
        [
            a * ordinary[1] + b * ordinary[3] + c * local[0],
            a * ordinary[2] + b * ordinary[4] + c * local[1],
            a * ordinary[5] + b * ordinary[7] + c * local[2],
            a * ordinary[6] + b * ordinary[8] + c * local[3],
        ]
        for ordinary, local in zip(radial, lo_radial)
    ]


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
    response_output: Path | None,
) -> dict[str, object]:
    minus = read_dump(minus_path)
    plus = read_dump(plus_path)
    compatible_mesh(minus, plus, "plus")
    center = read_dump(center_path) if center_path is not None else None
    if center is not None:
        compatible_mesh(minus, center, "center")

    minus_lo, minus_coefficients = read_lo_dump(minus_path)
    plus_lo, plus_coefficients = read_lo_dump(plus_path)
    center_lo, center_coefficients = (
        read_lo_dump(center_path) if center_path is not None else ({}, {})
    )
    if minus_lo.keys() != plus_lo.keys():
        raise RadialResponseError("plus: LO channels do not match minus")
    if center_path is not None and minus_lo.keys() != center_lo.keys():
        raise RadialResponseError("center: LO channels do not match minus")

    channels = []
    response_rows: list[str] = []
    for atom, ell in sorted(minus):
        radius, minus_values = minus[(atom, ell)]
        plus_values = plus[(atom, ell)][1]
        response = [
            [(p - m) / (2.0 * strain) for m, p in zip(mrow, prow)]
            for mrow, prow in zip(minus_values, plus_values)
        ]
        for radial_index, (radius_value, row) in enumerate(zip(radius, response), 1):
            response_rows.append(
                ":DR "
                f"{atom} {ell} {radial_index} {radius_value:.16e} "
                + " ".join(f"{row[index]:.16e}" for index in range(1, 9))
            )
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

    lo_channels = []
    for key in sorted(minus_lo):
        atom, ell, local_orbital = key
        ordinary_key = (atom, ell)
        minus_radius, minus_lo_radial = minus_lo[key]
        plus_radius, plus_lo_radial = plus_lo[key]
        if minus_radius != plus_radius:
            raise RadialResponseError(f"plus: LO radial mesh differs for {key}")
        if center_path is not None and minus_radius != center_lo[key][0]:
            raise RadialResponseError(f"center: LO radial mesh differs for {key}")

        complete_minus = complete_lo(
            minus[ordinary_key][1], minus_lo_radial, minus_coefficients[key]
        )
        complete_plus = complete_lo(
            plus[ordinary_key][1], plus_lo_radial, plus_coefficients[key]
        )
        complete_response = [
            [(p - m) / (2.0 * strain) for m, p in zip(mrow, prow)]
            for mrow, prow in zip(complete_minus, complete_plus)
        ]
        if center_path is None:
            complete_base = [
                [0.5 * (m + p) for m, p in zip(mrow, prow)]
                for mrow, prow in zip(complete_minus, complete_plus)
            ]
            nonlinearity = None
        else:
            complete_base = complete_lo(
                center[ordinary_key][1], center_lo[key][1], center_coefficients[key]
            )
            numerator = sum(
                (p + m - 2.0 * c) ** 2
                for mrow, prow, crow in zip(complete_minus, complete_plus, complete_base)
                for m, p, c in zip(mrow, prow, crow)
            )
            denominator = sum(
                (p - m) ** 2
                for mrow, prow in zip(complete_minus, complete_plus)
                for m, p in zip(mrow, prow)
            )
            nonlinearity = math.sqrt(numerator / denominator) if denominator else 0.0

        for radial_index, (radius_value, row) in enumerate(
            zip(minus_radius, complete_response), 1
        ):
            response_rows.append(
                ":DL "
                f"{atom} {ell} {local_orbital} {radial_index} "
                f"{radius_value:.16e} "
                + " ".join(f"{value:.16e}" for value in row)
            )
        lo_channels.append(
            {
                "atom": atom,
                "l": ell,
                "local_orbital": local_orbital,
                "lapw": center_coefficients.get(key, minus_coefficients[key])[0],
                "response_norm": math.sqrt(
                    max(
                        0.0,
                        trapz_square(
                            minus_radius,
                            [row[0] for row in complete_response],
                            [row[1] for row in complete_response],
                        ),
                    )
                ),
                "value_boundary_response": complete_response[-1][0],
                "slope_boundary_response": complete_response[-1][2],
                "central_nonlinearity_ratio": nonlinearity,
            }
        )

    if response_output is not None:
        response_output.write_text(
            "# EDMFT_RADIAL_RESPONSE_V1\n"
            f"# strain {strain:.16e}\n"
            "# :DR atom l ir r du_large du_small dudot_large dudot_small "
            "ddu_dr_large ddu_dr_small ddudot_dr_large ddudot_dr_small\n"
            "# :DL atom l jlo ir r dL_large dL_small dL_dr_large dL_dr_small\n"
            + "\n".join(response_rows)
            + "\n"
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
        "local_orbital_channels": lo_channels,
        "response_output": None if response_output is None else str(response_output),
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minus", type=Path, required=True)
    parser.add_argument("--plus", type=Path, required=True)
    parser.add_argument("--strain", type=float, required=True)
    parser.add_argument("--center", type=Path)
    parser.add_argument("--json", type=Path, help="write complete machine-readable output")
    parser.add_argument(
        "--response-output",
        type=Path,
        help="write the radial and complete-LO response consumed by dmft2 diagnostics",
    )
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    try:
        if not math.isfinite(options.strain) or options.strain <= 0.0:
            raise RadialResponseError("--strain must be finite and positive")
        result = analyze(
            options.minus,
            options.plus,
            options.strain,
            options.center,
            options.response_output,
        )
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
        if result["local_orbital_channels"]:
            print(" complete local-orbital response (boundary value/slope should vanish)")
            print(" atom   l  jlo       ||dL/deps||      dL(R)       dL'(R)  nonlinearity")
            for row in result["local_orbital_channels"]:
                nonlinearity = row["central_nonlinearity_ratio"]
                nonlinear_text = (
                    "       n/a" if nonlinearity is None else f"{nonlinearity:10.3e}"
                )
                print(
                    f"{row['atom']:5d} {row['l']:3d} {row['local_orbital']:4d} "
                    f"{row['response_norm']:16.8e} "
                    f"{row['value_boundary_response']:11.3e} "
                    f"{row['slope_boundary_response']:11.3e} {nonlinear_text}"
                )
        if options.response_output is not None:
            print(f"wrote {options.response_output}")
        if options.json is not None:
            options.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
            print(f"wrote {options.json}")
        return 0
    except RadialResponseError as error:
        print(f"radial_response_fd.py: error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
