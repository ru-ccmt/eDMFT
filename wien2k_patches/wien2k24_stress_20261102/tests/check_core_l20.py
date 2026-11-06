#!/usr/bin/env python3
"""Check the tetragonal l=2,m=0 core-stress angular algebra.

The lcore executable must first be run with WIEN_STRESS_DEBUG=1.  Its
output then contains both the raw radial integrals from Blaha et al.
Eqs. (E4)-(E5) and the normalized Cartesian contributions accumulated by
the code.  For a real Y20 potential and a spherical core density, direct
Cartesian integration of Eq. (E1) gives

    stored I1_x = -M sqrt(5)/(30 sqrt(pi)) I1
    stored I2_x = -M sqrt(5)/(10 sqrt(pi)) I2
    stored I{1,2}_z = -2 stored I{1,2}_x,

where M is the atomic multiplicity.  "stored" is emphasized because mixer
forms STR_DENS-STR_VAL-STR_CORE; the physical core correction has the
opposite sign.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


RADIAL_RE = re.compile(
    r"^:DBG_STR_CORE_RADIAL:\s+(\d+)\s+(-?\d+)\s+(-?\d+)"
    r"\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)"
)
ATOM_RE = re.compile(
    r"^:DBG_STR_CORE_ATOM:\s+(\d+)\s+(\d+)"
    r"\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)"
)
MULT_RE = re.compile(r"\bMULT=\s*(\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path,
                        help="lcore output file containing DBG_STR_CORE records")
    parser.add_argument("--struct", required=True, type=Path,
                        help="matching WIEN2k structure file")
    parser.add_argument("--tolerance", type=float, default=1.0e-12,
                        help="absolute tolerance in stress-times-volume units")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    multiplicities = [int(x) for x in MULT_RE.findall(args.struct.read_text())]
    if not multiplicities:
        raise SystemExit(f"no MULT records found in {args.struct}")

    radial: dict[int, tuple[float, float]] = {}
    cartesian: dict[tuple[int, int], tuple[float, float]] = {}
    for line in args.output.read_text().splitlines():
        match = RADIAL_RE.match(line)
        if match:
            atom, ell, emm = map(int, match.group(1, 2, 3))
            if (ell, emm) == (2, 0):
                radial[atom] = (float(match.group(4)), float(match.group(5)))
            continue
        match = ATOM_RE.match(line)
        if match:
            atom, index = map(int, match.group(1, 2))
            cartesian[atom, index] = (float(match.group(4)),
                                      float(match.group(5)))

    coefficient_i1 = math.sqrt(5.0) / (30.0 * math.sqrt(math.pi))
    coefficient_i2 = math.sqrt(5.0) / (10.0 * math.sqrt(math.pi))
    maximum_error = 0.0

    print("atom axis       predicted I1          code I1        predicted I2          code I2")
    for atom in sorted(radial):
        if atom > len(multiplicities):
            raise SystemExit(f"atom {atom} has no MULT record in {args.struct}")
        if (atom, 1) not in cartesian or (atom, 9) not in cartesian:
            raise SystemExit(f"missing Cartesian core records for atom {atom}")
        mult = multiplicities[atom - 1]
        integral_i1, integral_i2 = radial[atom]
        predicted_x = (-mult * coefficient_i1 * integral_i1,
                       -mult * coefficient_i2 * integral_i2)
        predicted_z = (-2.0 * predicted_x[0], -2.0 * predicted_x[1])
        for axis, index, predicted in (("x", 1, predicted_x),
                                       ("z", 9, predicted_z)):
            observed = cartesian[atom, index]
            maximum_error = max(maximum_error,
                                abs(observed[0] - predicted[0]),
                                abs(observed[1] - predicted[1]))
            print(f"{atom:4d} {axis:>4s} {predicted[0]:20.12e} {observed[0]:16.12e}"
                  f" {predicted[1]:20.12e} {observed[1]:16.12e}")

    if not radial:
        raise SystemExit("no l=2,m=0 DBG_STR_CORE_RADIAL records found")
    print(f"maximum absolute error: {maximum_error:.6e}")
    if maximum_error > args.tolerance:
        print(f"FAIL: error exceeds tolerance {args.tolerance:.6e}")
        return 1
    print("PASS: l=2,m=0 core angular factors agree with direct Cartesian integration")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

