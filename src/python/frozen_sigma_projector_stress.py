#!/usr/bin/env python3
"""Contract a material-potential response with the normalized-projector adjoint.

The NPZ input contains potential[k,i,j] (Ry per hydrostatic linear strain),
response[k,i,j] (the SUM of the window and normalization adjoints, including
all k/spin/frequency weights), and volume (Bohr^3).  Both arrays must use the
same central orthonormal band gauge.  This program reports the hydrostatic
MEAN tensile-stress correction, already divided by three.  It does not turn
one hydrostatic derivative into a general anisotropic stress tensor.

The matrices must already contain the material-coordinate transformation.
Raw differences of WIEN VNS coefficients are not this potential operator.
"""
import argparse
import json
from pathlib import Path

import numpy as np

RY_BOHR3_TO_GPA = 14710.507641078


def contract(potential, response, volume):
    if (
        potential.ndim != 3
        or potential.shape != response.shape
        or potential.shape[-1] != potential.shape[-2]
    ):
        raise ValueError("potential and adjoint must have matching k/band dimensions")
    if not np.isfinite(volume) or volume <= 0:
        raise ValueError("volume must be positive")
    for name, matrix in [("potential", potential), ("adjoint", response)]:
        if not np.all(np.isfinite(matrix)):
            raise ValueError(f"nonfinite {name}")
        error = np.max(abs(matrix - matrix.swapaxes(-1, -2).conj()))
        if error > 1e-9 * max(1.0, float(np.max(abs(matrix)))):
            raise ValueError(f"non-Hermitian {name}")
    derivative = float(np.einsum("kij,kji->", response, potential).real)
    correction = derivative / (3 * volume) * RY_BOHR3_TO_GPA
    return {
        "projector_energy_derivative_Ry": derivative,
        "mean_tensile_correction_GPa": correction,
        "pressure_correction_GPa": -correction,
        "volume_Bohr3": volume,
        "convention": "hydrostatic mean; already divided by three",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrices", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--analytic", type=float, help="baseline mean tensile stress in GPa"
    )
    parser.add_argument(
        "--numeric", type=float, help="reference mean tensile stress in GPa"
    )
    args = parser.parse_args()
    with np.load(args.matrices, allow_pickle=False) as data:
        result = contract(data["potential"], data["response"], float(data["volume"]))
    if args.analytic is not None:
        corrected = args.analytic + result["mean_tensile_correction_GPa"]
        result["corrected_mean_tensile_GPa"] = corrected
        if args.numeric is not None:
            result["corrected_minus_numeric_GPa"] = corrected - args.numeric
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
