"""Diagnostic only: compare matrix P^T rho with native spatial-density pairing.

The potential contains rV in spherical MT, V_lm in other MT channels, and
already-masked per-member IR coefficients. Density PW entries are star SUMS;
there is no additional star multiplicity. Only the raw VALENCE density is
compared here, before native mixer normalization and before adding core.
This is not an implementation of the complete Klein double-counting term.
"""

from pathlib import Path
import json
import sys
import numpy as np
from check_material import root
from lapw_potential_core_state import load_potential_core_state
from lapw_hartree_periodic import native_charge_weights
from lapw_charge_response import native_charge_covectors, native_mixer_inversion
from wien_response_density import Density

case = sys.argv[1]
source = Path("/private/tmp/edmft_full_k_response_20261108") / case / "charge/0"
pc = load_potential_core_state(
    "/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/"
    + case
    + "_potential_core.npz"
)
rho = pc.input_density.zeros().add(
    Density(source / (case + ".clmval"), source / (case + ".struct"))
)
g = np.zeros(pc.hartree.size)
for (a, l, m), sl in zip(pc.hartree.keys, pc.hartree.slices):
    radius = pc.hartree.integrations[a].r
    g[sl] = pc.structure.mult[a] * native_charge_weights(radius) * rho.rad[a, l, m]
    if (l, m) == (0, 0):
        g[sl] /= radius
volume = float(pc.structure.a * pc.structure.b * pc.structure.c)
g[pc.hartree.pw_offset :] = volume * np.array(
    [rho.pw[tuple(gv)] for gv in pc.hartree.stars.representatives], complex
).view(float)
field = np.zeros_like(g)
for (a, l, m), sl in zip(pc.hartree.keys, pc.hartree.slices):
    if (l, m) == (0, 0):
        field[sl] = pc.hartree.integrations[a].r
constant = np.zeros(pc.hartree.stars.nstar, complex)
constant[0] = 1
field[pc.hartree.pw_offset :] = pc.hartree.reanalysis.action(constant).view(float)
mt, ir = native_charge_covectors(
    rho, pc.structure, inversion=native_mixer_inversion(pc.structure)
)
native_number = float((mt + ir) @ rho.vector())
field_number = float(g @ field)
assert abs(native_number - field_number) < 1e-10, (native_number, field_number)
matrix = np.load(root / "three_k" / case / "physical_potential_covectors.npz")[
    "physical"
]
number = json.loads((root / "three_k" / case / "physical_potential.json").read_text())[
    "number"
]
direction = np.load(
    "/private/tmp/edmft_hartree_periodic_20261108/" + case + "_direction.npz"
)["density"]
potential = pc.potential_operator.matvec(direction)
induced = np.load(
    Path("/private/tmp/edmft_strain_charge_20261109/three_k")
    / case
    / "analytic_induced_potential.npy"
)
report = dict(
    case=case,
    raw_spatial_valence_number=native_number,
    constant_field_pairing=field_number,
    band_number=number,
    band_minus_spatial_number=number - native_number,
    scope="raw valence spatial/matrix pairing diagnostic, not full Klein cancellation",
    comparison={},
)
for label, v in [("potential_probe", potential), ("screened_induced", induced)]:
    report["comparison"][label] = dict(
        matrix_Ry=float(matrix @ v),
        spatial_Ry=float(g @ v),
        difference_Ry=float((matrix - g) @ v),
        MT_difference_Ry=float(
            (matrix - g)[: pc.hartree.pw_offset] @ v[: pc.hartree.pw_offset]
        ),
        IR_difference_Ry=float(
            (matrix - g)[pc.hartree.pw_offset :] @ v[pc.hartree.pw_offset :]
        ),
    )
np.savez_compressed(
    root / "three_k" / case / "spatial_pairing_covectors.npz",
    spatial=g,
    matrix_minus_spatial=matrix - g,
)
(root / "three_k" / case / "spatial_pairing.json").write_text(
    json.dumps(report, indent=2) + "\n"
)
print(json.dumps(report), flush=True)
