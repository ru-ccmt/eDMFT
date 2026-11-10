"""Disentangle nonspherical wavefunction cutoff from radial-bank mismatch.

No native Hamiltonian or material input is changed. Evaluate potential
covectors on the SAME central KS coefficient density with l<=4 or l<=8,
using either the H/O radial bank or the ATPARN density bank. This local
matrix experiment diagnoses cancellation conventions, not a new stress.
"""

from pathlib import Path
import json
import sys
import numpy as np
from check_material import root
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_nonspherical_response import PreparedNonspherical

case = sys.argv[1]
source = Path("/private/tmp/edmft_full_k_response_20261108") / case / "charge/0"
pc = load_potential_core_state(
    "/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/"
    + case
    + "_potential_core.npz"
)
spatial = np.load(root / "three_k" / case / "spatial_pairing_covectors.npz")["spatial"]
matrix = np.load(root / "three_k" / case / "physical_potential_covectors.npz")[
    "physical"
]
induced = np.load(
    Path("/private/tmp/edmft_strain_charge_20261109/three_k")
    / case
    / "analytic_induced_potential.npy"
)
with load_native_valence_state(
    source,
    case,
    complex_vectors=(case == "SrVO3"),
    energy_policy="fixed",
    matrix_quadrature="wien2k24",
    library="/Users/haule/dbin/liblapw_response.so",
) as v:
    select, _ = pc.hartree.native_potential_output(
        tuple(v.bands.spherical.slices),
        {a: op.layout.harmonics for a, (op, ks) in v.bands.ns.items()},
        v.bands.ir[0].representatives,
        complex_interstitial=v.bands.input.complex_interstitial,
    )
    rho = v.model.central_density()
    mask = np.zeros_like(spatial)
    for key, sl in zip(pc.hartree.keys, pc.hartree.slices):
        if key[1:] != (0, 0):
            mask[sl] = 1
    report = dict(
        case=case,
        scope="fixed central coefficients; changed diagnostic contraction only",
        rows={},
    )
    for name, bank in [("matrix", v.contexts), ("density", v.projector_contexts)]:
        for cutoff in [4, 8]:
            ns = {}
            for a, (old, ks) in v.bands.ns.items():
                phi = {l: bank[a][l].primitive for l in range(cutoff + 1)}
                lo = {l: bank[a][l].lo_coefficients for l in phi}
                lapw = {l: bank[a][l].lapw for l in phi}
                entries = []
                assert len(v.direct.matching[a]) == len(ks)
                for channels in v.direct.matching[a]:
                    ch = {l: channels[l] for l in phi}
                    z = np.concatenate(
                        [
                            ch[l]
                            .value(phi[l][-1], lo[l], lapw=lapw[l])
                            .reshape(-1, v.model.eps.shape[1])
                            for l in phi
                        ]
                    )
                    entries.append((z, ch, lo, lapw))
                op = PreparedNonspherical(
                    phi,
                    v.contexts[a][0].overlap_metric,
                    old.layout.harmonics,
                    old.v,
                    entries,
                )
                ns[a] = op.potential_pullback([rho[k] for k in ks])
            native = v.bands.input.pack(
                np.zeros(v.bands.input.spherical_size),
                ns,
                np.zeros(v.bands.input.nstar, complex if case == "SrVO3" else float),
            )
            g = select.rmatvec(native)
            difference = g - mask * spatial
            label = f"{name}_l{cutoff}"
            report["rows"][label] = dict(
                difference_norm=float(np.linalg.norm(difference)),
                induced_difference_Ry=float(difference @ induced),
                induced_contraction_Ry=float(g @ induced),
                versus_native_matrix_norm=float(np.linalg.norm(g - mask * matrix)),
            )
            print(case, label, report["rows"][label], flush=True)
    report["spatial_nonspherical_induced_Ry"] = float((mask * spatial) @ induced)
    report["spherical_rows"] = {}
    spherical_mask = np.zeros_like(mask)
    for key, sl in zip(pc.hartree.keys, pc.hartree.slices):
        if key[1:] == (0, 0):
            spherical_mask[sl] = 1
    for name, bank, cutoff in [
        ("matrix", v.contexts, 10),
        ("matrix", v.contexts, 8),
        ("density", v.projector_contexts, 8),
    ]:
        sb = np.zeros(v.bands.input.spherical_size)
        for a, l, k, op in v.bands.spherical.spherical:
            if l > cutoff:
                continue
            ctx = bank[a][l]
            z = op.channel.value(ctx.primitive[-1], ctx.lo_coefficients, lapw=ctx.lapw)[
                ..., : v.model.eps.shape[1]
            ]
            d = np.einsum("min,np,mjp->ij", z, rho[k], z.conj(), optimize=True).real
            g = (
                np.einsum(
                    "rci,ij,rcj,rc->r",
                    ctx.primitive,
                    d,
                    ctx.primitive,
                    op.metric,
                    optimize=True,
                )
                / op.radius
            )
            sb[v.bands.spherical.slices[a]] += g
        native = v.bands.input.pack(
            sb,
            {a: np.zeros_like(op.v) for a, (op, ks) in v.bands.ns.items()},
            np.zeros(v.bands.input.nstar, complex if case == "SrVO3" else float),
        )
        g = select.rmatvec(native)
        difference = g - spherical_mask * spatial
        label = f"{name}_l{cutoff}"
        report["spherical_rows"][label] = dict(
            difference_norm=float(np.linalg.norm(difference)),
            induced_difference_Ry=float(difference @ induced),
            versus_native_matrix_norm=float(
                np.linalg.norm(g - spherical_mask * matrix)
            ),
        )
        print(case, "spherical", label, report["spherical_rows"][label], flush=True)
    (root / "three_k" / case / "angular_cutoff.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
