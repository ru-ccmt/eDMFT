"""Material checks of the matrix-defined fixed-basis potential subtraction."""

from pathlib import Path
import json
import sys
import numpy as np
from check_material import root
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from lapw_physical_potential import FixedBasisPhysicalPotential
from wien_response_density import Density

case = sys.argv[1]
source = Path("/private/tmp/edmft_full_k_response_20261108") / case / "charge/0"
with load_native_valence_state(
    source,
    case,
    complex_vectors=(case == "SrVO3"),
    energy_policy="fixed",
    matrix_quadrature="wien2k24",
    library="/Users/haule/dbin/liblapw_response.so",
) as v:
    pc = load_potential_core_state(
        "/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/"
        + case
        + "_potential_core.npz"
    )
    core = Density(
        Path("/private/tmp/edmft_xc_core_response_20261108")
        / case
        / "hxc_core/0"
        / (case + ".clmcor"),
        source / (case + ".struct"),
    )
    state = NativeScreenedState(
        v, pc, core, target_charge={"SrVO3": 85, "FeSe": 120}[case], fixed_number=True
    )
    p = FixedBasisPhysicalPotential(v.bands)
    operator = p.operator @ state.potential_selector
    rho = v.model.central_density()
    native_gradient = p.density_gradient(rho)
    field_gradient = state.potential_selector.rmatvec(native_gradient)
    direction = np.load(
        "/private/tmp/edmft_hartree_periodic_20261108/" + case + "_direction.npz"
    )["density"]
    potential = pc.potential_operator.matvec(direction)
    h, o, raw, hq, oq = p.output.unpack(operator.matvec(potential))
    analytic = float(field_gradient @ potential)
    direct = float(np.vdot(rho, h).real)
    assert abs(analytic - direct) < 1e-10
    report = dict(
        case=case,
        kpoints=len(v.records),
        analytic_Ry=analytic,
        direct_Ry=direct,
        steps=[],
        stress_applied=False,
        scope="matrix-defined fixed-basis physical-potential scalar, not full Klein cancellation",
    )
    for step in [0.001, 0.0005]:
        # Hold normalized U fixed; ONLY the physical H changes. This is the
        # ordinary potential term, so window/normalization must not be added.
        values = [v.model.lattice_scalar(h, sign * step) for sign in [-1, 1]]
        fd = (values[1] - values[0]) / (2 * step)
        report["steps"].append(dict(step=step, fd_Ry=fd, error_Ry=fd - analytic))
    fd = (4 * report["steps"][1]["fd_Ry"] - report["steps"][0]["fd_Ry"]) / 3
    report["richardson"] = dict(
        fd_Ry=fd,
        error_Ry=fd - analytic,
        relative_error=abs(fd - analytic) / max(abs(analytic), 1e-30),
    )
    rng = np.random.default_rng(826)
    y = rng.normal(size=p.output.size)
    lhs = float(y @ operator.matvec(potential))
    rhs = float(operator.rmatvec(y) @ potential)
    report["transpose_relative_error"] = abs(lhs - rhs) / max(1, abs(lhs), abs(rhs))
    field = np.zeros(pc.hartree.size)
    for key, sl in zip(pc.hartree.keys, pc.hartree.slices):
        if key[1:] == (0, 0):
            field[sl] = pc.hartree.integrations[key[0]].r
    constant = np.zeros(pc.hartree.stars.nstar, complex)
    constant[0] = 1
    field[pc.hartree.pw_offset :] = pc.hartree.reanalysis.action(constant).view(float)
    hh, _, _, hq, _ = p.output.unpack(operator.matvec(field))
    report["constant_potential_H_RR_minus_I_norm"] = float(
        np.linalg.norm(hh - np.eye(hh.shape[1]))
    )
    report["constant_potential_H_QR_norm"] = float(
        np.sqrt(sum(np.linalg.norm(a) ** 2 for a in hq))
    )
    report["number"] = float(np.trace(rho, axis1=1, axis2=2).real.sum())
    report["constant_potential_number_derivative"] = float(field_gradient @ field)
    old = np.load(root / "three_k" / case / "covectors.npz")
    remainder = old["lattice_potential"] - field_gradient
    radial = old["lattice_potential"] - old["projector_potential"] - field_gradient
    response = np.load(
        Path("/private/tmp/edmft_strain_charge_20261109/three_k")
        / case
        / "analytic_induced_potential.npy"
    )
    report["induced_potential_split_Ry"] = dict(
        lattice=float(old["lattice_potential"] @ response),
        ordinary_physical=float(field_gradient @ response),
        radial_basis=float(radial @ response),
        projector=float(old["projector_potential"] @ response),
        matrix_remainder=float(remainder @ response),
    )
    out = root / "three_k" / case
    np.savez_compressed(
        out / "physical_potential_covectors.npz",
        physical=field_gradient,
        radial=radial,
        matrix_remainder=remainder,
    )
    (out / "physical_potential.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)
