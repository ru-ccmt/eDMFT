"""Internal-coordinate derivative of the frozen-Sigma functional on a central LAPW state.

For an orbit displacement u the input is the same array of native MT radial
and interstitial Fourier density coefficients as in force_screened_diagnostic.
The atom-centred grids translate, RMT and local axes stay fixed. Differentiate
F + lambda . (n_out - n_in) at this *fixed* array and fixed central lambda.
Its negative derivative divided by orbit multiplicity is the force per atom.
One readout is shared by all modes. No displaced native run is required.
"""
import numpy as np
from scipy.special import spherical_jn

from force_response_modes import internal_modes
from lapw_charge_response import native_mixer_inversion
from lapw_atomic_geometry import atomic_geometry_readouts
from lapw_lattice import NativeLattice
from lapw_strain_field import diagonal_strain_potential
from lapw_strain_potential import masked_potential_strain_source


def atomic_charge_integral_derivative(density, structure, velocities):
    """Derivative of MIXER's IR integration weights; MT radial weights stay fixed.

    Native density PW entries are star sums. Integrating a translating sphere
    multiplies its transform by exp(-i G.R); its derivative enters the IR
    with the opposite sign. At G=0 the interstitial volume is unchanged.
    """
    labels = np.array(list(density.pw))
    q = NativeLattice.from_structure(structure).momenta(labels)
    g = np.linalg.norm(q, axis=1)
    mask = np.zeros(len(g), complex)
    for pos, radius, velocity in zip(structure.pos, np.repeat(structure.rmt, structure.mult), velocities):
        if not np.any(velocity):
            continue
        x = g * radius
        f = np.divide(spherical_jn(1, x), x, out=np.full_like(x, 1/3), where=x != 0)
        phase = np.exp(-2j*np.pi*(labels @ pos))
        mask += 1j * (q @ velocity) * phase * (4*np.pi*radius**3) * f
    if native_mixer_inversion(structure):
        mask.imag = 0
    return mask.view(float)


def atomic_force_report(total, readout, atoms, progress):
    """Evaluate both partials separately so the native +/- sources can check each.

    The band and spatial geometry terms are rank-local. Potential, normalization,
    radial energy searches and density energy are already globally reduced or
    replicated and enter only once.
    """
    state = total.state
    pc = state.potential_core
    s = pc.structure
    if hasattr(state.valence, "branches") or getattr(state.valence.model, "electronic_model", "dmft") != "dmft":
        raise NotImplementedError("analytic atomic sources require nonmagnetic frozen-Sigma eDMFT")
    modes = internal_modes(s, atoms)
    fixed_cell = np.zeros((3, 3))
    functional = total.functional
    density = functional.density
    electrostatic_weight = density.boundary_gradient - 0.5 * density.pairing.gradient(density.density)
    norm = state.normalization
    results = []
    force = np.zeros((s.nat, 3))
    local = atomic_geometry_readouts(state.valence, total.directory, total.case, modes,
                                      functional.lattice.joint_gradient, readout.band_covector,
                                      readout.raw_covector)
    fixed_geometry, induced_geometry = readout.parallel.sum(local)
    for i, (atom, direction, pattern) in enumerate(modes):
        progress(f"Analytic atomic source on atom {atom} along {direction}")
        potential, parts = diagonal_strain_potential(pc, fixed_cell, displacements=pattern)
        # RMT, cell volume and radial quadrature remain fixed. In MT coordinates
        # all position dependence of C_es is in the regular electrostatic field.
        # Both IR XC energy per electron and IR XC potential acquire the same
        # moving-mask derivative. Local XC fields at fixed coefficients stay fixed.
        ir = slice(pc.hartree.pw_offset, None)
        dxc_energy = masked_potential_strain_source(
            pc.hartree.reanalysis, density.xc.field_response.ir.central_stars, parts["mask_derivative"])
        dcoul = float(electrostatic_weight @ parts["coulomb_geometry"])
        dxc = pc.hartree.volume * float(density.density[ir] @ (dxc_energy.view(float) - parts["xc_geometry"][ir]))
        d_search = total.energy.A @ potential
        explicit = (dcoul + dxc + fixed_geometry[i]
                    + functional.potential_gradient @ potential + total.energy.g @ d_search)
        # The normalization rule integrates *raw output* charge with a moving
        # interstitial mask, even though the stored input array is unchanged.
        d_ir = atomic_charge_integral_derivative(pc.input_density, s, pattern)
        d_norm = norm.explicit_source(ir_geometry=float(d_ir @ norm.base[ir]))
        induced = (induced_geometry[i]
                   + readout.potential_covector @ potential
                   + readout.adjoint @ (d_norm + total.energy.P @ d_search))
        factor = -1000.0 / s.mult[atom-1]  # Ry -> mRy; per equivalent atom
        component = float(factor * (explicit + induced))
        results.append(dict(atom=atom, direction_cartesian=direction.tolist(),
                            explicit_force_mRy_per_Bohr=float(factor*explicit),
                            screened_force_mRy_per_Bohr=float(factor*induced),
                            force_component_mRy_Bohr=component))
        force[atom-1] += component * direction
        progress(f"Analytic screened force on atom {atom} along {direction}: {component:.9f} mRy/Bohr")
    return dict(observable="symmetry-preserving frozen-Sigma screened atomic force",
                scope="nonmagnetic frozen-Sigma atomic derivatives at fixed input charge, with one shared adjoint",
                selected_atoms=sorted(set(atoms or range(1, s.nat+1))),
                case=total.case, structure=str(total.directory / (total.case+".struct")),
                force_mRy_Bohr=force.tolist(), components=results, source_method="analytic")
