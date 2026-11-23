"""Diagnostic pullbacks to native radial linearization energies.

These are covectors, not a prescription for the strain derivative of an
energy-search algorithm. They expose the energy inputs discarded under the
fixed-energy policy; both native matrix and ATPARN projector/density banks
are retained. All coefficients are Ry per Ry of the corresponding energy.
"""

def band_energy_gradient(spherical, vector):
    """Energy covectors of the compact spherical bands in their two radial banks."""
    return spherical.pullback(vector, with_energies=True)[1]


def density_energy_gradient(direct, covector):
    """Energy covectors of the moving-basis density at fixed band density."""
    return direct.pullback(covector, with_energies=True)[1]


def screened_energy_gradient(state, observable, adjoint):
    """dGamma/dE, including the induced self-consistent density via one lambda.

    E denotes the independent native radial-energy inputs here. The result
    contracts the TOTAL E_nu from any separately specified energy-search
    policy. It must not be added if that policy already entered B and J.
    """
    ch = state.chain
    normalized = ch.normalization.rmatvec(adjoint)
    joint = observable.lattice.joint_gradient + ch.electronic.rmatvec(
        ch.density.rmatvec(normalized)
    )
    ge = band_energy_gradient(state.valence.bands.spherical, joint)
    for key, value in density_energy_gradient(state.valence.direct, normalized).items():
        ge[key] += value
    return ge
