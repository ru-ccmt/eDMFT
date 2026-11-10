"""Complete frozen-self-energy Klein tangent in the prepared native coordinates.

The scalar differentiated here is C_DEN[n,cell]+sum_core f E+L[H,O,U].
The impurity/DC constants are fixed. Its derivative is the Klein tangent at
that fixed local self-energy; it is not the derivative of SUM or the DMFT
internal energy. All three pieces must be combined before assessing stress.

This module is a readout, not a production publisher. It assumes the fixed
RMT/radial-energy/frame and nonmagnetic policies of NativeScreenedState.
Actual fixed-point convergence and agreement of radial-energy policies with
finite-difference reference runs remain caller obligations.
"""

import numpy as np
from lapw_density_response import checked
from lapw_strain_functional import PreparedStrainFunctional
from lapw_functional_cancellation import CoreEigenvalueFunctional
from lapw_density_functional import PreparedDensityFunctional


class PreparedKleinFunctional:
    """One scalar and one adjoint for all supported strain directions.

    Write b_nu=b_geo,nu+B v_exp,nu, K=partial V/partial n, and
      (I-J) q_nu=s_nu.
    The complete derivative is
      D_nu=a_nu+c.q_nu,
      a_nu=g_L.b_geo,nu+(B^T g_L+g_core).v_exp,nu+C_exp,nu,
      c=K^T(B^T g_L+g_core)+grad_n C_DEN.
    Equivalently solve (I-J)^T lambda=c once and form a_nu+lambda.s_nu.

    The complete discrete C_DEN gradient is used, not the continuum
    cancellation -K^T W n. That replacement can miss quadrature, XC, nuclear
    boundary and gauge terms. No standalone M_Sigma overlap term is appended:
    the normalized projector/metric chain already occurs in g_L.
    """

    def __init__(self, state, *, lattice=None):
        self.state = state
        self.lattice = PreparedStrainFunctional(state) if lattice is None else lattice
        if self.lattice.state is not state or self.lattice.component != "lattice":
            raise ValueError("the complete lattice covector of this state is required")
        self.core = CoreEigenvalueFunctional.from_state(state.potential_core)
        self.density = PreparedDensityFunctional(state.potential_core)
        self.potential_gradient = (
            self.lattice.potential_gradient + self.core.potential_gradient
        )
        self.charge_gradient = state.observable_gradient(
            self.density.density_gradient, self.potential_gradient
        )
        self.potential_gradient.setflags(write=False)
        self.charge_gradient.setflags(write=False)

    def explicit(self, geometry, potentials, strains):
        geo, pot_lattice = self.lattice._explicit(geometry, potentials)
        strains = checked(strains, (len(geo), 3, 3), "strain directions", real=True)
        core = self.core.potential_gradient @ potentials
        density = np.array([self.density.explicit(e)[0] for e in strains])
        return geo + pot_lattice + core + density, dict(
            lattice_geometry_Ry=geo.tolist(),
            lattice_explicit_potential_Ry=pot_lattice.tolist(),
            core_explicit_Ry=core.tolist(),
            density_explicit_Ry=density.tolist(),
        )

    def _report(self, explicit, induced, parts):
        return dict(
            observable="complete frozen-Sigma native Klein tangent",
            units="Ry per supplied perturbation parameter",
            explicit_Ry=explicit.tolist(),
            induced_Ry=induced.tolist(),
            parts=parts,
            full_klein_tangent=True,
            stress_applied=False,
            baseline_subtracted=False,
            radial_energy_policy="fixed",
            central=dict(self.state.metadata),
        )

    def forward_contractions(self, geometry, potentials, strains, density_responses):
        explicit, parts = self.explicit(geometry, potentials, strains)
        q = checked(
            density_responses,
            (len(self.charge_gradient), len(explicit)),
            "screened density responses",
            real=True,
        )
        induced = self.charge_gradient @ q
        return explicit + induced, self._report(explicit, induced, parts)

    def adjoint_contractions(self, geometry, potentials, strains, sources, **controls):
        explicit, parts = self.explicit(geometry, potentials, strains)
        sources = checked(
            sources,
            (len(self.charge_gradient), len(explicit)),
            "normalized charge sources",
            real=True,
        )
        values, lam, solver = self.state.screened.contractions(
            self.charge_gradient, sources, explicit=explicit, **controls
        )
        report = self._report(explicit, lam @ sources, parts)
        report["adjoint"] = solver
        return values, lam, report
