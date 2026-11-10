"""The implemented native LAPW0 density scalar and its discrete derivatives.

Nonmagnetic LDA/PBE, primitive orthogonal cells, fixed RMT, neutral native
input policy. All electrostatic constants and the native potential gauge are
retained. This is a functional component, not a separately additive stress.
"""

from types import SimpleNamespace
import numpy as np
from lapw_density_response import checked
from lapw_functional_cancellation import SpatialPotentialPairing, PreparedXCEnergy
from lapw_hydrostatic_potential import (
    hydrostatic_coulomb_geometry,
    hydrostatic_xc_geometry,
)
from wien_response_operators import hydrostatic_transport_source
from lapw_strain_field import diagonal_strain_potential


class PreparedDensityFunctional:
    """Differentiate WIEN ENERGY's DEN scalar, including nuclear constants.

    Write v_tilde=H*n+v_ion,boundary, i.e. omit the invariant local nuclear
    Dirichlet term -2Z(1-r/R) from r*V_c. Let b_a=v_tilde_a(R)/R. Then

      C_es=-n^T W v_tilde/2 - sum_a m_a Z_a^2/R_a
                                    - sum_a m_a Z_a b_a/2.

    This is ENERGY's VCOUL2 integral and ECONST combined algebraically;
    the singular -int n Z/r and +Z Q_MT/R terms cancel the omitted local
    Dirichlet contribution exactly. No self-adjointness of H is assumed.
    With t selecting -m_a Z_a/(2R_a) at the last spherical radial point,

      grad C_es=-W v_tilde/2 + H^T(-W n/2+t).

    C_xc and its gradient are the actual discrete XC scalar. Native radial
    integration, angular projection, Fourier truncation/mask and vacuum
    policy are inherited unchanged. Charged-cell background/E-field and
    nonlocal/meta-GGA energy terms require separate implementations.
    """

    def __init__(self, state):
        self.state = state
        h, s = state.hartree, state.structure
        self.pairing = SpatialPotentialPairing(h, s.mult)
        self.density = self.pairing.density_vector(state.input_density).copy()
        n = self.density
        self.xc = PreparedXCEnergy(state)
        _, vc, _ = hydrostatic_coulomb_geometry(state)
        self.coulomb = vc
        self.regular_coulomb = vc.copy()
        self.boundary_gradient = np.zeros(h.size)
        constant = 0.0
        for (a, l, m), sl in zip(h.keys, h.slices):
            if (l, m) != (0, 0):
                continue
            r = h.integrations[a].r
            z = s.Znuc[a]
            mult = s.mult[a]
            R = s.rmt[a]
            self.regular_coulomb[sl] += 2 * z * (1 - r / R)
            self.boundary_gradient[sl.stop - 1] = -mult * z / (2 * R)
            constant -= mult * z * z / R
        wn = self.pairing.gradient(n)
        self.electrostatic_value_ry = (
            -0.5 * self.pairing.value(n, self.regular_coulomb)
            + constant
            + self.boundary_gradient @ self.regular_coulomb
        )
        self.electrostatic_gradient = -0.5 * self.pairing.gradient(
            self.regular_coulomb
        ) + h.pullback(-0.5 * wn + self.boundary_gradient)
        self.value_ry = float(
            self.electrostatic_value_ry + self.xc.density_term_value_ry
        )
        self.density_gradient = (
            self.electrostatic_gradient + self.xc.density_term_gradient
        )
        self.variational_defect = (
            self.density_gradient + state.potential_operator.rmatvec(wn)
        )

    def action(self, density_direction):
        d = checked(
            density_direction, self.density.shape, "density derivative", real=True
        )
        return float(self.density_gradient @ d)

    def geometry(self, strain):
        """Partial derivative at fixed native density coefficients; RMT fixed.

        W' is zero in MT and 3W in IR. The Coulomb and XC field sources include
        moving masks and reciprocal vectors. Nuclear local Dirichlet terms
        and -Z^2/R are fixed. This is not yet the material-density convention.
        """
        pc = self.state
        h = pc.hartree
        n = self.density
        _, parts = diagonal_strain_potential(pc, strain, transport_rate=0)
        dc, dv = parts["coulomb_geometry"], parts["xc_geometry"]
        energy_state = SimpleNamespace(
            hartree=h,
            xc=self.xc.field_response,
            input_density=pc.input_density,
            structure=pc.structure,
            potential_operator=pc.potential_operator,
        )
        _, energy_parts = diagonal_strain_potential(
            energy_state, strain, transport_rate=0
        )
        de = energy_parts["xc_geometry"]
        trace = np.trace(strain)
        ir = slice(h.pw_offset, None)
        es = -0.5 * self.pairing.value(n, dc) + self.boundary_gradient @ dc
        es -= 0.5 * trace * h.volume * (n[ir] @ self.regular_coulomb[ir])
        xc = self.pairing.value(n, de - dv)
        xc += (
            trace
            * h.volume
            * (n[ir] @ (self.xc.energy_field - self.xc.potential_field)[ir])
        )
        return float(es + xc), dict(electrostatic=float(es), xc=float(xc))

    def explicit(self, strain, *, transport_rate=None):
        """Geometry plus the same scalar-density transport used in (I-J)q=s."""
        geometry, parts = self.geometry(strain)
        rate = np.trace(strain) / 3 if transport_rate is None else float(transport_rate)
        transport = rate * self.action(
            hydrostatic_transport_source(self.state.input_density).vector()
        )
        return geometry + transport, parts | dict(
            transport=transport, geometry=geometry
        )

    def hydrostatic_explicit(self):
        return self.explicit(np.eye(3))
