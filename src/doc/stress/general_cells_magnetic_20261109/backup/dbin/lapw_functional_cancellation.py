"""Core and spatial potential readouts for auditing Klein cancellations.

These are fixed-geometry scalar building blocks, not an additive stress.
In continuum notation the density part C=E_Hxc-int(n V_Hxc)+E_ion has
dC=-int(n dV_Hxc) at fixed geometry. The discrete LAPW0 energy need not
satisfy that identity exactly. We therefore expose the spatial pairing and
the core eigenvalue derivative separately; neither is silently replaced by
a band-matrix trace. The input charge, raw output and normalized output are
different objects until convergence and must be supplied deliberately.

The current periodic adapter supports the same primitive orthogonal cells,
uncompressed harmonics, nonmagnetic LDA/PBE and fixed-RMT policy as the native
potential/core state. No density response, eigenvalue search or native process
is run when these covectors are applied.
"""

from __future__ import annotations

import numpy as np

from lapw_density_response import checked
from lapw_hartree_periodic import native_charge_weights
from stress_response_adjoint import PairedOperator


class SpatialPotentialPairing:
    """Native int(n dV) at fixed geometry and its two exact derivatives.

    Density coordinates contain 4*pi*r^2*n_sph, r^2*n_lm and complex PW
    STAR SUMS. Potential coordinates contain r*V_sph (Ry), V_lm (Ry) and
    already interstitial-masked per-member Fourier coefficients. Thus the
    real-coordinate pairing is n^T W v, with W diagonal:

      W_sph = multiplicity * CHARGE_weight / r,
      W_lm  = multiplicity * CHARGE_weight,
      W_PW  = cell volume (on both real and imaginary coordinates).

    No additional PW star multiplicity or muffin-tin mask belongs here.
    Cubic compressed-harmonic conventions are rejected by the layout builder.
    This map does not normalize charge or truncate it to the H angular cutoff.
    """

    def __init__(self, layout, multiplicities):
        mult = checked(
            multiplicities, (len(layout.integrations),), "multiplicities", real=True
        )
        if np.any(mult < 1) or np.any(mult != np.floor(mult)):
            raise ValueError("positive integer atom multiplicities required")
        self.size = layout.size
        self.keys = layout.keys
        self.pw_keys = tuple(map(tuple, layout.stars.representatives))
        self.grids = tuple(np.log(integral.r) for integral in layout.integrations)
        self.weights = np.full(self.size, layout.volume)
        for (a, l, m), sl in zip(layout.keys, layout.slices):
            r = layout.integrations[a].r
            self.weights[sl] = mult[a] * native_charge_weights(r)
            if (l, m) == (0, 0):
                self.weights[sl] /= r
        self.weights.setflags(write=False)
        self.operator = PairedOperator(
            self.size,
            self.size,
            self.gradient,
            self.gradient,
            "fixed-geometry native spatial density/potential pairing",
        )

    def density_vector(self, density):
        """Reject reordered/missing channels instead of silently pairing them."""
        if tuple(density.rad) != self.keys or tuple(density.pw) != self.pw_keys:
            raise ValueError("density coordinates differ from potential layout")
        if len(density.grids) != len(self.grids) or any(
            np.shape(actual) != expected.shape
            or not np.allclose(actual, expected, atol=1e-12, rtol=0)
            for actual, expected in zip(density.grids, self.grids)
        ):
            raise ValueError("density radial grid differs from potential layout")
        return checked(density.vector(), (self.size,), "density", real=True)

    def gradient(self, vector):
        """W times density (potential gradient), or potential (density gradient)."""
        return self.weights * checked(vector, (self.size,), "pairing input", real=True)

    def value(self, density, potential):
        """Bilinear scalar in Ry; inputs are explicit native real vectors."""
        potential = checked(potential, (self.size,), "potential", real=True)
        return float(self.gradient(density) @ potential)


class CoreEigenvalueFunctional:
    """Occupation-weighted native core sum in Ry, with its potential covector.

    LCORE matches F(E,V)=0, so dE=-F_V dV/F_E. MatchedCoreState supplies
    exactly that discrete derivative, including the origin and outer boundary.
    Its energies/potentials are Hartree, whereas the periodic input is rV in
    Ry. For atom a, with extension operator T_a and fixed occupations f_ai,

      E_c = sum_ai m_a f_ai (2 E_ai + shift_a),
      g_c[a] = 2 m_a T_a^T sum_i f_ai grad_V E_ai.

    The factor 2 converts Hartree to Ry; it is NOT a spin factor. Occupations
    already include the nonmagnetic degeneracy. The fixed inc shift has zero
    derivative. Multiplicity is included here, not inside atom core density.
    Native MIXER scfana.f uses precisely these occupations/multiplicities.

    T includes the extrapolated/clamped potential outside RMT. Omitting it
    loses a genuine derivative even though CLMCOR is written only inside RMT.
    interior_gradient/exterior_gradient split this effect diagnostically.
    This readout does not replace the explicit geometric core stress (E3).
    """

    def __init__(self, core, multiplicities, shifts_ry):
        nat = len(core.layout.integrations)
        mult = checked(multiplicities, (nat,), "multiplicities", real=True)
        shifts = checked(shifts_ry, (nat,), "fixed core shifts", real=True)
        if np.any(mult < 1) or np.any(mult != np.floor(mult)):
            raise ValueError("positive integer atom multiplicities required")
        self.value_ry = 0.0
        self.potential_gradient = np.zeros(core.layout.size)
        self.interior_gradient = np.zeros_like(self.potential_gradient)
        self.core_number = 0.0
        for a, atom in core.atoms.items():
            bar = np.zeros(len(atom.extension.r))
            for state in atom.states:
                self.value_ry += mult[a] * state.occ * (2 * state.e + shifts[a])
                self.core_number += mult[a] * state.occ
                bar += 2 * mult[a] * state.occ * state.energy_gradient
            sl = core.slices[a]
            self.potential_gradient[sl] = atom.extension.pullback(bar)
            interior = bar.copy()
            interior[atom.n :] = 0
            self.interior_gradient[sl] = atom.extension.pullback(interior)
        self.exterior_gradient = self.potential_gradient - self.interior_gradient
        for g in [
            self.potential_gradient,
            self.interior_gradient,
            self.exterior_gradient,
        ]:
            g.setflags(write=False)

    @classmethod
    def from_state(cls, state):
        return cls(
            state.core,
            state.structure.mult,
            [atom["shift_Ry"] for atom in state.metadata["core"]],
        )

    def action(self, potential):
        return float(
            self.potential_gradient
            @ checked(
                potential,
                self.potential_gradient.shape,
                "potential response",
                real=True,
            )
        )


def stationary_potential_readout(state, lattice_gradient):
    """Return the candidate induced-potential cancellation, explicitly scoped.

    g = g_lattice + g_core - W*n_INPUT; c = K^T*g.

    This is the continuum stationary-density identity, evaluated with native
    spatial coordinates. It is not an assertion that native LAPW0's scalar
    satisfies that identity. A separate native-energy check must measure its
    discrete defect. Explicit geometric derivatives of the density functional
    and subtraction of the existing Fortran baseline also remain downstream.
    Use n_INPUT, not an output density or P^T*rho, for the LAPW0 subtraction.
    No normalization pullback belongs here: normalization already enters q.
    """
    pairing = SpatialPotentialPairing(state.hartree, state.structure.mult)
    core = CoreEigenvalueFunctional.from_state(state)
    lattice = checked(lattice_gradient, (pairing.size,), "lattice gradient", real=True)
    spatial = pairing.gradient(pairing.density_vector(state.input_density))
    potential = lattice + core.potential_gradient - spatial
    return dict(
        lattice=lattice.copy(),
        core=core.potential_gradient.copy(),
        input_density_pairing=spatial,
        potential=potential,
        charge=state.potential_operator.rmatvec(potential),
        full_klein_functional=False,
        native_density_functional_validated=False,
        baseline_subtracted=False,
        stress_applied=False,
    )


class PreparedXCEnergy:
    """Differentiate the implemented spatial XC scalar, not its continuum limit.

    Let e[n] be native XC energy PER ELECTRON (Ry) projected/masked into the
    same field coordinates as V_xc, and W the spatial pairing. Then

      E_xc[n] = n^T W e[n],
      g_E = W e + (D e)^T W n,
      d(E_xc - int n V_xc)/dn = g_E - W V_xc - K_xc^T W n.

    The continuum cancellation g_E = W V_xc need not hold after native radial
    splines, angular projection, Fourier truncation and gradient stencils.
    variational_defect is this discrepancy in density coordinates. It must
    not be thrown away when differentiating the actually implemented scalar.

    Construction reuses geometry; complex-step local derivatives are prepared
    once. Subsequent scalar contractions need only a dot product. This is
    the XC portion at fixed geometry, not the full LAPW0 density functional.
    """

    def __init__(self, state):
        from lapw_xc_response import (
            MuffinTinXCResponse,
            InterstitialXCResponse,
            PeriodicXCResponse,
        )

        layout = state.hartree
        pair = SpatialPotentialPairing(layout, state.structure.mult)
        density = pair.density_vector(state.input_density)
        mt = []
        for a, source in enumerate(state.xc.mt):
            central = np.column_stack(
                [state.input_density.rad[k] for k in layout.keys if k[0] == a]
            )
            extra = (
                {}
                if source.functional == "LDA"
                else dict(
                    sine=source.sine.ravel(),
                    cosine=source.cosine.ravel(),
                    angular_derivatives=source.angular[1:],
                )
            )
            mt.append(
                MuffinTinXCResponse(
                    source.r,
                    source.ells,
                    central,
                    source.Y,
                    source.weights,
                    functional=source.functional,
                    harmonics=source.harmonics,
                    quantity="energy_per_electron",
                    **extra,
                )
            )
        ir = InterstitialXCResponse(
            layout.reanalysis,
            np.array(list(state.input_density.pw.values())),
            [state.structure.a, state.structure.b, state.structure.c],
            functional=state.metadata["functional"],
            rho_min=state.metadata["rho_min"],
            quantity="energy_per_electron",
        )
        energy_map = PeriodicXCResponse(layout, mt, ir)
        self.field_response = energy_map
        self.energy_operator = energy_map.operator
        self.energy_field = np.zeros(layout.size)
        self.potential_field = np.zeros(layout.size)
        for i, slices in enumerate(energy_map.atom_slices):
            for j, sl in enumerate(slices):
                self.energy_field[sl] = mt[i].central[:, j]
                self.potential_field[sl] = state.xc.mt[i].central[:, j]
        self.energy_field[layout.pw_offset :] = ir.central.view(float)
        self.potential_field[layout.pw_offset :] = state.xc.ir.central.view(float)
        wn = pair.gradient(density)
        self.value_ry = pair.value(density, self.energy_field)
        self.density_term_value_ry = pair.value(
            density, self.energy_field - self.potential_field
        )
        self.energy_gradient = pair.gradient(
            self.energy_field
        ) + self.energy_operator.rmatvec(wn)
        self.variational_defect = self.energy_gradient - pair.gradient(
            self.potential_field
        )
        self.density_term_gradient = (
            self.variational_defect - state.xc.operator.rmatvec(wn)
        )
        for field in [
            self.energy_field,
            self.potential_field,
            self.energy_gradient,
            self.variational_defect,
            self.density_term_gradient,
        ]:
            field.setflags(write=False)
