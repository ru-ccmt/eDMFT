"""Strain readout of the regularized frozen-Sigma LATTICE scalar.

This is a diagnostic building block, not the full Klein functional and not an
additive correction to production stress. It keeps the direct Hamiltonian and
normalized-projector terms separate. Potential double counting, core, DFT,
impurity/DC and existing analytic-stress terms must be reconciled downstream.

All contractions use real isometric RR/QR/raw coordinates. K-point, spin and
Matsubara weights are already present in the electronic covectors. The empty
complement and supported material policies are inherited from the state.
"""

from __future__ import annotations

import numpy as np

from lapw_density_response import checked
from lapw_projector_functional import projector_functional_gradient


def lattice_functional_gradients(complement):
    """Return separated joint covectors and their common coordinate system.

    In the symmetric orthonormal gauge, at the central KS eigenbasis,

      dL = Re Tr[rho (H' - {O',epsilon}/2)] + <g_U,dU>.

    The installed projector pullback differentiates U=Pi_W phi X, INCLUDING
    the metric term -O' phi/2, window motion, overlap symmetrization and X'.
    Consequently the direct contribution here is precisely

      g_H = rho,  g_O = -{rho,epsilon}/2,  g_phi = g_H,QR = g_O,QR = 0.

    Q is thermally empty and outside the DMFT window; its central density and
    density cross blocks vanish under the prepared complement policy. Its
    projector response is nevertheless retained in the projector covector.

    Do NOT append an isolated self-energy-weighted overlap moment here. Such
    a term can occur in a rearranged covariant formula, but the corresponding
    projector-metric contribution must then be rearranged too. The complete
    expression below is tested under a general moving-basis coordinate change.
    A density-response output metric term belongs to d(rho), not to this scalar.
    """
    projector, coordinates = projector_functional_gradient(complement)
    model = complement.model
    rho = model.central_density()
    overlap = -0.5 * (rho * model.eps[:, None, :] + model.eps[:, :, None] * rho)
    zero_qr = [np.zeros((s.nq, s.nr), complex) for s in complement.solvers]
    direct = coordinates.pack(rho, overlap, np.zeros_like(model.raw), zero_qr, zero_qr)
    gradients = dict(projector=projector, direct=direct, lattice=projector + direct)
    for gradient in gradients.values():
        gradient.setflags(write=False)
    return gradients, coordinates


class PreparedStrainFunctional:
    """One scalar covector and one adjoint, reusable for many strain sources.

    Required explicit inputs prevent accidentally omitting strain geometry:
      b_nu = b_geometry,nu + B v_explicit,nu,
      (I-J)q_nu = s_nu,
      D_nu = g_b b_nu + c q_nu,  c = K^T B^T g_b.
    Equivalently (I-J)^T lambda=c gives D_nu=g_b b_nu+lambda^T s_nu.

    Select component='projector', 'direct', or 'lattice'. Only 'lattice' is
    invariant under a general basis change; the two separate pieces may cancel.
    Neither component nor their sum is the full stress. Keep native contexts
    alive for the lifetime of this object. No new native process is launched.
    """

    def __init__(self, state, *, component="lattice", gradients=None):
        if component not in {"projector", "direct", "lattice"}:
            raise ValueError("component must be projector, direct, or lattice")
        self.state = state
        self.component = component
        if gradients is None:
            gradients, self.coordinates = lattice_functional_gradients(
                state.valence.complement
            )
        else:
            # An internal shared preparation avoids repeating frequency work
            # when separate direct/projector reports are requested together.
            gradients, self.coordinates = gradients
        self.joint_gradient = gradients[component]
        if state.chain.basis.output_size != self.coordinates.size:
            raise ValueError("functional and basis coordinates differ")
        self.potential_gradient = state.chain.basis.rmatvec(self.joint_gradient)
        self.charge_gradient = state.chain.potential.rmatvec(self.potential_gradient)
        self.potential_gradient.setflags(write=False)
        self.charge_gradient.setflags(write=False)

    def _explicit(self, geometry, potentials):
        geometry = np.asarray(geometry)
        if geometry.ndim != 2 or geometry.shape[1] < 1:
            raise ValueError("geometry requires one or more direction columns")
        geometry = checked(
            geometry,
            (self.coordinates.size, geometry.shape[1]),
            "explicit band geometry",
            real=True,
        )
        potentials = checked(
            potentials,
            (len(self.potential_gradient), geometry.shape[1]),
            "explicit potential sources",
            real=True,
        )
        return self.joint_gradient @ geometry, self.potential_gradient @ potentials

    def _report(self, geometry, potential, induced):
        return dict(
            observable=self.component + " regularized frozen-Sigma lattice term",
            units="Ry per supplied perturbation parameter",
            geometry_derivatives=geometry.tolist(),
            explicit_potential_derivatives=potential.tolist(),
            induced_derivatives=induced.tolist(),
            stress_applied=False,
            full_klein_functional=False,
            baseline_subtracted=False,
        )

    def forward_contractions(self, geometry, explicit_potentials, density_responses):
        """Contract already solved q, including both explicit strain sources."""
        geo, pot = self._explicit(geometry, explicit_potentials)
        q = checked(
            density_responses,
            (len(self.charge_gradient), len(geo)),
            "screened density responses",
            real=True,
        )
        induced = self.charge_gradient @ q
        return geo + pot + induced, self._report(geo, pot, induced)

    def adjoint_contractions(self, geometry, explicit_potentials, sources, **controls):
        """Solve one adjoint and contract every supplied complete charge source."""
        geo, pot = self._explicit(geometry, explicit_potentials)
        sources = checked(
            sources,
            (len(self.charge_gradient), len(geo)),
            "normalized charge sources",
            real=True,
        )
        values, lam, solver = self.state.screened.contractions(
            self.charge_gradient, sources, explicit=geo + pot, **controls
        )
        report = self._report(geo, pot, lam @ sources)
        report["adjoint"] = solver
        return values, lam, report
