"""Complete-basis readout of the extra normalized-projector lattice term.

This differentiates the regularized lattice scalar ONLY through its normalized
projector U=Pi_W phi X. It excludes the direct density/Hamiltonian derivative,
impurity/DC functional terms and the already assembled baseline stress. It is
therefore a building block for the physical stress correction, not a total
Klein-functional or stress evaluator.

All covectors use the same real isometric coordinates as the screened charge
operator. K/spin/frequency weights and normalization symmetry are already in
the retained projector gradient and must not be applied again here.
"""

from __future__ import annotations
import numpy as np
from frozen_sigma_response import herm
from lapw_density_response import checked
from lapw_complement_response import (
    ComplementMatrixCoordinates,
    ComplementPerturbationCoordinates,
)


def projector_functional_gradient(complement):
    """Return (packed covector, coordinates) for joint H/O/raw-R inputs.

    Let (h_bar,p_bar) be the retained orthonormal projector-only covectors.
    The moving metric gives
      h_eff = H_RR' - {O_RR',epsilon}/2,
      p_eff = phi_R' - O_RR' phi_R/2 + T^dagger phi_Q,
      H_Q T - T epsilon = -(H_QR' - O_QR' epsilon).
    Reverse the SAME prepared shifted solves used by the charge response:
      T_bar = phi_Q p_bar^dagger,
      Z = shifted_solve(T_bar),
      H_bar_QR = Z/2, O_bar_QR = -Z epsilon/2.
    The halves follow the Hermitian RR+QR+RQ pairing. Raw overlaps use an
    ordinary complex Frobenius pairing and have no such half.

    There is no density-metric readout term -(y rho+rho y)/2: that belongs to
    the derivative of density coefficients, whereas this is a scalar readout.
    The thermally empty, gapped complement policy is inherited from the
    prepared EmptyComplementResponse. No high-state eigensolve or frequency
    sum is performed after the retained gradient has been prepared.
    """
    model = complement.model
    hbar, pbar = model.projector_gradients()
    obar = -0.5 * (hbar * model.eps[:, None, :] + model.eps[:, :, None] * hbar)
    obar -= 0.5 * herm(np.einsum("knic,kmic->knm", pbar, model.raw.conj()))
    hq, oq = [], []
    for k, solver in enumerate(complement.solvers):
        tb = (
            complement.raw_q[k].reshape(solver.nq, -1)
            @ pbar[k].reshape(solver.nr, -1).conj().T
        )
        z = solver.action(tb)
        hq.append(z / 2)
        oq.append(-z * solver.energies / 2)
    nk, nr = model.eps.shape
    matrix = ComplementMatrixCoordinates(
        (nk, nr, nr), [s.nq for s in complement.solvers]
    )
    coordinates = ComplementPerturbationCoordinates(matrix, model.raw.shape)
    gradient = coordinates.pack(hbar, herm(obar), pbar, hq, oq)
    gradient.setflags(write=False)
    return gradient, coordinates


class PreparedProjectorFunctional:
    """One common projector-only potential/charge gradient for many sources.

    Construct from NativeScreenedState. Its prepared native basis operator B
    includes H/O/raw-projector radial, nonspherical and interstitial response
    and the native coordinate selector. Then
      g_v = B^T g_joint,  c = K^T g_v.
    These two pullbacks occur ONCE. Density-response iterations use the same J
    as the screened calculation. Keep the owning valence contexts alive.

    This does not decide which pieces overlap a production stress baseline.
    Every returned report identifies the scalar and keeps stress_applied=False.
    """

    def __init__(self, screened_state):
        self.state = screened_state
        self.joint_gradient, self.coordinates = projector_functional_gradient(
            screened_state.valence.complement
        )
        basis = screened_state.chain.basis
        if basis.output_size != self.coordinates.size:
            raise ValueError("projector and native basis coordinates differ")
        self.potential_gradient = basis.rmatvec(self.joint_gradient)
        self.charge_gradient = screened_state.chain.potential.rmatvec(
            self.potential_gradient
        )
        self.potential_gradient.setflags(write=False)
        self.charge_gradient.setflags(write=False)

    def potential_contraction(self, potential):
        """Projector-only derivative in Ry for a TOTAL supplied potential.

        The caller supplies v_explicit+Kq in prepared native field coordinates.
        Explicit geometry changes of H/O/phi at fixed potential are separate;
        this operation cannot manufacture or silently add them.
        """
        v = checked(
            potential, self.potential_gradient.shape, "potential direction", real=True
        )
        return float(self.potential_gradient @ v)

    def adjoint_contractions(self, sources, explicit_potentials, **controls):
        """One adjoint: g_v^T v_explicit,nu + lambda^T s_nu for each source.

        Sources have shape (density coordinates, directions) and include the
        complete explicit NORMALIZED charge source. explicit_potentials has
        shape (prepared potential coordinates, directions); it is required so
        omission of the affine/ionic source cannot silently mean zero.
        Output values are energy derivatives in Ry, without volume division.
        This method returns only the potential-dependent projector contribution,
        not the full a_nu or stress tensor. Other explicit band geometry and
        baseline terms must be assembled separately with no double counting.
        """
        sources = np.asarray(sources)
        if sources.ndim != 2 or sources.shape[1] < 1:
            raise ValueError("sources require one or more direction columns")
        sources = checked(
            sources,
            (self.state.jacobian.input_size, sources.shape[1]),
            "density sources",
            real=True,
        )
        potentials = checked(
            explicit_potentials,
            (len(self.potential_gradient), sources.shape[1]),
            "explicit potential sources",
            real=True,
        )
        explicit = self.potential_gradient @ potentials
        values, lam, report = self.state.screened.contractions(
            self.charge_gradient, sources, explicit=explicit, **controls
        )
        report = dict(
            report,
            observable="extra normalized-projector regularized lattice term",
            units="Ry per supplied perturbation parameter",
            directions=sources.shape[1],
            explicit_derivatives=explicit.tolist(),
            induced_derivatives=(values - explicit).tolist(),
            stress_applied=False,
            includes_direct_density_hamiltonian_term=False,
        )
        return values, lam, report
