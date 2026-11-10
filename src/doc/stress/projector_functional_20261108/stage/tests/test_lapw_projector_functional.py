"""Independent nonlinear scalar, metric/complement and screened readout checks."""

from dataclasses import replace
from types import SimpleNamespace as NS
import unittest
import numpy as np
from frozen_sigma_response import herm
from frozen_sigma_tangent import CachedFrozenSigmaResponse, real_inner
from lapw_complement_response import EmptyComplementResponse
from lapw_projector_functional import (
    projector_functional_gradient,
    PreparedProjectorFunctional,
)
from stress_response_adjoint import PairedOperator, ScreenedResponse
from test_lapw_complement_response import fixture
from test_generalized_sigma_response import endpoint


def dense(a):
    return PairedOperator(a.shape[1], a.shape[0], lambda x: a @ x, lambda y: a.T @ y)


class ProjectorFunctionalTests(unittest.TestCase):
    def setUp(self):
        self.small, self.full, self.m, self.t, self.rng = fixture()
        self.gradient, self.coordinates = projector_functional_gradient(self.m)
        self.nr = self.small.eps.shape[1]

    def perturbation(self):
        n = self.full.eps.shape[1]
        rng = self.rng
        h = 0.15 * herm(rng.normal(size=(2, n, n)) + 1j * rng.normal(size=(2, n, n)))
        o = 0.07 * herm(rng.normal(size=h.shape) + 1j * rng.normal(size=h.shape))
        p = 0.07 * (
            rng.normal(size=self.full.raw.shape)
            + 1j * rng.normal(size=self.full.raw.shape)
        )
        return h, o, p

    def packed(self, h, o, p):
        td = self.t.swapaxes(1, 2).conj()
        h = td @ h @ self.t
        o = td @ o @ self.t
        p = np.einsum("knm,kmic->knic", td, p)
        n = self.nr
        return self.coordinates.pack(
            h[:, :n, :n], o[:, :n, :n], p[:, :n], h[:, n:, :n], o[:, n:, :n]
        )

    def scalar(self, h, o, p, t):
        # Full independent nonlinear eigensolve, metric root, Pi_W and X.
        # Hold the physical central Hamiltonian fixed in the scalar so this
        # isolates its projector derivative, excluding direct Tr(rho dh).
        _, _, u = endpoint(self.full, h, o, p, t)
        return self.full.lattice_scalar(u=u)

    def test_joint_gradient_against_nonlinear_full_scalar(self):
        h, o, p = self.perturbation()
        expected = self.gradient @ self.packed(h, o, p)
        errors = []
        for dt in [4e-4, 2e-4]:
            fd = (self.scalar(h, o, p, dt) - self.scalar(h, o, p, -dt)) / (2 * dt)
            errors.append(abs(fd - expected))
        self.assertLess(errors[-1], 2e-8)
        self.assertGreater(errors[0] / errors[1], 3.0)

    def test_separate_overlap_raw_and_complement_legs(self):
        h, o, p = self.perturbation()
        n = self.nr
        cases = [(h, 0 * o, 0 * p), (0 * h, o, 0 * p), (0 * h, 0 * o, p)]
        for index, (hh, oo, pp) in enumerate(cases):
            # Isolate the cross blocks in the first two cases. This catches
            # both the QR pairing halves and the retained-energy sign in O.
            hh = hh.copy()
            oo = oo.copy()
            for x in [hh, oo]:
                x[:, :n, :n] = 0
                x[:, n:, n:] = 0
            dt = 2e-5
            fd = (self.scalar(hh, oo, pp, dt) - self.scalar(hh, oo, pp, -dt)) / (2 * dt)
            expected = self.gradient @ self.packed(hh, oo, pp)
            np.testing.assert_allclose(fd, expected, rtol=3e-6, atol=2e-9)

    def test_retaining_all_states_gives_same_gradient(self):
        h, o, p = self.perturbation()
        hb, pb = self.full.projector_gradients()
        ob = -0.5 * (hb * self.full.eps[:, None, :] + self.full.eps[:, :, None] * hb)
        ob -= 0.5 * herm(np.einsum("knic,kmic->knm", pb, self.full.raw.conj()))
        expected = real_inner(hb, h) + real_inner(ob, o) + real_inner(pb, p)
        np.testing.assert_allclose(
            self.gradient @ self.packed(h, o, p), expected, rtol=3e-12, atol=3e-13
        )

    def test_zero_sigma_gives_zero_functional_gradient(self):
        old = self.small
        meta = replace(
            old.meta,
            sigma_inf=np.zeros_like(old.meta.sigma_inf),
            sigma=np.zeros_like(old.meta.sigma),
        )
        small = CachedFrozenSigmaResponse(
            meta, old.eps, old.raw, old.apply_x(old.raw, old.X), old.kw, old.W
        )
        m = EmptyComplementResponse(small, self.m.solvers, self.m.raw_q)
        g, _ = projector_functional_gradient(m)
        np.testing.assert_array_equal(g, np.zeros_like(g))

    def test_gauge_normalization_and_basis_covariance(self):
        f = self.full
        n = f.eps.shape[1]
        nk = len(f.eps)
        zero = np.zeros((nk, n, n), complex)
        identity = np.broadcast_to(np.eye(n), zero.shape)
        self.assertLess(
            abs(self.gradient @ self.packed(identity, zero, 0 * f.raw)), 1e-14
        )
        # Global raw-amplitude change must be cancelled by derivative of X.
        self.assertLess(abs(self.gradient @ self.packed(zero, zero, f.raw)), 2e-13)
        phases = np.exp(1j * self.rng.normal(size=(nk, n)))
        # Infinitesimal diagonal phases commute with central H.
        a = np.array([np.diag(1j * np.angle(p)) for p in phases])
        ad = a.swapaxes(1, 2).conj()
        raw = np.einsum("knm,kmic->knic", ad, f.raw)
        self.assertLess(abs(self.gradient @ self.packed(zero, zero, raw)), 2e-13)
        # A general basis change need not cancel in the projector part alone.
        # The independently computed direct density/H contribution cancels it.
        a = 0.07 * (
            self.rng.normal(size=zero.shape) + 1j * self.rng.normal(size=zero.shape)
        )
        ad = a.swapaxes(1, 2).conj()
        h = ad * f.eps[:, None, :] + f.eps[:, :, None] * a
        o = ad + a
        raw = np.einsum("knm,kmic->knic", ad, f.raw)
        heff = h - 0.5 * (o * f.eps[:, None, :] + f.eps[:, :, None] * o)
        total = self.gradient @ self.packed(h, o, raw) + real_inner(
            f.central_density(), heff
        )
        self.assertLess(abs(total), 2e-12)

    def test_native_pullback_and_one_adjoint_multiple_sources(self):
        rng = self.rng
        nq = 4
        nv = 5
        b = rng.normal(size=(len(self.gradient), nv)) * 0.03
        k = rng.normal(size=(nv, nq)) * 0.1
        j = rng.normal(size=(nq, nq)) * 0.1
        state = NS(
            valence=NS(complement=self.m),
            chain=NS(basis=dense(b), potential=dense(k)),
            jacobian=dense(j),
            screened=ScreenedResponse(dense(j)),
        )
        obs = PreparedProjectorFunctional(state)
        s = rng.normal(size=(nq, 3))
        v = rng.normal(size=(nv, 3))
        actual, lam, report = obs.adjoint_contractions(s, v, rtol=1e-12, atol=1e-14)
        q = np.linalg.solve(np.eye(nq) - j, s)
        expected = self.gradient @ b @ (v + k @ q)
        np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-14)
        self.assertFalse(report["stress_applied"])
        self.assertFalse(report["includes_direct_density_hamiltonian_term"])
        for nu in range(3):
            self.assertAlmostEqual(
                obs.potential_contraction(v[:, nu] + k @ q[:, nu]),
                actual[nu],
                places=12,
            )
        with self.assertRaises(ValueError):
            obs.adjoint_contractions(s, np.zeros((nv, 2)))
        with self.assertRaises(ValueError):
            obs.adjoint_contractions(s, None)
        with self.assertRaises(ValueError):
            obs.potential_contraction(np.zeros(nv + 1))


if __name__ == "__main__":
    unittest.main()
