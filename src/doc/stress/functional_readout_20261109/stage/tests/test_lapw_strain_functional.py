"""Independent nonlinear lattice endpoints and screened strain contractions."""

from dataclasses import replace
from types import SimpleNamespace as NS
import unittest

import numpy as np

from frozen_sigma_response import herm
from frozen_sigma_tangent import CachedFrozenSigmaResponse, real_inner
from lapw_complement_response import EmptyComplementResponse
from lapw_strain_functional import (
    lattice_functional_gradients,
    PreparedStrainFunctional,
)
from stress_response_adjoint import ScreenedResponse
from test_generalized_sigma_response import endpoint
from test_lapw_complement_response import fixture
from test_lapw_projector_functional import dense


class StrainFunctionalTests(unittest.TestCase):
    def setUp(self):
        self.small, self.full, self.m, self.t, self.rng = fixture()
        self.gradients, self.coordinates = lattice_functional_gradients(self.m)
        self.nr = self.small.eps.shape[1]

    def packed(self, h, o, p):
        td = self.t.swapaxes(1, 2).conj()
        h, o = td @ h @ self.t, td @ o @ self.t
        p = np.einsum("knm,kmic->knic", td, p)
        n = self.nr
        return self.coordinates.pack(
            h[:, :n, :n], o[:, :n, :n], p[:, :n], h[:, n:, :n], o[:, n:, :n]
        )

    def random_direction(self):
        shape = (len(self.full.eps), self.full.eps.shape[1], self.full.eps.shape[1])
        random = lambda s: self.rng.normal(size=s) + 1j * self.rng.normal(size=s)
        return (
            0.15 * herm(random(shape)),
            0.07 * herm(random(shape)),
            0.07 * random(self.full.raw.shape),
        )

    def scalar(self, h, o, p, step):
        _, dh, u = endpoint(self.full, h, o, p, step)
        # Unlike the projector-only oracle, BOTH the physical orthonormal H
        # and normalized U change. No analytic covectors enter this endpoint.
        return self.full.lattice_scalar(dh, 1.0, u)

    def test_full_lattice_gradient_matches_nonlinear_endpoint(self):
        h, o, p = self.random_direction()
        analytic = self.gradients["lattice"] @ self.packed(h, o, p)
        errors = []
        for step in [4e-4, 2e-4]:
            fd = (self.scalar(h, o, p, step) - self.scalar(h, o, p, -step)) / (2 * step)
            errors.append(abs(fd - analytic))
        self.assertLess(errors[-1], 3e-8)
        self.assertGreater(errors[0] / errors[1], 3)

    def test_direct_piece_matches_independent_full_density(self):
        h, o, p = self.random_direction()
        e = self.full.eps
        heff = h - 0.5 * (o * e[:, None, :] + e[:, :, None] * o)
        expected = real_inner(self.full.central_density(), heff)
        np.testing.assert_allclose(
            self.gradients["direct"] @ self.packed(h, o, p),
            expected,
            rtol=2e-12,
            atol=1e-13,
        )

    def test_pure_basis_change_cancels_direct_and_projector(self):
        h, o, _ = self.random_direction()
        a = h + 1j * o
        ad = a.swapaxes(1, 2).conj()
        e = self.full.eps
        dh = ad * e[:, None, :] + e[:, :, None] * a
        do = ad + a
        dp = np.einsum("knm,kmic->knic", ad, self.full.raw)
        b = self.packed(dh, do, dp)
        self.assertGreater(abs(self.gradients["projector"] @ b), 1e-5)
        self.assertLess(abs(self.gradients["lattice"] @ b), 3e-12)

    def test_zero_sigma_is_dft_direct_term(self):
        old = self.small
        meta = replace(
            old.meta,
            sigma_inf=np.zeros_like(old.meta.sigma_inf),
            sigma=np.zeros_like(old.meta.sigma),
        )
        small = CachedFrozenSigmaResponse(
            meta, old.eps, old.raw, old.apply_x(old.raw, old.X), old.kw, old.W
        )
        complement = EmptyComplementResponse(small, self.m.solvers, self.m.raw_q)
        g, _ = lattice_functional_gradients(complement)
        np.testing.assert_array_equal(g["projector"], 0)
        np.testing.assert_array_equal(g["lattice"], g["direct"])
        rho = small.central_density()
        f = (
            2
            * small.kw[:, None]
            / (1 + np.exp((small.eps - small.meta.mu) / small.meta.temperature))
        )
        np.testing.assert_allclose(
            rho, np.array([np.diag(row) for row in f]), atol=1e-14
        )

    def test_uniform_energy_shift_derivative_is_particle_number(self):
        n = self.full.eps.shape[1]
        identity = np.broadcast_to(np.eye(n), (len(self.full.eps), n, n))
        b = self.packed(identity, 0 * identity, 0 * self.full.raw)
        number = np.trace(self.full.central_density(), axis1=1, axis2=2).real.sum()
        self.assertAlmostEqual(self.gradients["lattice"] @ b, number, places=12)
        self.assertLess(abs(self.gradients["projector"] @ b), 1e-14)

    def synthetic_state(self):
        b = 0.03 * self.rng.normal(size=(self.coordinates.size, 5))
        k = 0.1 * self.rng.normal(size=(5, 4))
        j = 0.1 * self.rng.normal(size=(4, 4))
        state = NS(
            valence=NS(complement=self.m),
            chain=NS(basis=dense(b), potential=dense(k)),
            jacobian=dense(j),
            screened=ScreenedResponse(dense(j)),
        )
        return state, b, k, j

    def test_six_explicit_sources_forward_and_single_adjoint(self):
        state, b, k, j = self.synthetic_state()
        geo = 0.01 * self.rng.normal(size=(self.coordinates.size, 6))
        pot = self.rng.normal(size=(5, 6))
        sources = self.rng.normal(size=(4, 6))
        q = np.linalg.solve(np.eye(4) - j, sources)
        values = {}
        for component in ["projector", "direct", "lattice"]:
            obs = PreparedStrainFunctional(state, component=component)
            forward, report = obs.forward_contractions(geo, pot, q)
            adjoint, lam, report = obs.adjoint_contractions(
                geo, pot, sources, rtol=1e-12, atol=1e-14
            )
            expected = self.gradients[component] @ (geo + b @ (pot + k @ q))
            np.testing.assert_allclose(forward, expected, rtol=2e-12, atol=2e-14)
            np.testing.assert_allclose(adjoint, expected, rtol=2e-12, atol=2e-14)
            self.assertFalse(report["full_klein_functional"])
            self.assertFalse(report["baseline_subtracted"])
            self.assertFalse(report["stress_applied"])
            values[component] = forward
        np.testing.assert_allclose(
            values["lattice"], values["direct"] + values["projector"], atol=1e-14
        )

    def test_missing_geometry_or_wrong_layout_is_rejected(self):
        state, _, _, _ = self.synthetic_state()
        obs = PreparedStrainFunctional(state)
        geo = np.zeros((self.coordinates.size, 1))
        for g, v, q in [
            (None, np.zeros((5, 1)), np.zeros((4, 1))),
            (geo, None, np.zeros((4, 1))),
            (geo, np.zeros((5, 2)), np.zeros((4, 1))),
            (geo, np.zeros((5, 1)), np.zeros((4, 2))),
        ]:
            with self.assertRaises(ValueError):
                obs.forward_contractions(g, v, q)
        with self.assertRaises(ValueError):
            PreparedStrainFunctional(state, component="total_stress")


if __name__ == "__main__":
    unittest.main()
