"""Independent finite differences, reverse contractions, and cache checks."""

import unittest
from dataclasses import replace
from unittest.mock import patch

import numpy as np

from frozen_sigma_response import Metadata, FrozenSigmaResponse, invsqrt, herm
from frozen_sigma_tangent import CachedFrozenSigmaResponse, real_inner
from test_frozen_sigma_response import model


def complex_model(cls=CachedFrozenSigmaResponse, cache_mb=4):
    """Two inequivalent k points, complex orbitals and nontrivial real symmetry.

    Deliberately use a non-self-adjoint Hermiticity-preserving symmetry map.
    An implementation that silently assumes every exported map is an
    orthogonal projector will fail this stronger algebraic control.
    """
    rng = np.random.default_rng(91)
    t = 0.07
    om = np.pi * t * (2 * np.arange(45) + 1)
    weights = rng.uniform(0.6, 1.8, len(om))
    raw = 0.2 * (rng.normal(size=(2, 5, 2, 1)) + 1j * rng.normal(size=(2, 5, 2, 1)))
    rotation = np.array([[1, 0.2j], [0.1, 0.8]], complex)

    def sym(a):
        return herm(rotation @ herm(a) @ rotation.conj().T)

    re, im = [], []
    for j in range(4):
        unit = np.zeros(4, complex)
        unit[j] = 1
        re.append(sym(unit.reshape((2, 2), order="F")).ravel(order="F"))
        im.append(sym(1j * unit.reshape((2, 2), order="F")).ravel(order="F"))
    kw = np.array([0.3, 0.7])
    x = invsqrt(
        sym(np.einsum("k,kni,knj->ij", kw, raw[:, 1:4, :, 0].conj(), raw[:, 1:4, :, 0]))
    )
    si = np.array([[0.18, 0.03j], [-0.03j, -0.1]])
    pole = np.array([[0.13, 0.02 + 0.015j], [0.02 - 0.015j, 0.08]])
    sigma = si[:, :, None] + pole[:, :, None] / (1j * om - 0.5)
    meta = Metadata(
        0.1,
        t,
        om,
        weights,
        np.array([2]),
        np.array([0]),
        np.array([2]),
        np.array([[0], [1]]),
        np.ones((2, 1), bool),
        x[:, :, None],
        si[:, :, None],
        sigma[:, :, None, :],
        np.array(re).T,
        np.array(im).T,
    )
    energies = np.array([[-1.4, -0.3, 0.08, 0.5, 1.5], [-1.1, -0.2, 0.1, 0.6, 1.8]])
    normalized = np.einsum("knic,ij->knjc", raw, x)
    args = (meta, energies, raw, normalized, kw, [1, 2, 3])
    return cls(
        *args, **({"cache_mb": cache_mb} if cls is CachedFrozenSigmaResponse else {})
    )


class TangentTests(unittest.TestCase):
    def setUp(self):
        self.m = complex_model()
        rng = np.random.default_rng(73)
        shape = (2, 5, 5)
        self.x = herm(rng.normal(size=shape) + 1j * rng.normal(size=shape))
        self.y = herm(rng.normal(size=shape) + 1j * rng.normal(size=shape))

    def test_cached_forward_matches_independent_reference(self):
        reference = complex_model(FrozenSigmaResponse)
        for fixed_n in [False, True]:
            for embedding in [False, True]:
                expected, a = reference.response(self.x, fixed_n, embedding)
                actual, b = self.m.response(self.x, fixed_n, embedding)
                np.testing.assert_allclose(actual, expected, atol=3e-14, rtol=3e-13)
                self.assertAlmostEqual(a["dmu"], b["dmu"], places=12)

    def test_complex_normalized_endpoint(self):
        actual, _ = self.m.response(self.x, fixed_number=False)
        errors = []
        for h in [1e-4, 5e-5]:
            plus = self.m.density(self.x, h, self.m.embedding_endpoint(self.x, h))
            minus = self.m.density(self.x, -h, self.m.embedding_endpoint(self.x, -h))
            errors.append(np.linalg.norm((plus - minus) / (2 * h) - actual))
        self.assertLess(errors[-1], 3e-6)
        self.assertGreater(errors[0] / errors[1], 3.8)

    def test_real_linear_symmetry_transpose(self):
        rng = np.random.default_rng(20)
        x = rng.normal(size=self.m.X.shape) + 1j * rng.normal(size=self.m.X.shape)
        y = rng.normal(size=x.shape) + 1j * rng.normal(size=x.shape)
        self.assertAlmostEqual(
            real_inner(y, self.m.meta.sym(x)),
            real_inner(self.m.symmetry_adjoint(y), x),
            places=13,
        )

    def test_embedding_transpose(self):
        rng = np.random.default_rng(28)
        y = rng.normal(size=self.m.uc.shape) + 1j * rng.normal(size=self.m.uc.shape)
        forward, _, _ = self.m.embedding_response(self.x)
        self.assertAlmostEqual(
            real_inner(y, forward),
            real_inner(self.m.embedding_adjoint(y), self.x),
            places=12,
        )

    def test_density_transpose_all_constraints(self):
        for fixed_n in [False, True]:
            for embedding in [False, True]:
                f, _ = self.m.response(self.x, fixed_n, embedding)
                b = self.m.adjoint(self.y, fixed_n, embedding)
                np.testing.assert_allclose(
                    real_inner(self.y, f), real_inner(b, self.x), atol=2e-12, rtol=2e-12
                )

    def test_fixed_n_charge_covector_is_null(self):
        identity = np.broadcast_to(np.eye(5), self.x.shape)
        np.testing.assert_allclose(self.m.adjoint(identity), 0, atol=1e-14)
        rho, report = self.m.response(identity)
        np.testing.assert_allclose(rho, 0, atol=2e-14)
        self.assertAlmostEqual(report["dmu"], 1.0, places=12)

    def test_warm_actions_do_not_diagonalize_or_invert(self):
        self.m.prepare()
        with patch(
            "frozen_sigma_tangent.eigh",
            side_effect=AssertionError("unexpected diagonalization"),
        ), patch(
            "numpy.linalg.inv", side_effect=AssertionError("unexpected inversion")
        ):
            self.m.response(self.x)
            self.m.adjoint(self.y)
        self.assertLessEqual(self.m.cache_info()["resolvent_bytes"], self.m.cache_limit)

    def test_no_cache_agrees_in_both_directions(self):
        uncached = complex_model(cache_mb=0)
        for action in ["response", "adjoint"]:
            a = getattr(self.m, action)(self.x)
            b = getattr(uncached, action)(self.x)
            np.testing.assert_allclose(
                a[0] if action == "response" else a,
                b[0] if action == "response" else b,
                atol=1e-14,
            )
        self.assertEqual(uncached.cache_info()["resolvent_bytes"], 0)

    def test_empty_complement(self):
        original = model()
        args = (
            original.meta,
            original.eps,
            original.raw,
            original.apply_x(original.raw, original.X),
            original.kw,
            np.arange(4),
        )
        m = CachedFrozenSigmaResponse(*args)
        x, y = self.x[:, :4, :4], self.y[:, :4, :4]
        du, _, _ = m.embedding_response(x)
        np.testing.assert_array_equal(du, 0)
        f, _ = m.response(x)
        self.assertAlmostEqual(real_inner(y, f), real_inner(m.adjoint(y), x), places=12)

    def test_scalar_projector_gradient_without_direct_density_term(self):
        expected = real_inner(self.m.projector_potential_gradient(), self.x)
        h = 1e-5
        plus = self.m.lattice_scalar(u=self.m.embedding_endpoint(self.x, h))
        minus = self.m.lattice_scalar(u=self.m.embedding_endpoint(self.x, -h))
        self.assertAlmostEqual((plus - minus) / (2 * h), expected, places=7)

    def test_complete_scalar_gradient(self):
        gradient = self.m.density() + self.m.projector_potential_gradient()
        h = 1e-5
        plus = self.m.lattice_scalar(self.x, h, self.m.embedding_endpoint(self.x, h))
        minus = self.m.lattice_scalar(self.x, -h, self.m.embedding_endpoint(self.x, -h))
        self.assertAlmostEqual(
            (plus - minus) / (2 * h), real_inner(gradient, self.x), places=7
        )

    def test_central_arrays_are_read_only(self):
        for array in [self.m.eps, self.m.sigma, self.m.X, self.m.meta.weights]:
            with self.assertRaises(ValueError):
                array.flat[0] = 0

    def test_changed_mu_rejects_stale_context(self):
        self.m.prepare()
        self.m.meta.mu += 0.01
        with self.assertRaisesRegex(ValueError, "construct a new response context"):
            self.m.response(self.x)

    def test_partial_cache_matches_full_cache(self):
        partial = complex_model(cache_mb=0.0035)
        a, _ = partial.response(self.x)
        b, _ = self.m.response(self.x)
        np.testing.assert_allclose(a, b, atol=3e-14)
        np.testing.assert_allclose(
            partial.adjoint(self.y), self.m.adjoint(self.y), atol=3e-14
        )
        self.assertEqual(partial.cache_info()["precontracted_kpoints"], 1)
        self.assertLessEqual(
            partial.cache_info()["total_cache_bytes"], partial.cache_limit
        )

    def test_zero_self_energy_removes_all_projector_response(self):
        m = self.m
        meta = replace(
            m.meta,
            sigma=np.zeros_like(m.meta.sigma),
            sigma_inf=np.zeros_like(m.meta.sigma_inf),
        )
        zero = CachedFrozenSigmaResponse(
            meta, m.eps, m.raw, m.apply_x(m.raw, m.X), m.kw, m.W
        )
        np.testing.assert_allclose(zero.projector_potential_gradient(), 0, atol=1e-15)
        a, _ = zero.response(self.x, embedding=True)
        b, _ = zero.response(self.x, embedding=False)
        np.testing.assert_allclose(a, b, atol=1e-15)
        # With Sigma=0, equal k weights are not required: per-k scalar
        # weights commute with the ordinary independent-electron response.
        np.testing.assert_allclose(zero.adjoint(self.x), a, atol=3e-14)


if __name__ == "__main__":
    unittest.main()
