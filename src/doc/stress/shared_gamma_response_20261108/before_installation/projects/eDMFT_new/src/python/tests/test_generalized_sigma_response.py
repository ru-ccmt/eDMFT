"""Independent nonlinear endpoints and moving-basis covariance controls."""

import unittest
from unittest.mock import patch

import numpy as np

from frozen_sigma_response import herm, invsqrt
from frozen_sigma_tangent import real_inner
from test_frozen_sigma_tangent import complex_model


def endpoint(model, dh, do, draw, step):
    """Nonlinear Loewdin endpoint, independently rediagonalized and normalized."""
    nb = model.eps.shape[1]
    transform = np.array([invsqrt(np.eye(nb) + step * o) for o in do])
    h = transform @ (np.array([np.diag(e) for e in model.eps]) + step * dh) @ transform
    raw = np.einsum("knm,kmic->knic", transform, model.raw + step * draw)
    q = np.zeros_like(raw)
    for k in range(len(h)):
        _, v = np.linalg.eigh(h[k])
        window = v[:, model.W]
        q[k] = np.einsum("nm,mic->nic", window @ window.conj().T, raw[k])
    overlap = model.overlap(q) + model.offset
    x = np.zeros_like(model.X)
    for c, n in enumerate(model.meta.dimensions):
        x[:n, :n, c] = np.eye(n)
        active = np.flatnonzero(model.meta.active[:n, c])
        ix = np.ix_(active, active)
        x[:, :, c][ix] = invsqrt(overlap[:, :, c][ix])
    u = model.meta.compact(model.apply_x(q, x))
    delta = h - np.array([np.diag(e) for e in model.eps])
    density = model.density(delta, 1.0, u)
    return transform @ density @ transform, delta, u


class GeneralizedResponseTests(unittest.TestCase):
    def setUp(self):
        self.model = complex_model()
        rng = np.random.default_rng(386)
        self.h = herm(rng.normal(size=(2, 5, 5)) + 1j * rng.normal(size=(2, 5, 5)))
        self.o = herm(rng.normal(size=(2, 5, 5)) + 1j * rng.normal(size=(2, 5, 5)))
        self.y = herm(rng.normal(size=(2, 5, 5)) + 1j * rng.normal(size=(2, 5, 5)))
        self.raw = 0.2 * (
            rng.normal(size=self.model.raw.shape)
            + 1j * rng.normal(size=self.model.raw.shape)
        )

    def test_central_density_matches_full_resolvent_oracle(self):
        np.testing.assert_allclose(
            self.model.central_density(), self.model.density(), atol=2e-15
        )

    def test_explicit_raw_endpoint_and_pullback(self):
        m = self.model
        forward, _ = m.response(self.h, fixed_number=False, raw_derivative=self.raw)
        zero = np.zeros_like(self.o)
        h = 1e-5
        plus, _, _ = endpoint(m, self.h, zero, self.raw, h)
        minus, _, _ = endpoint(m, self.h, zero, self.raw, -h)
        np.testing.assert_allclose(
            (plus - minus) / (2 * h), forward, atol=2e-7, rtol=2e-7
        )
        for fixed_number in (False, True):
            f, _ = m.response(
                self.h, fixed_number=fixed_number, raw_derivative=self.raw
            )
            hb, rb = m.adjoint(self.y, fixed_number=fixed_number, return_raw=True)
            np.testing.assert_allclose(
                real_inner(self.y, f),
                real_inner(hb, self.h) + real_inner(rb, self.raw),
                atol=1e-11,
            )

    def test_outside_window_explicit_raw_is_projected_out(self):
        raw = self.raw.copy()
        raw[:, self.model.W] = 0
        f, _, _ = self.model.embedding_response(np.zeros_like(self.h), raw)
        np.testing.assert_array_equal(f, 0)

    def test_generalized_nonlinear_endpoint_quadratic_convergence(self):
        m = self.model
        f, _ = m.generalized_response(self.h, self.o, self.raw, fixed_number=False)
        errors = []
        for h in (1e-4, 5e-5):
            plus, _, _ = endpoint(m, self.h, self.o, self.raw, h)
            minus, _, _ = endpoint(m, self.h, self.o, self.raw, -h)
            errors.append(np.linalg.norm((plus - minus) / (2 * h) - f))
        self.assertLess(errors[-1], 1e-5)
        self.assertGreater(errors[0] / errors[1], 3.8)

    def test_generalized_adjoint_both_ensembles_and_cache_modes(self):
        for cache_mb in (0, 0.0035, 4):
            m = complex_model(cache_mb=cache_mb)
            for fixed_number in (False, True):
                f, info = m.generalized_response(self.h, self.o, self.raw, fixed_number)
                hb, ob, rb = m.generalized_adjoint(self.y, fixed_number)
                np.testing.assert_allclose(
                    real_inner(self.y, f),
                    real_inner(hb, self.h)
                    + real_inner(ob, self.o)
                    + real_inner(rb, self.raw),
                    atol=2e-11,
                    rtol=2e-12,
                )
                if fixed_number:
                    self.assertLess(abs(info["dN"]), 1e-13)
                    self.assertGreater(abs(info["coefficient_trace"]), 1e-3)

    def test_pure_basis_change_cancels_physical_density(self):
        m = self.model
        # Any complex A gives B(t)=B(0)(I+tA). This changes coefficients,
        # overlap, Hamiltonian, and raw projector, but no physical operator.
        a = 0.3 * self.h + 0.2j * self.o
        ad = a.swapaxes(-1, -2).conj()
        dh = ad * m.eps[:, None, :] + m.eps[:, :, None] * a
        do = ad + a
        draw = np.einsum("knm,kmic->knic", ad, m.raw)
        rho = m.central_density()
        for fixed_number in (False, True):
            f, info = m.generalized_response(dh, do, draw, fixed_number)
            np.testing.assert_allclose(f + a @ rho + rho @ ad, 0, atol=2e-13)
            self.assertLess(abs(info["dN"]), 2e-13)
            self.assertLess(abs(info["dmu"]), 2e-13)

    def test_generalized_warm_action_has_no_eigensolve(self):
        m = self.model
        m.prepare()
        with patch(
            "frozen_sigma_tangent.eigh",
            side_effect=AssertionError("eigh in warm action"),
        ), patch(
            "numpy.linalg.inv", side_effect=AssertionError("inverse in warm action")
        ):
            m.generalized_response(self.h, self.o, self.raw)
            m.generalized_adjoint(self.y)

    def test_common_projector_gradient_includes_explicit_raw(self):
        m = self.model
        hb, rb = m.projector_gradients()
        h = 1e-5
        _, _, up = endpoint(m, self.h, np.zeros_like(self.o), self.raw, h)
        _, _, um = endpoint(m, self.h, np.zeros_like(self.o), self.raw, -h)
        difference = (m.lattice_scalar(u=up) - m.lattice_scalar(u=um)) / (2 * h)
        np.testing.assert_allclose(
            difference,
            real_inner(hb, self.h) + real_inner(rb, self.raw),
            atol=1e-8,
            rtol=1e-7,
        )

    def test_bad_raw_derivative_rejected(self):
        with self.assertRaisesRegex(ValueError, "raw-projector"):
            self.model.response(self.h, raw_derivative=self.raw[..., 0])
        raw = self.raw.copy()
        raw[0, 0, 0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "raw-projector"):
            self.model.generalized_response(self.h, self.o, raw)


if __name__ == "__main__":
    unittest.main()
