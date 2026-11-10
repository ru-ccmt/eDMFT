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

    def inverse_endpoint(self, step, k, frequencies, dmu):
        """Independent nonlinear overlap/window/normalization reconstruction."""
        m = self.model
        _, _, u = endpoint(m, self.h, self.o, self.raw, step)
        metric = np.eye(m.eps.shape[1]) + step * self.o[k]
        e, v = np.linalg.eigh(metric)
        root = (v * np.sqrt(e)) @ v.conj().T
        covariant_u = root @ u[k]
        z = m.meta.mu + step * dmu + 1j * m.meta.omega[frequencies]
        return (
            z[:, None, None] * metric
            - (np.diag(m.eps[k]) + step * self.h[k])
            - (covariant_u @ m.sigma[frequencies]) @ covariant_u.conj().T
        )

    def test_inverse_green_endpoint_including_metric_and_normalization(self):
        m = self.model
        frequencies = np.array([0, 4, 17, 44])
        dmu = 0.21
        prepared = m.prepare_inverse_response(self.h, self.o, self.raw, dmu=dmu)
        for k in range(len(m.eps)):
            actual = prepared.block(k, frequencies)
            # Sigma(iw) is not Hermitian. Projection onto Hermitian output
            # would erase part of the overlap and dynamical embedding terms.
            self.assertGreater(np.linalg.norm(actual - herm(actual)), 1.0)
            errors = []
            for h in [1e-4, 5e-5]:
                plus = self.inverse_endpoint(h, k, frequencies, dmu)
                minus = self.inverse_endpoint(-h, k, frequencies, dmu)
                errors.append(np.linalg.norm((plus - minus) / (2 * h) - actual))
            self.assertLess(errors[-1], 1e-5)
            self.assertGreater(errors[0] / errors[1], 3.7)
            center = self.inverse_endpoint(0.0, k, frequencies, dmu)
            g = np.linalg.inv(center)
            analytic_g = -g @ actual @ g
            fd_g = (np.linalg.inv(plus) - np.linalg.inv(minus)) / (2 * h)
            np.testing.assert_allclose(analytic_g, fd_g, rtol=2e-6, atol=2e-6)

    def test_inverse_green_pure_basis_change_has_correct_covariance(self):
        m = self.model
        a = 0.3 * self.h + 0.2j * self.o
        ad = a.swapaxes(-1, -2).conj()
        dh = ad * m.eps[:, None, :] + m.eps[:, :, None] * a
        do = ad + a
        raw = np.einsum("knm,kmic->knic", ad, m.raw)
        prepared = m.prepare_inverse_response(dh, do, raw)
        frequencies = np.array([0, 9, 44])
        for k, u in enumerate(m.uc):
            z = m.meta.mu + 1j * m.meta.omega[frequencies]
            central = (
                z[:, None, None] * np.eye(len(m.eps[k]))
                - np.diag(m.eps[k])
                - (u @ m.sigma[frequencies]) @ u.conj().T
            )
            expected = ad[k] @ central + central @ a[k]
            np.testing.assert_allclose(
                prepared.block(k, frequencies), expected, rtol=2e-13, atol=2e-13
            )
            # Pulling A1 back to the original basis removes a pure coordinate
            # change: no physical electron-phonon perturbation remains.
            np.testing.assert_allclose(
                prepared.block(k, frequencies) - expected, 0, atol=2e-13
            )

    def test_inverse_green_fixed_number_uniform_shift(self):
        m = self.model
        identity = np.broadcast_to(np.eye(5), self.h.shape)
        _, info = m.generalized_response(
            identity, np.zeros_like(self.o), fixed_number=True
        )
        readout = m.prepare_inverse_response(
            identity, np.zeros_like(self.o), dmu=info["dmu"]
        )
        for k in range(len(m.eps)):
            np.testing.assert_allclose(readout.block(k, slice(None)), 0, atol=3e-14)

    def test_inverse_green_readout_is_prepared_once_and_does_not_invert(self):
        m = self.model
        h, o, raw = self.h.copy(), self.o.copy(), self.raw.copy()
        with patch(
            "frozen_sigma_tangent.eigh", side_effect=AssertionError("eigh in readout")
        ), patch("numpy.linalg.inv", side_effect=AssertionError("inverse in readout")):
            readout = m.prepare_inverse_response(h, o, raw)
            expected = readout.block(0, slice(None))
            # Caller-owned perturbations must not change a prepared readout.
            h[:] = 0
            o[:] = 0
            raw[:] = 0
            with patch.object(
                m,
                "embedding_response",
                side_effect=AssertionError("repeated projector preparation"),
            ):
                actual = np.concatenate(
                    [readout.block(0, slice(0, 20)), readout.block(0, slice(20, None))]
                )
        np.testing.assert_array_equal(actual, expected)
        self.assertFalse(m._ready)  # No costly density-response cache needed.

    def test_inverse_green_invalid_frequency_and_mu_are_rejected(self):
        m = self.model
        with self.assertRaises(ValueError):
            m.prepare_inverse_response(self.h, self.o, dmu=1j)
        prepared = m.prepare_inverse_response(self.h, self.o)
        for index in [[], [-1], [45], [1.5], [[1, 2]]]:
            with self.subTest(index=index), self.assertRaises(ValueError):
                prepared.block(0, index)
        with self.assertRaises(ValueError):
            prepared.block(2, [0])
        m.meta.mu += 0.1
        with self.assertRaisesRegex(ValueError, "central mu"):
            m.prepare_inverse_response(self.h, self.o)

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
