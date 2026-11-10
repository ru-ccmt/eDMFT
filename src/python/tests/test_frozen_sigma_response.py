"""Small independent derivative and screening checks (no WIEN installation)."""

import unittest
import numpy as np
from frozen_sigma_response import Metadata, FrozenSigmaResponse
from frozen_sigma_screening import solve_response


def model():
    temperature = 0.05
    omega = (2 * np.arange(160) + 1) * np.pi * temperature
    raw = np.array(
        [[[[0.3]], [[0.7]], [[0.5j]], [[0.2]]], [[[0.4j]], [[0.6]], [[0.4]], [[0.3j]]]],
        complex,
    )
    x = np.array([[[1 / np.sqrt(np.sum(abs(raw[:, 1:3]) ** 2) / 2)]]])
    sigma_inf = np.array([[[0.12]]], complex)
    sigma = (0.12 + 0.2 / (1j * omega - 0.6)).reshape(1, 1, 1, -1)
    meta = Metadata(
        0.0,
        temperature,
        omega,
        np.ones_like(omega),
        np.array([1]),
        np.array([0]),
        np.array([1]),
        np.array([[0]]),
        np.array([[True]]),
        x,
        sigma_inf,
        sigma,
        np.ones((1, 1), complex),
        np.zeros((1, 1), complex),
    )
    energies = np.array([[-1.2, -0.2, 0.3, 1.3], [-1.1, -0.15, 0.4, 1.4]])
    return FrozenSigmaResponse(meta, energies, raw, raw * x, [0.5, 0.5], [1, 2])


class ResponseTests(unittest.TestCase):
    def test_resolvent_and_normalization_derivative(self):
        m = model()
        rng = np.random.default_rng(72)
        a = rng.normal(size=(2, 4, 4)) + 1j * rng.normal(size=(2, 4, 4))
        perturbation = 0.1 * (a + a.swapaxes(-1, -2).conj())
        analytic, _ = m.response(perturbation, fixed_number=False)
        errors = []
        for h in [2e-4, 1e-4]:
            plus = m.density(perturbation, h, m.embedding_endpoint(perturbation, h))
            minus = m.density(perturbation, -h, m.embedding_endpoint(perturbation, -h))
            errors.append(np.linalg.norm((plus - minus) / (2 * h) - analytic))
        self.assertLess(errors[1], 2e-7)
        self.assertGreater(errors[0] / errors[1], 3.8)

    def test_uniform_shift_and_fixed_particle_number(self):
        m = model()
        perturbation = np.broadcast_to(np.eye(4), (2, 4, 4))
        du, ds, dx = m.embedding_response(perturbation)
        self.assertLess(np.max(abs(du)), 1e-14)
        rho, report = m.response(perturbation)
        self.assertLess(np.max(abs(rho)), 1e-13)
        self.assertAlmostEqual(report["dmu"], 1.0, places=12)

    def test_reject_nonhermitian_perturbation(self):
        m = model()
        perturbation = np.zeros((2, 4, 4), complex)
        perturbation[0, 0, 1] = 1
        with self.assertRaisesRegex(ValueError, "Hermitian"):
            m.response(perturbation)

    def test_screening_unstable_under_plain_mixing(self):
        # One response eigenvalue makes beta=.3 simple mixing divergent.
        a = np.array([[-7.0, 0.4, 0.0], [0.1, 0.3, 0.2], [0.0, -0.2, -0.5]])
        source = np.array([1.0, 2.0, -1.0])
        answer, report = solve_response(
            source, lambda x: a @ x, absolute_tolerance=1e-10, max_iterations=12
        )
        self.assertTrue(report["converged"])
        np.testing.assert_allclose(
            answer, np.linalg.solve(np.eye(3) - a, source), atol=1e-11
        )

    def test_nonconvergence_returns_verified_iterate(self):
        source = np.array([1.0])
        answer, report = solve_response(source, lambda x: 2 * x, max_iterations=1)
        self.assertFalse(report["converged"])
        self.assertEqual(answer[0], 1.0)
        self.assertEqual(report["iterations"][-1]["residual_norm"], 2.0)


if __name__ == "__main__":
    unittest.main()
