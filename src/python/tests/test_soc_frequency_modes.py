"""Signed-frequency contraction, budgets and the full-operator acceptance check."""

import unittest
from contextlib import ExitStack
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from frozen_sigma_response import herm
from lapw_frequency_modes import prepare_frequency_modes
from lapw_soc_screened import SOCScreenedForce
from response_parallel import SERIAL
from stress_response_adjoint import PairedOperator
from test_frozen_sigma_tangent import complex_model


class FrequencyModeTests(unittest.TestCase):
    def test_signed_complex_adjoint_including_mixed_compression(self):
        model = complex_model()
        model.precontract = False
        weights = model.meta.weights.copy()
        weights[[3, 8, 21]] *= -1
        weights[-1] = 0
        model.meta = replace(model.meta, weights=weights)
        model.prepare()
        rng = np.random.default_rng(304)
        y = herm(rng.normal(size=(2, 5, 5)) + 1j * rng.normal(size=(2, 5, 5)))
        with ExitStack() as resources:
            tables = prepare_frequency_modes(model, 1e-7, 2**24, resources)
            self.assertEqual(set(tables), {0, 1})
            self.assertIs(tables[0].stream, tables[1].stream)
            self.assertLess(tables[0].offset, tables[1].offset)
            for table in tables.values():
                self.assertIn(-1, table.signs)
                self.assertLessEqual(table.error, 1e-7)
            for selection in [tables, {1: tables[1]}]:
                for fixed_n in (False, True):
                    actual = model.adjoint(y, fixed_number=fixed_n, return_raw=True,
                                           frequency_modes=selection)
                    expected = model.adjoint(y, fixed_number=fixed_n, return_raw=True)
                    for got, ref in zip(actual, expected):
                        np.testing.assert_allclose(got, ref, atol=2e-12, rtol=2e-12)
        self.assertTrue(all(t.stream.closed for t in tables.values()))

    def test_insufficient_resources_or_rank_uses_full_frequencies(self):
        model = complex_model()
        model.precontract = False
        model.prepare()
        with ExitStack() as resources:
            self.assertFalse(prepare_frequency_modes(model, 1e-6, 0, resources))
            # Room for a single mode, insufficient for this metallic response.
            model.cache_limit = 16 * 40
            self.assertFalse(prepare_frequency_modes(model, 1e-7, 2**24, resources))

    def test_disk_budget_only_needs_retained_modes(self):
        model = complex_model()
        model.precontract = False
        model.prepare()
        w, q = len(model.W), len(model.Q)
        width = w*w + q + 2*w*model.uc.shape[-1]
        full_grid_bytes = len(model.meta.weights) * width * 16
        with ExitStack() as resources:
            budget = full_grid_bytes - width*16
            tables = prepare_frequency_modes(model, 1e-6, budget, resources)
            self.assertEqual(set(tables), {0, 1})
            import os
            self.assertLessEqual(os.fstat(tables[0].stream.fileno()).st_size, budget)

    def test_independent_residual_retries_with_full_frequencies(self):
        # Force both fallback branches. The reduced operator gives a false
        # solution when J=0, or cannot converge when J=I. The true J and its
        # fixed point are known analytically in this control.
        for approximate_j in (0.0, 1.0):
            with self.subTest(approximate_j=approximate_j), ExitStack() as resources:
                f = SOCScreenedForce.__new__(SOCScreenedForce)
                f.parallel = SERIAL
                f.progress = lambda line: None
                f.records = [SimpleNamespace(index=1)]
                f.scratch_limit, f.cache_scratch, f.n = 1024, 0, 2
                f.model = None
                identity = PairedOperator(2, 2, lambda x: x, lambda x: x)
                zero = PairedOperator(2, 2, np.zeros_like, np.zeros_like)
                f.pc = SimpleNamespace(potential_operator=identity, core=SimpleNamespace(operator=zero))
                f.norm = SimpleNamespace(operator=identity)
                f.core = SimpleNamespace(potential_gradient=np.zeros(2))
                f.den = SimpleNamespace(density_gradient=np.array([1., 2.]))
                true_j = np.diag([0.2, -0.3])

                def readout(y=None, *, geometry_file=None, frequency_modes=None):
                    if y is None:
                        return np.zeros(2)
                    return approximate_j * y if frequency_modes else true_j @ y

                f.potential_readout = readout
                modes = {0: SimpleNamespace(signs=[1], error=1e-7)}
                with patch("lapw_soc_screened.prepare_frequency_modes", return_value=modes):
                    sol, report, _, _, _ = f._screened_solution(
                        None, None, resources, 1, 1e-8, 2, 1e-6)
                self.assertTrue(report["full_frequency_fallback"])
                self.assertTrue(report["full_frequency_check"])
                self.assertTrue(report["converged"])
                np.testing.assert_allclose(sol, np.linalg.solve(np.eye(2)-true_j,
                                                               f.den.density_gradient), atol=1e-12)


if __name__ == "__main__":
    unittest.main()
