"""Partition ownership and shared adjoint source contractions, independent matrices."""

from types import SimpleNamespace as NS
from pathlib import Path
import hashlib
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from response_parallel import ResponseParallel, SERIAL
from stress_response_adjoint import PairedOperator
from lapw_strain_charge import PreparedChargeSourceContraction


class ParallelTests(unittest.TestCase):
    def test_partition_covers_once_and_balances_cubic_cost(self):
        costs = np.arange(1, 20, dtype=float) ** 3
        ownership = []
        loads = []
        for rank in range(4):
            p = ResponseParallel(NS(Get_rank=lambda: rank, Get_size=lambda: 4))
            ids = p.partition(costs)
            ownership.extend(ids.tolist())
            loads.append(costs[ids].sum())
        self.assertEqual(sorted(ownership), list(range(len(costs))))
        self.assertLess(max(loads) - min(loads), costs.max())
        with self.assertRaisesRegex(ValueError, "exceed"):
            p.partition([1, 1])

    def test_fingerprints_identical_to_native_sha256(self):
        with tempfile.TemporaryDirectory() as folder:
            f = Path(folder) / "vector"
            f.write_bytes(b"native basis" * 500000)
            self.assertEqual(
                SERIAL.fingerprints({"v": f}),
                {"v": hashlib.sha256(f.read_bytes()).hexdigest()},
            )

    def test_shared_adjoint_replaces_full_charge_source_exactly(self):
        rng = np.random.default_rng(733)

        def op(a):
            return PairedOperator(
                a.shape[1], a.shape[0], lambda x: a @ x, lambda y: a.T @ y
            )

        B = rng.normal(size=(7, 4))
        E = rng.normal(size=(8, 7))
        R = rng.normal(size=(5, 8))
        D = rng.normal(size=(5, 4))
        C = rng.normal(size=(5, 4))
        N = rng.normal(size=(5, 5))
        state = NS(
            jacobian=NS(input_size=5),
            chain=NS(
                basis=op(B),
                electronic=op(E),
                density=op(R),
                basis_density=op(D),
                core=op(C),
                normalization=op(N),
            ),
        )
        lam = rng.normal(size=5)
        readout = PreparedChargeSourceContraction(state, lam)
        # The actual native normalization/transport offset is tested separately;
        # this test verifies that its contribution is added once, outside N.
        offset = rng.normal(size=5)
        with patch(
            "lapw_strain_charge.material_source_from_raw", return_value=(offset, {})
        ):
            for unused in range(6):
                b = rng.normal(size=7)
                v = rng.normal(size=4)
                g = rng.normal(size=5)
                expected = lam @ (
                    N @ (R @ E @ (b + B @ v) + D @ v + C @ v + g) + offset
                )
                actual = readout.contraction(b, g, v, np.eye(3))
                self.assertAlmostEqual(actual, expected, places=11)
