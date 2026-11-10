"""Snapshot policy, lifecycle and stale-state rejection controls."""

from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch
import unittest
import numpy as np
from lapw_valence_state import load_native_valence_state, _prepare_snapshot
from lapw_potential_core_state import fingerprint_files
from test_frozen_sigma_tangent import complex_model


class ValenceStateTests(unittest.TestCase):
    def _files(self, path):
        names = [
            "case" + s
            for s in [
                ".struct",
                ".indmfl",
                ".vector",
                ".vsp",
                ".vns",
                ".clmval",
                ".in1",
            ]
        ]
        names += [
            "projectorw.dat",
            "response.meta",
            "response.k00000001.bin",
            "native.k00000001.bin",
            "native.ham.k00000001.bin",
        ]
        for name in names:
            (path / name).write_bytes(b"policy-test placeholder")
        return {name: path / name for name in names}

    def test_native_matrix_reader_endian_and_length(self):
        from lapw_native_snapshot import read_basis_matrix

        rng = np.random.default_rng(42)
        m = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))
        m = (m + m.conj().T) / 2
        labels = np.arange(15).reshape(3, 5)
        k = np.array([0.1, 0.2, 0.3])
        packed = np.concatenate([m[: j + 1, j] for j in range(5)])
        with TemporaryDirectory() as d:
            p = Path(d) / "matrix.bin"
            for endian in ["<", ">"]:
                with p.open("wb") as f:
                    f.write(b"EDMFTHV1")
                    np.array([0x01020304], endian + "i4").tofile(f)
                    np.array([1, 7, 0, 5, 2, 2, 1, 1], endian + "i4").tofile(f)
                    np.array(k, endian + "f8").tofile(f)
                    np.array(labels, endian + "i4").ravel(order="F").tofile(f)
                    np.array(packed, endian + "c16").tofile(f)
                actual = read_basis_matrix(p, b"EDMFTHV1")
                np.testing.assert_array_equal(actual.matrix, m)
                np.testing.assert_array_equal(actual.labels, labels)
                np.testing.assert_array_equal(actual.reciprocal_k, k)
                self.assertEqual((actual.kpoint, actual.nlo), (7, 2))
                with self.assertRaises(ValueError):
                    read_basis_matrix(p, b"EDMFTOV1")
                with p.open("ab") as f:
                    f.write(b"x")
                with self.assertRaisesRegex(ValueError, "length"):
                    read_basis_matrix(p, b"EDMFTHV1")

    def test_explicit_policy_and_missing_snapshot(self):
        with TemporaryDirectory() as d:
            for kwargs, error in [
                (
                    {"energy_policy": "automatic", "complex_vectors": False},
                    NotImplementedError,
                ),
                ({"energy_policy": "fixed", "complex_vectors": None}, ValueError),
                ({"energy_policy": "fixed", "complex_vectors": False}, ValueError),
            ]:
                with self.assertRaises(error):
                    load_native_valence_state(
                        d, "case", matrix_quadrature="wien2k24", **kwargs
                    )

    def test_fingerprint_check_precedes_any_native_preparation(self):
        with TemporaryDirectory() as d:
            p = Path(d)
            paths = self._files(p)
            expected = fingerprint_files(paths)
            (p / "case.vsp").write_bytes(b"changed potential")
            with patch("lapw_valence_state._prepare_snapshot") as prepare:
                with self.assertRaisesRegex(ValueError, "fingerprints"):
                    load_native_valence_state(
                        p,
                        "case",
                        complex_vectors=False,
                        energy_policy="fixed",
                        matrix_quadrature="wien2k24",
                        expected_fingerprints=expected,
                    )
                prepare.assert_not_called()

    def test_preparation_failure_releases_all_contexts(self):
        closed = []

        def fail(*args, **kwargs):
            args[-1].callback(lambda: closed.append(True))
            raise ValueError("deliberate downstream failure")

        with TemporaryDirectory() as d:
            self._files(Path(d))
            with patch("lapw_valence_state._prepare_snapshot", side_effect=fail):
                with self.assertRaisesRegex(ValueError, "downstream"):
                    load_native_valence_state(
                        d,
                        "case",
                        complex_vectors=False,
                        energy_policy="fixed",
                        matrix_quadrature="wien2k24",
                    )
        self.assertEqual(closed, [True])

    def test_stale_global_normalization_is_rejected(self):
        # Exported normalized overlaps remain mutually consistent with X;
        # only X's dependence on the global raw-overlap sum is deliberately stale.
        m = complex_model()
        m.meta.normalization = m.meta.normalization * 1.03
        m.X = m.meta.normalization.copy()
        given = m.apply_x(m.raw, m.X)
        records = [
            SimpleNamespace(
                energies=e,
                coefficients=np.eye(5)[:, :, None],
                labels=np.zeros((3, 5), int),
                nlo=0,
                projector=given[k],
            )
            for k, e in enumerate(m.eps)
        ]
        with patch(
            "lapw_valence_state._snapshot_inputs",
            return_value=SimpleNamespace(s=None, nlo_total=0),
        ), patch(
            "lapw_valence_state.CachedFrozenSigmaResponse.from_export",
            return_value=(m, records),
        ):
            from contextlib import ExitStack

            with ExitStack() as stack, self.assertRaisesRegex(
                ValueError, "normalization is stale"
            ):
                _prepare_snapshot(
                    "case",
                    Path("/does-not-need-to-exist"),
                    5,
                    False,
                    None,
                    False,
                    stack,
                )


if __name__ == "__main__":
    unittest.main()
