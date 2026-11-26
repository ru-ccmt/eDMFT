"""WIEN's energy cutoff can give a different stored band count at each k."""
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from frozen_sigma_response import FrozenSigmaResponse, Metadata
from lapw_dft_vector import load_vectors
from response_parallel import SERIAL
from test_frozen_sigma_response import model


def write_export(prefix, k, n, window=(2, 3)):
    """One native DMFT2 export. Bands are numbered from one in the file."""
    nbasis = 6
    path = prefix.with_suffix(f".k{k:08d}.bin")
    with path.open("wb") as fh:
        fh.write(b"EDMRPK01")
        np.array([1, k, nbasis, 1, n, 1, n, *window, 1, 1, 0], "<i4").tofile(fh)
        np.array([1.0], "<f8").tofile(fh)  # each of the two k points has weight 1/2
        np.arange(n, dtype="<f8").tofile(fh)
        np.eye(nbasis, n, dtype="<c16").ravel(order="F").tofile(fh)
        np.arange(n, dtype="<c16").tofile(fh)
        np.zeros((nbasis, 3), "<i4").tofile(fh)


class VariableBandTests(unittest.TestCase):
    def test_common_explicit_block_preserves_window_and_original_records(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(Metadata, "read", return_value=model().meta):
            prefix = Path(directory) / "response"
            write_export(prefix, 1, 5)
            write_export(prefix, 2, 4)
            result, records = FrozenSigmaResponse.from_export(prefix, lapw_complement=True)
            self.assertEqual(result.eps.shape, (2, 4))
            np.testing.assert_array_equal(result.W, [1, 2])
            self.assertEqual([len(r.energies) for r in records], [5, 4])
            with self.assertRaisesRegex(ValueError, "require the native LAPW"):
                FrozenSigmaResponse.from_export(prefix)
            with self.assertRaisesRegex(ValueError, "misses DMFT window"):
                FrozenSigmaResponse.from_export(prefix, lapw_complement=True, retained=2)
            write_export(prefix, 2, 4, window=(2, 4))
            with self.assertRaisesRegex(ValueError, "fixed-index DMFT window"):
                FrozenSigmaResponse.from_export(prefix, lapw_complement=True)

    def test_dft_vector_keeps_each_k_points_own_band_count(self):
        def record(fh, data):
            size = struct.pack("<i", len(data))
            fh.write(size + data + size)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.vector"
            with path.open("wb") as fh:
                for _ in range(2):  # one atom, two radial header records
                    record(fh, b"radial")
                for nbands in [4, 3]:
                    record(fh, struct.pack("<ddd10siid", 0, 0, 0, b"K", 6, nbands, 1))
                    record(fh, np.zeros(18, "<i4").tobytes())
                    for i in range(nbands):
                        record(fh, struct.pack("<id", i + 1, float(i)))
                        record(fh, np.eye(6, dtype="<f8")[i].tobytes())
            records, weights = load_vectors(path, 1, 0, SERIAL)
            self.assertEqual([r.nstored for r in records], [4, 3])
            np.testing.assert_array_equal(weights, [0.5, 0.5])
