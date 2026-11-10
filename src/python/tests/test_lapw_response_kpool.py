"""Native pool stream preservation and failure checks; no WIEN license needed."""

from pathlib import Path
import struct
import tempfile
import unittest
from lapw_response_kpool import klist_records, merge_streams, remap_matrix


def record(payload):
    marker = struct.pack("<i", len(payload))
    return marker + payload + marker


class PoolTests(unittest.TestCase):
    def test_klist_preserves_global_weight_and_energy_limits(self):
        for width in (5, 10):
            with self.subTest(width=width), tempfile.TemporaryDirectory() as folder:
                rows = [
                    f"{k:10d}"
                    + "".join(f"{v:{width}d}" for v in (1, 3, 5, 20))
                    + f"{weight:5.1f}"
                    for k, weight in [(1, 8), (2, 4)]
                ]
                suffix = " -7.0  1.5"
                path = Path(folder) / "case.klist"
                path.write_text(rows[0] + suffix + "\n" + rows[1] + "\nEND\n")
                parsed, end, limits = klist_records(path)
                self.assertEqual(end, 15 + 4 * width)
                self.assertEqual(parsed[1][:end] + limits, rows[1] + suffix)

    def test_vector_header_is_kept_once_without_transforming_coefficients(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            header = record(struct.pack("<2d", -0.3, 0.5)) + record(
                struct.pack("<d", 995.0)
            )
            paths = []
            for i in range(2):
                p = root / str(i)
                p.write_bytes(header + record(bytes([i]) * 50))
                paths.append(p)
            out = root / "vector"
            merge_streams(paths, out, 1, binary=True)
            self.assertEqual(
                out.read_bytes(), header + record(bytes(50)) + record(bytes([1]) * 50)
            )
            paths[1].write_bytes(
                record(struct.pack("<2d", -0.4, 0.5))
                + record(struct.pack("<d", 995.0))
                + record(b"x")
            )
            previous = out.read_bytes()
            with self.assertRaisesRegex(ValueError, "different radial"):
                merge_streams(paths, out, 1, binary=True)
            self.assertEqual(out.read_bytes(), previous)
            self.assertFalse(Path(str(out) + ".merging").exists())

    def test_energy_header_and_k_order(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            paths = []
            for k in (1, 2):
                p = root / str(k)
                p.write_bytes(b"energy\nlo\n" + f"k {k}\n".encode())
                paths.append(p)
            out = root / "energy"
            merge_streams(paths, out, 1, binary=False)
            self.assertEqual(out.read_bytes(), b"energy\nlo\nk 1\nk 2\n")

    def test_remap_changes_only_integer_global_provenance(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            p = root / "local"
            out = root / "global"
            payload = (
                b"EDMFTHV1"
                + struct.pack("<I8i", 0x01020304, 1, 2, 0, 10, 1, 2, 1, 1)
                + b"unchanged matrix payload"
            )
            p.write_bytes(payload)
            remap_matrix(p, out, 75, 3)
            self.assertEqual(out.read_bytes()[:16], payload[:16])
            self.assertEqual(out.read_bytes()[24:], payload[24:])
            self.assertEqual(struct.unpack("<ii", out.read_bytes()[16:24]), (75, 3))
