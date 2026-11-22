"""Scalar, spinor, energy and norms have to describe the same per-k spaces."""

from pathlib import Path
from tempfile import TemporaryDirectory
import struct
import unittest

from lapw_soc_rank import _put, enforce_ranks, ranks, parent_ranks, parent_lo_filter


class RankTests(unittest.TestCase):
    def test_parent_matching_policy(self):
        with TemporaryDirectory() as temp:
            p = Path(temp) / "Ce.scfso"
            self.assertTrue(parent_lo_filter(p))
            p.write_text("unmarked old spectrum\n")
            self.assertTrue(parent_lo_filter(p))
            for marker, expected in (
                ("LINEAR", False),
                ("LEGACY_COMPONENT_CUTOFF", True),
            ):
                p.write_text(":SOC_LO_MATCHING: " + marker + "\n")
                self.assertEqual(parent_lo_filter(p), expected)
            p.write_text(":SOC_LO_MATCHING: UNKNOWN\n")
            with self.assertRaisesRegex(ValueError, "unknown"):
                parent_lo_filter(p)

    def make_vector(self, folder, name, counts, spacing=0.25):
        vector, energy, norm = [
            folder / (name + ext) for ext in (".vector", ".energy", ".norm")
        ]
        with vector.open("wb") as v, energy.open("w") as e, norm.open("w") as n:
            for _ in range(2):
                _put(v, b"radial")
                e.write("radial\n")
            for index, count in enumerate(counts):
                k = (0.25 * index, 0.0, 0.0)
                _put(v, struct.pack("<ddd10siid", *k, b"K", 3, count, 1.0))
                e.write(
                    "".join(f"{q:19.12f}" for q in k)
                    + "K".ljust(10)
                    + f"{3:6}{count:6}{1.:5.1f}\n"
                )
                _put(v, struct.pack("<9i", *range(9)))
                for b in range(1, count + 1):
                    _put(v, struct.pack("<id", b, spacing * b))
                    _put(v, b"\0" * 48)
                    e.write(f"{b} {spacing*b}\n")
                for i in range(0, count, 4):
                    n.write(
                        "".join(f"{i+j:20.9f}" for j in range(min(4, count - i))) + "\n"
                    )
        return vector, energy, norm

    def test_actual_rank_at_each_k_and_norm_blocks(self):
        with TemporaryDirectory() as temp:
            d = Path(temp)
            scalar = self.make_vector(d, "sc", (4, 3))
            up = self.make_vector(d, "up", (2, 5))
            down = self.make_vector(d, "dn", (2, 5))
            manifest = parent_ranks(scalar[0], up[0], down[0], 1)
            pool = self.make_vector(d, "pool", (8, 8))
            enforce_ranks(pool[0], pool[1], 1, manifest["spinor"], pool[2])
            self.assertEqual(ranks(pool[0], 1), manifest["spinor"])
            self.assertEqual(pool[1].read_text(), up[1].read_text())
            self.assertEqual(pool[2].read_text(), up[2].read_text())
            with self.assertRaisesRegex(ValueError, "too small"):
                enforce_ranks(pool[0], pool[1], 1, manifest["scalar"], pool[2])

    def test_signed_reference_checks_cut_before_replacing_files(self):
        with TemporaryDirectory() as temp:
            d = Path(temp)
            target = ranks(self.make_vector(d, "target", (2,))[0], 1)
            v, e, n = self.make_vector(d, "pool", (4,), spacing=0)
            original = [f.read_bytes() for f in (v, e, n)]
            with self.assertRaisesRegex(ValueError, "degenerate"):
                enforce_ranks(v, e, 1, target, n, check_cut=True)
            self.assertEqual(original, [f.read_bytes() for f in (v, e, n)])
            enforce_ranks(v, e, 1, target, n)
            self.assertEqual(ranks(v, 1), target)

    def test_mismatched_spinors_are_rejected(self):
        with TemporaryDirectory() as temp:
            d = Path(temp)
            one = self.make_vector(d, "one", (4, 4))[0]
            two = self.make_vector(d, "two", (4, 3))[0]
            with self.assertRaisesRegex(ValueError, "meshes differ"):
                parent_ranks(one, one, two, 1)
