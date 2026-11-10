"""Portable native-context archive: physical policy, provenance and paired maps."""

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import json, tempfile, unittest
import numpy as np
from lapw_potential_core_state import (
    save_potential_core_state,
    load_potential_core_state,
    fingerprint_files,
)
from wien_response_density import Density


def inputs():
    r = np.geomspace(1e-4, 1.3, 101)
    d = Density.__new__(Density)
    d.lines = []
    d.spans = {}
    d.pw_lines = {}
    d.valence_format = False
    d.grids = [np.log(r)]
    d.mult = [1]
    d.rad = {(0, 0, 0): 4 * np.pi * r * r * 0.03}
    d.pw = {(0, 0, 0): 0.03 + 0j}
    s = SimpleNamespace(
        nat=1,
        lattyp="P",
        mode="NREL",
        a=5.0,
        b=6.0,
        c=7.0,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        mult=np.array([1]),
        jrj=np.array([len(r)]),
        r0=np.array([r[0]]),
        rmt=np.array([r[-1]]),
        pos=np.array([[0.1, 0.2, 0.3]]),
        timat=np.eye(3, dtype=int)[None],
        tau=np.zeros((1, 3)),
        rotloc=np.eye(3)[None],
        iatnr=np.array([-1]),
        Znuc=np.array([1.0]),
    )
    return dict(
        density=d,
        structure=s,
        fft_shape=(8, 8, 8),
        inversion=False,
        functional="LDA",
        angular_data=[
            dict(
                angular=np.array([[1 / np.sqrt(4 * np.pi)]]),
                weights=np.array([4 * np.pi]),
            )
        ],
        core_data=[dict(vsp=-2 * np.ones(len(r)), shift_Ry=0.0, states=[])],
    )


class StateTests(unittest.TestCase):
    def test_roundtrip_and_exact_transposes_without_native_launch(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "central.npz"
            args = inputs()
            central = save_potential_core_state(path, **args)
            loaded = load_potential_core_state(path, expected_fingerprints={})
            np.testing.assert_array_equal(
                loaded.input_density.vector(), args["density"].vector()
            )
            rng = np.random.default_rng(8)
            x = rng.normal(size=central.hartree.size)
            np.testing.assert_array_equal(
                loaded.potential_operator.matvec(x),
                central.potential_operator.matvec(x),
            )
            self.assertLess(loaded.potential_operator.check_transpose(samples=2), 1e-13)
            np.testing.assert_array_equal(loaded.core.action(x), x * 0)
            self.assertLess(loaded.core.operator.check_transpose(samples=1), 1e-14)
            normal = loaded.normalization(loaded.input_density, 6.3)
            self.assertLess(abs((normal.mt + normal.ir) @ normal.action(x)), 1e-12)

    def test_source_fingerprints_reject_stale_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / "clmsum"
            source.write_text("central")
            files = {"density": source}
            path = tmp / "central.npz"
            save_potential_core_state(path, source_files=files, **inputs())
            load_potential_core_state(
                path, expected_fingerprints=fingerprint_files(files)
            )
            source.write_text("new iteration")
            with self.assertRaisesRegex(ValueError, "fingerprints differ"):
                load_potential_core_state(
                    path, expected_fingerprints=fingerprint_files(files)
                )

    def test_failed_save_keeps_previous_archive(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "central.npz"
            args = inputs()
            save_potential_core_state(path, **args)
            before = path.read_bytes()
            args["angular_data"][0]["weights"][0] = float("nan")
            with self.assertRaises(ValueError):
                save_potential_core_state(path, **args)
            self.assertEqual(before, path.read_bytes())
            self.assertEqual([p.name for p in Path(tmp).iterdir()], ["central.npz"])

    def test_rejects_unsupported_conventions_and_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "central.npz"
            save_potential_core_state(path, **inputs())
            with np.load(path, allow_pickle=False) as z:
                base = {k: z[k].copy() for k in z.files}
            for field, value, error in [
                ("version", 2, ValueError),
                ("mode", "UNKNOWN", NotImplementedError),
                ("inversion", True, ValueError),
                ("functional", "SCAN", NotImplementedError),
                ("fft_shape", [8.1, 8, 8], ValueError),
                ("policy", {}, NotImplementedError),
            ]:
                data = deepcopy(base)
                meta = json.loads(str(data["metadata"]))
                meta[field] = value
                data["metadata"] = np.array(json.dumps(meta))
                np.savez(path, **data)
                with self.subTest(field=field), self.assertRaises(error):
                    load_potential_core_state(path)
            data = deepcopy(base)
            del data["angular_0"]
            np.savez(path, **data)
            with self.assertRaisesRegex(ValueError, "arrays"):
                load_potential_core_state(path)
            data = deepcopy(base)
            data["density"] = data["density"].astype(complex)
            np.savez(path, **data)
            with self.assertRaisesRegex(ValueError, "real coordinates"):
                load_potential_core_state(path)
            data = deepcopy(base)
            data["metadata"] = np.array("[]")
            np.savez(path, **data)
            with self.assertRaisesRegex(ValueError, "JSON object"):
                load_potential_core_state(path)

    def test_rejects_fractional_fft_shape_before_writing(self):
        args = inputs()
        args["fft_shape"] = (8.1, 8, 8)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "central.npz"
            with self.assertRaisesRegex(ValueError, "integer FFT"):
                save_potential_core_state(path, **args)
            self.assertFalse(path.exists())

    def test_rejects_density_from_different_radial_grid(self):
        args = inputs()
        args["density"].grids[0][-1] += 0.01
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "radial grids differ"):
                save_potential_core_state(Path(tmp) / "central.npz", **args)
