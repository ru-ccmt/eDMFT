"""Regression tests for the opt-in internal coordinate and cell optimizer."""

import os
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import structure_relaxation as relax


@contextmanager
def private_case():
    with tempfile.TemporaryDirectory() as name:
        old = Path.cwd()
        os.chdir(name)
        try:
            yield Path.cwd()
        finally:
            os.chdir(old)


def cubic():
    return SimpleNamespace(
        lattyp="P",
        a=5.0,
        b=5.0,
        c=5.0,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        timat=np.array([np.eye(3, dtype=int), np.diag([-1, 1, -1])]),
        nat=1,
        rotloc=np.eye(3)[None],
        tau=np.zeros((2, 3)),
    )


class TestStructuralRelaxation(unittest.TestCase):
    def test_opt_in_defaults_and_limits(self):
        self.assertFalse(relax.options(relax.DEFAULTS)["relax"])
        for key, val in [
            ("relax", 1),
            ("relax_max_steps", 0),
            ("relax_stress_tol", np.nan),
            ("relax_bulk_modulus", -1),
        ]:
            with self.assertRaises(ValueError):
                relax.options(dict(relax.DEFAULTS, **{key: val}))

    def test_pressure_sign_shear_projection_and_cap(self):
        s = cubic()
        sigma = np.array([[-5.0, 0.3, 0.4], [0.3, -5.0, -0.2], [0.4, -0.2, -5.0]])
        residual = relax.stress_residual(s, sigma, 5.0)
        np.testing.assert_allclose(residual, [[0, 0, 0.4], [0, 0, 0], [0.4, 0, 0]], atol=1e-13)
        eps = relax.cell_step(-10 * np.eye(3), relax.DEFAULTS)
        np.testing.assert_allclose(eps, np.eye(3) * relax.DEFAULTS["relax_max_strain"])
        self.assertLess(relax.cell_step(residual, relax.DEFAULTS)[0, 2], 0)

    def test_elastic_secant_and_bad_curvature(self):
        old = np.diag([4.0, 5.0, 6.0])
        eps = np.diag([0.001, 0.002, -0.003])
        new = np.diag(np.exp(np.diag(eps))) @ old
        change = 50 * eps
        opts = relax.DEFAULTS
        prior = dict(
            cell=old,
            next_cell=new,
            stress=np.zeros((3, 3)),
            compliance=relax.elastic_compliance(opts),
            controls=[opts[k] for k in relax.ELASTIC_CONTROLS],
        )
        h, applied = relax.updated_compliance(new, change, opts, prior)
        self.assertTrue(applied)
        np.testing.assert_allclose(h @ change.ravel(), eps.ravel(), atol=2e-14)
        self.assertGreater(np.min(np.linalg.eigvalsh(h)), 0)
        self.assertFalse(relax.updated_compliance(new, -change, opts, prior)[1])

    def test_published_force_correction_in_local_axes(self):
        s = cubic()
        s.rotloc = np.array([[[0, 0, 1], [0, 1, 0], [-1, 0, 0]]])
        delta = np.array([[0.0, 0.0, -0.9512]])
        rows = "".join(f":{tag}001:   1.ATOM 0 0 0 0\n" for tag in ("FCO", "FHF", "FSU"))
        with patch(
            "force_response_modes.internal_modes",
            return_value=[(1, np.array([0.0, 0.0, 1.0]), None)],
        ):
            np.testing.assert_allclose(
                relax.native_forces(rows + relax.correction_rows(delta, s), s), delta
            )
            with self.assertRaisesRegex(ValueError, "incomplete"):
                relax.native_forces(rows, s)

    def test_restore_modes_after_failure_and_edit(self):
        with private_case() as root:
            saved = {
                "X.in2": b"TOT             other\ntail\n",
                "X.in0": b"TOT   XC_PBE\n",
                "X.inm": b"MSR1a   0.0 YES\n",
            }
            for f, v in saved.items():
                Path(f).write_bytes(v)
            with self.assertRaisesRegex(RuntimeError, "abort"):
                with relax.managed_headers("X", root):
                    self.assertEqual(Path("X.in2").read_text().split()[0], "FOR")
                    self.assertEqual(Path("X.in0").read_text().split()[1], "XC_PBE")
                    with Path("X.in2").open("ab") as f:
                        f.write(b"later\n")
                    raise RuntimeError("abort")
            for f, v in saved.items():
                self.assertEqual(Path(f).read_bytes(), v + (b"later\n" if f == "X.in2" else b""))
            self.assertFalse((root / "input_headers.json").exists())

    def test_half_applied_cell_transaction(self):
        with private_case() as root:
            staged = []
            for name in ("struct", "clmsum"):
                src, dst = Path("next." + name), Path("X." + name)
                src.write_bytes(b"new" + name.encode())
                dst.write_bytes(b"old")
                staged.append(
                    dict(
                        dest=str(dst),
                        src=str(src),
                        before=relax.digest(dst),
                        after=relax.digest(src),
                    )
                )
            relax.atomic_json(root / "cell_commit.json", staged)
            Path("X.struct").write_bytes(Path("next.struct").read_bytes())
            Path("X.clmsum").write_bytes(b"external change")
            with self.assertRaisesRegex(RuntimeError, "external change"):
                relax.commit_cell(root)
            Path("X.clmsum").write_bytes(b"old")
            relax.commit_cell(root)
            self.assertEqual(Path("X.clmsum").read_bytes(), Path("next.clmsum").read_bytes())
            self.assertFalse((root / "cell_commit.json").exists())

    def test_adjoint_transport_uses_ordered_harmonics_and_g_labels(self):
        with private_case() as root:
            s = cubic()
            s.Znuc = [1]
            s.mult = [1]
            s.jrj = [2]
            s.r0 = [0.01]
            s.rmt = [1.0]
            d = SimpleNamespace(rad={(0, 0, 0): np.zeros(2)}, pw={(0, 0, 0): 0, (1, 1, 0): 0})
            pc = SimpleNamespace(structure=s, input_density=d)
            seed = relax.AdjointSeed(pc, "scalar")
            seed.save(root, np.arange(6.0))
            f = root / "adjoint_seed.npz"
            d.pw = {(1, 1, 0): 0, (2, 0, 0): 0, (0, 0, 0): 0}
            np.testing.assert_equal(
                relax.AdjointSeed(pc, "scalar").load(f), [0, 1, 4, 5, 0, 0, 2, 3]
            )
            self.assertIsNone(relax.AdjointSeed(pc, "soc").load(f))
            f.write_bytes(b"short")
            self.assertIsNone(seed.load(f))

    def test_coordinate_records_must_agree(self):
        with private_case() as root:
            s = SimpleNamespace(pos=np.array([[0.2, 0.4, 0.3]]))
            p = root / "X.struct"
            p.write_text("head\nPrecise positions\n .2 .4 .3000000041\n")
            with patch.object(relax, "structure", return_value=s):
                relax.agree_on_positions(p)
                np.testing.assert_equal(
                    np.fromstring(p.read_text().splitlines()[-1], sep=" "), s.pos[0]
                )
                p.write_text("head\nPrecise positions\n .2 .4 .31\n")
                with self.assertRaisesRegex(ValueError, "inconsistent"):
                    relax.agree_on_positions(p)

    def test_first_force_line_shows_actual_charge_residual(self):
        with private_case() as root:
            block = relax.ElectronicBlock.__new__(relax.ElectronicBlock)
            block.case = "X"
            block.directory = root / "0001"
            block.inner_iteration = 1
            block.last_mode = "MSR1a"
            block.frozen = False
            block.force_kind = "FOR"
            block.peak = 32.13878
            block.predictor = None
            block.hold_for_refresh = False
            Path("X.scfm").write_text(":DIS  : CHARGE DISTANCE (0.3 for atom 1) 0.000832\n")
            params = {"cc": 2e-5, "ec": 1e-6}
            block.check_ordinary(1.0, 1.0, params)
            line = (root / "history").read_text().splitlines()[-1]
            self.assertRegex(
                line,
                r"^\[\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\] 0001\.001 MSR1a \(int\.-relax\) "
                r"FOR +32\.13878 mRy/Bohr, rho 8\.32e-4 +\(<=2e-5,no\)$",
            )
            block.inner_iteration = 2
            block.frozen = True
            block.last_mode = "MSR1"
            block.force_kind = "FOR+scrn"
            block.peak = 5.04273
            block.check_ordinary(params["cc"], 3e-6, params)
            self.assertTrue(block.charge_met)
            held = (root / "history").read_text().splitlines()[-1]
            self.assertTrue(" ".join(held.split()).endswith(
                "0001.002 MSR1 (posi.-held) FOR+scrn 5.04273 mRy/Bohr, rho 2.00e-5 (<=2e-5,yes)"
            ))
            for marker in ("(", "FOR", "mRy/Bohr", "rho", "(<="):
                self.assertEqual(line.index(marker), held.index(marker), marker)
            self.assertEqual(line.index(".13878"), held.index(".04273"))
            for charge, expected in [(9.4e-6, "9.40e-6"), (1e-12, "1.00e-12"), (0, "0.00e+0")]:
                block.check_ordinary(charge, 3e-6, params)
                row = (root / "history").read_text().splitlines()[-1]
                self.assertIn("rho " + expected, row)
                self.assertEqual(line.index("(<="), row.index("(<="))

    def test_only_charge_and_force_converged_structures_move(self):
        s = cubic()
        s.aname = ["A"]
        s.pos = np.zeros((1, 3))
        s.mult = [1]
        with private_case() as root:
            for suffix, content in [
                ("in2", "TOT  rest"),
                ("inm", "MSR1  data"),
                ("in0", "TOT  XC_PBE"),
                ("struct", "old geometry"),
                ("clmsum", "charge"),
            ]:
                Path("X." + suffix).write_text(content + "\n")

            def measured(cmd, env, **kw):
                step = Path(env["EDMFT_RELAX_BLOCK"])
                index = int(step.name)
                geo = step / "measured.struct"
                geo.write_bytes(Path("X.struct").read_bytes())
                relax.atomic_json(
                    step / "measurement.json",
                    dict(
                        structure=str(geo),
                        charge_distance=0.00001,
                        charge_converged=index > 1,
                        force_mRy_Bohr=[[0, 0, 1.0 if index == 2 else 0.0]],
                        stress_GPa=(np.diag([-0.164949, -0.164949, -0.440762])
                                    if index < 4 else -2 * np.eye(3)).tolist(),
                        correction=[[0, 0, 0]],
                    ),
                )

            control = dict(
                relax.DEFAULTS, relax=True, relax_pressure=2.0, relax_max_steps=1,
                finish=1, cc=0.0001, ec=0.00001
            )
            with patch.object(relax, "structure", return_value=s), patch.object(
                relax.subprocess, "run", side_effect=measured
            ), patch.object(relax, "advance") as move:
                for _ in range(3):
                    with self.assertRaisesRegex(RuntimeError, "relax_max_steps reached"):
                        relax.run("X", control, "script")
                    move.assert_not_called()
                self.assertTrue((root / "X.relax/pending_step.json").exists())
                relax.run("X", control, "script")
                move.assert_called_once()
                self.assertTrue(relax.read_json(root / "X.relax/status.json")["converged"])
                history = (root / "X.relax/history").read_text().splitlines()
                z = next(i for i, line in enumerate(history) if "0003 tensile stress z" in line)
                self.assertTrue(history[z + 1].endswith("0003 pressure (GPa) 0.256887"))
                self.assertTrue(any(line.endswith("0004 pressure (GPa) 2.000000") for line in history))
                self.assertEqual(Path("X.in2").read_text(), "TOT  rest\n")

    def test_mixer_history_survives_small_internal_moves(self):
        with private_case() as root:
            step = root / "X.relax" / "0001"
            step.mkdir(parents=True)
            Path("X.struct").write_text("test structure without precise positions")
            Path("X.inm").write_text("MSR1   0.0 YES\n")
            Path("X.scf").write_text("")
            other = Path("Another.broyd1")
            other.write_text("another case")
            b = relax.ElectronicBlock("X", step, relax.DEFAULTS)
            files = [Path("X.broyd1"), Path("X.broyd2001")]

            def mix():
                for f in files:
                    f.write_text("history")

            with patch.object(relax, "structure", return_value=cubic()), patch.object(
                relax, "native_forces"
            ) as forces:
                forces.return_value = np.array([[0.0, 0.0, 2.0]])
                b.mix(0, mix)
                b.mix(0, lambda: self.assertTrue(all(f.exists() for f in files)))
                forces.return_value = np.zeros((1, 3))
                b.mix(0, lambda: self.assertFalse(any(f.exists() for f in files)))
            self.assertTrue(other.is_file())
            self.assertIn("2 files", (step.parent / "history").read_text())

    def test_current_overlap_cache_roundtrip(self):
        from lapw_soc_electronic import saved_normalization_error

        data = np.array([[[1.0], [0.2 + 0.1j]], [[0.2 - 0.1j], [1.4]]])
        meta = SimpleNamespace(dimensions=[2], normalization=data)
        with private_case() as root:
            step = root / "X.relax" / "0001"
            step.mkdir(parents=True)
            Path("SOlapm.dat").write_bytes(b"old projector")
            relax.refresh_projector(step)
            self.assertFalse(Path("SOlapm.dat").exists())
            self.assertEqual(next(step.glob("SOlapm.input.*")).read_bytes(), b"old projector")
            with patch("frozen_sigma_response.Metadata.read", return_value=meta):
                relax.save_measured_normalization(step, step, ("",))
            self.assertLess(saved_normalization_error(step / "SOlapm.measured.dat", meta), 1e-14)

    def test_small_moving_force_waits_for_charge(self):
        with private_case() as root:
            step = root / "X.relax" / "0001"
            step.mkdir(parents=True)
            Path("X.struct").write_text("no precision tail")
            Path("X.inm").write_text("MSR1   0.0 YES\n")
            Path("X.scf").write_text("")
            b = relax.ElectronicBlock("X", step, relax.DEFAULTS)
            with patch.object(relax, "structure", return_value=cubic()), patch.object(
                relax, "native_forces", return_value=np.array([[0, 0, 0.1]])
            ):
                for _ in range(4):
                    b.mix(0, lambda: None)
                    self.assertFalse(b.frozen)
                b.charge_met = True
                b.mix(0, lambda: None)
                self.assertTrue(b.frozen)


if __name__ == "__main__":
    unittest.main()
