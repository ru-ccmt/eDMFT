"""Boundary tests for isolated analytic preparation and total publication."""

from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import io, json, os, tempfile, unittest
import dmft2_step as step


@contextmanager
def fixture():
    old = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        os.chdir(d)
        try:
            root = Path(d).resolve()
            bin = root / "bin"
            wi = root / "wien"
            bin.mkdir()
            wi.mkdir()
            for base, names in [
                (
                    bin,
                    [
                        "dmft2",
                        "lapw0_response",
                        "lapw1_response",
                        "lcore_response",
                        "liblapw_response.so",
                    ],
                ),
                (wi, ["x", "lapw0", "lapw1", "lapw1c", "lcore"]),
            ]:
                for name in names:
                    p = base / name
                    p.write_text("#!/bin/sh\nexit 0\n")
                    p.chmod(0o755)
            for name, text in {
                "sample.indmf2": "c 1.25 1 1.0 .T. 0.02\n",
                "sample.in2": "STRT\n",
                "sample.in0": "STR XC_PBE\n",
                "sample.indmfl": "1 4 1 5\n",
                "sample.struct": "test structure",
                "sample.clmsum": "test charge",
                "EF.dat": "1.25",
                "projectorw.dat": "projector",
                "lapw1.def": "5, 'sample.in1', 'old', 'formatted', 0\n",
                "sample.in1": "WFFIL\n",
                "sample.inm": "MSR1 0.0 YES\n",
            }.items():
                Path(name).write_text(text)
            files = {
                2: "sample.indmf2",
                5: "sample.in2",
                7: "sample.indmfl",
                20: "sample.struct",
            }
            Path("dmft2.def").write_text(
                "".join(
                    f"{u}, '{f}', 'old', 'formatted', 0\n" for u, f in files.items()
                )
            )
            req = {
                "config": {"solver": {"wien_root": str(wi)}},
                "files": files,
                "definition": "dmft2.def",
                "mode": "STRT",
            }
            yield root, bin, wi, req
        finally:
            os.chdir(old)


class IntegrationTests(unittest.TestCase):
    def test_analytic_inherits_mpi_but_explicit_serial_still_wins(self):
        with fixture() as (root, bin, wi, req):
            structure = SimpleNamespace(
                a=4., b=5., c=6., lattyp="P", alpha=90.0, beta=90.0, gamma=90.0, ReadStruct=lambda _: None
            )
            with patch("wstruct.Struct", return_value=structure), patch(
                "lapw_total_stress.invariant_strains"
            ), patch("response_mpi.automatic_python_launch",
                     return_value=["python-mpiexec", "-n", "4"]):
                result = step.preflight("sample", req, bin, "mpirun -np 4")
                self.assertEqual(result["backend"], "analytic")
                self.assertEqual(result["mpi"], ["python-mpiexec", "-n", "4"])
                self.assertEqual(result["native_mpi"], ["mpirun", "-np", "4"])
                self.assertEqual(result["dmft2"], str(bin / "dmft2"))
                self.assertEqual(result["lapw1_worker"], str(bin / "lapw1_response"))
                req["config"]["solver"]["mpi"] = []
                result = step.preflight("sample", req, bin, "mpirun -np 4")
                self.assertEqual(result["mpi"], [])
                self.assertEqual(result["native_mpi"], ["mpirun", "-np", "4"])

    def test_response_executable_override(self):
        with fixture() as (_, bin, _, req):
            structure = SimpleNamespace(
                a=4., b=5., c=6., lattyp="P", alpha=90., beta=90., gamma=90.,
                ReadStruct=lambda _: None,
            )
            # An older installation can contain a separate response binary.
            alternate = bin / "dmft2_response"
            alternate.touch(mode=0o755)
            with patch("wstruct.Struct", return_value=structure), patch(
                "lapw_total_stress.invariant_strains"
            ):
                self.assertEqual(step.preflight("sample", req, bin, "")["dmft2"], str(bin / "dmft2"))
                req["config"]["solver"]["dmft2"] = str(alternate)
                self.assertEqual(step.preflight("sample", req, bin, "")["dmft2"], str(alternate))
                alternate.unlink()
                with self.assertRaisesRegex(ValueError, "missing executable"):
                    step.preflight("sample", req, bin, "")

    def test_screened_force_method_defaults_and_displaced_override(self):
        with fixture() as (_, bin, _, req):
            structure = SimpleNamespace(
                a=4., b=5., c=6., lattyp="P", alpha=90., beta=90., gamma=90.,
                nat=1, ReadStruct=lambda _: None,
            )
            with patch("wstruct.Struct", return_value=structure), patch(
                "lapw_total_stress.invariant_strains"
            ), patch("force_response_modes.internal_modes", return_value=[object()]), patch(
                "force_response_modes.check_force_frame"
            ):
                for mode in ("FORT", "STRFT"):
                    req["mode"] = mode
                    for source in (None, "displaced"):
                        with self.subTest(mode=mode, source=source):
                            solver = req["config"]["solver"]
                            if source is None:
                                solver.pop("force_source", None)
                            else:
                                solver["force_source"] = source
                            result = step.preflight("sample", req, bin, "")
                            self.assertEqual(result["force_source"], source or "analytic")
                            self.assertEqual(result["backend"], "analytic")

    def test_stress_ensemble_is_independent_of_central_fermi_search(self):
        with fixture() as (_, bin, _, req):
            structure = SimpleNamespace(
                a=4., b=5., c=6., lattyp="P", alpha=90., beta=90., gamma=90.,
                nat=1, ReadStruct=lambda _: None,
            )
            Path("sample.inm").write_text("MSR1 0.0 NO\n")
            settings = req["config"]["solver"]
            with patch("wstruct.Struct", return_value=structure), patch(
                "lapw_total_stress.invariant_strains"
            ), patch("force_response_modes.internal_modes", return_value=[]), patch(
                "force_response_modes.check_force_frame"
            ):
                for search, ensemble, legacy, expected in (
                    (0, None, None, "fixed-n"), (1, None, None, "fixed-n"),
                    (0, "fixed-mu", None, "fixed-mu"), (1, "fixed-mu", None, "fixed-mu"),
                    (0, None, False, "fixed-mu"), (0, "fixed-n", False, "fixed-n"),
                ):
                    with self.subTest(search=search, ensemble=ensemble, legacy=legacy):
                        Path("sample.indmf2").write_text(f"c 1.25 {search} 1.0 .T. 0.02\n")
                        if ensemble:
                            settings["stress_ensemble"] = ensemble
                        else:
                            settings.pop("stress_ensemble", None)
                        if legacy is not None:
                            settings["physical_policy"] = {"fixed_number": legacy}
                        else:
                            settings.pop("physical_policy", None)
                        for mode in ("STRT", "STRFT"):
                            req["mode"] = mode
                            result = step.preflight("sample", req, bin, "")
                            self.assertEqual(result["physical_policy"]["fixed_number"], bool(search))
                            self.assertFalse(result["physical_policy"]["normalization_enabled"])
                            self.assertEqual(result["stress_ensemble"], expected)
                        self.assertEqual(Path("EF.dat").read_text(), "1.25")
                req["mode"] = "FORT"
                settings.pop("stress_ensemble")
                settings.pop("physical_policy")
                self.assertNotIn("stress_ensemble", step.preflight("sample", req, bin, ""))
                req["mode"] = "STRT"
                settings["stress_ensemble"] = "fixed-ne"
                with self.assertRaisesRegex(ValueError, "stress_ensemble"):
                    step.preflight("sample", req, bin, "")

    def test_soc_stress_uses_saved_derivative_constraint_and_replays_older_snapshots(self):
        import lapw_soc_screened as soc

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            argv = ["lapw_soc_screened.py", "--case", "sample", "--snapshot", directory,
                    "--pressure", "--output", str(root / "out")]
            for central, derivative in ((False, True), (True, False), (False, None)):
                with self.subTest(central=central, derivative=derivative):
                    snapshot = dict(case="sample", spinor=True, physical_policy={"fixed_number": central})
                    if derivative is not None:
                        snapshot["response_policy"] = {"fixed_number": derivative}
                    (root / "preparation.json").write_text(json.dumps(snapshot))
                    with patch("sys.argv", argv), patch.object(soc, "SOCScreenedForce") as model:
                        model.return_value.run.return_value = {}
                        soc.main()
                    expected = central if derivative is None else derivative
                    self.assertEqual(model.call_args.kwargs["fixed_number"], expected)

    def test_preparation_and_driver_are_isolated_then_ready_not_applied(self):
        with fixture() as (root, bin, wi, req):
            before = {p: p.read_bytes() for p in root.iterdir() if p.is_file()}
            solver = dict(
                backend="analytic", stress_ensemble="fixed-n",
                physical_policy=dict(magnetic=False,ferro=False,fixed_number=False,normalization_enabled=True),
                wien_root=str(wi),
                dmft2=str(bin / "dmft2"),
                library=str(bin / "liblapw_response.so"),
            )
            for label in ["lapw0", "lapw1", "lcore"]:
                solver[label + "_worker"] = str(bin / (label + "_response"))
            calls = []

            def child(command, **kw):
                calls.append(command)
                self.assertEqual(kw["env"]["OPENBLAS_NUM_THREADS"], "1")
                script = Path(command[1]).name
                if script == "lapw_response_prepare.py":
                    cfg = json.loads(Path(command[2]).read_text())
                    self.assertFalse(cfg["physical_policy"]["fixed_number"])
                    self.assertEqual(cfg["stress_ensemble"], "fixed-n")
                    self.assertEqual(cfg["dmft2"], str(bin / "dmft2"))
                    self.assertEqual(Path(cfg["reference"]).name, "reference")
                    self.assertEqual(Path(cfg["workspace"]).name, "central")
                    self.assertNotEqual(Path(cfg["reference"]), root)
                else:
                    self.assertEqual(script, "lapw_total_stress_driver.py")
                    self.assertEqual(Path(command[command.index("--progress-file") + 1]),
                                     root / "sample.outputdmf2_response")
                    out = Path(command[command.index("--output") + 1])
                    out.mkdir()
                    (out / "report.json").write_text(
                        json.dumps(
                            {
                                "converged": True,
                                "adjoint": {"converged": True},
                                "full_tensor_within_symmetry": True,
                                "stress_Ry_Bohr3": [
                                    [0.1, 0, 0],
                                    [0, 0.1, 0],
                                    [0, 0, 0.2],
                                ],
                            }
                        )
                    )
                return SimpleNamespace(returncode=0)

            with patch.object(step, "preflight", return_value=solver), patch.object(
                step.subprocess, "run", side_effect=child
            ):
                result = step.compute_response(
                    "sample", req, bin, "", "run", io.StringIO()
                )
            self.assertEqual(len(calls), 2)
            self.assertTrue(result["stress_ready"])
            self.assertFalse(result["stress_applied"])
            self.assertNotIn("lapw1.def", result["parent_fingerprints"])
            self.assertIn("sample.clmsum", result["parent_fingerprints"])
            for p, data in before.items():
                self.assertEqual(p.read_bytes(), data)

    def test_failed_native_preparation_does_not_start_driver(self):
        with fixture() as (root, bin, wi, req):
            solver = dict(
                backend="analytic", stress_ensemble="fixed-n",
                physical_policy=dict(magnetic=False,ferro=False,fixed_number=True,normalization_enabled=True),
                wien_root=str(wi),
                dmft2=str(bin / "dmft2"),
                library=str(bin / "liblapw_response.so"),
            )
            for label in ["lapw0", "lapw1", "lcore"]:
                solver[label + "_worker"] = str(bin / (label + "_response"))
            with patch.object(step, "preflight", return_value=solver), patch.object(
                step.subprocess, "run", return_value=SimpleNamespace(returncode=1)
            ) as child:
                with self.assertRaisesRegex(RuntimeError, "prepare.log"):
                    step.compute_response("sample", req, bin, "", "run", io.StringIO())
                child.assert_called_once()


if __name__ == "__main__":
    unittest.main()


class MPIIntegrationTests(unittest.TestCase):
    def test_separate_python_native_launches_and_cache_budgets(self):
        with fixture() as (root, bin, wi, req):
            solver = dict(
                backend="analytic", stress_ensemble="fixed-n",
                physical_policy=dict(magnetic=False,ferro=False,fixed_number=True,normalization_enabled=True),
                wien_root=str(wi),
                dmft2=str(bin / "dmft2"),
                library=str(bin / "liblapw_response.so"),
                mpi=["python-mpiexec", "-n", "4"],
                native_mpi=["native-mpirun", "-np", "2"],
                cache_mb=64.0,
                basis_cache_mb=32.0, disk_matching=True,
            )
            for label in ["lapw0", "lapw1", "lcore"]:
                solver[label + "_worker"] = str(bin / (label + "_response"))
            calls = []

            def child(command, **kw):
                calls.append(command)
                if len(calls) == 1:
                    cfg = json.loads(Path(command[-1]).read_text())
                    self.assertEqual(cfg["native_mpi"], solver["native_mpi"])
                    self.assertNotIn("python-mpiexec", command)
                else:
                    self.assertEqual(command[:3], solver["mpi"])
                    self.assertIn("--mpi", command)
                    self.assertIn("--disk-matching", command)
                    self.assertEqual(command[command.index("--cache-mb") + 1], "64.0")
                    self.assertEqual(
                        command[command.index("--basis-cache-mb") + 1], "32.0"
                    )
                    out = Path(command[command.index("--output") + 1])
                    out.mkdir()
                    (out / "report.json").write_text(
                        json.dumps(
                            dict(
                                converged=True, adjoint={"converged": True},
                                full_tensor_within_symmetry=True,
                                stress_Ry_Bohr3=[[0.1, 0, 0], [0, 0.1, 0], [0, 0, 0.1]],
                            )
                        )
                    )
                return SimpleNamespace(returncode=0)

            with patch.object(step, "preflight", return_value=solver), patch.object(
                step.subprocess, "run", side_effect=child
            ):
                result = step.compute_response(
                    "sample", req, bin, "production-mpi", "run", io.StringIO()
                )
            self.assertTrue(result["stress_ready"])
            self.assertEqual(len(calls), 2)
