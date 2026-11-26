"""Execution-order and input-preservation tests for the response integration."""

import ast
from contextlib import contextmanager
import fcntl
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import dmft2_step as step


@contextmanager
def case_directory(mode="STRT", charge="c"):
    old = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        os.chdir(d)
        try:
            Path("sample.in2").write_text(mode + "\n -10 20 0.0\n")
            Path("sample.indmf2").write_text(
                "# header\n" + charge + " 1.25 1 1.0 .T. 0.02\n"
            )
            Path("dmft2.def").write_text(
                " 2, 'sample.indmf2', 'old', 'formatted', 0\n 5, 'sample.in2', 'old', 'formatted', 0\n"
            )
            yield Path(d)
        finally:
            os.chdir(old)


class ExecutionTests(unittest.TestCase):
    def test_ordinary_run_does_not_respond(self):
        with case_directory("STR"), patch.object(
            step, "run_core"
        ) as core, patch.object(step, "compute_response") as response:
            result = step.execute("sample", "/unused")
            self.assertEqual(result["status"], "disabled")
            core.assert_called_once()
            response.assert_not_called()
            self.assertFalse(Path("sample.stress_response.status.json").exists())

    def test_live_request_and_defer_until_final(self):
        with case_directory(), patch.object(step, "run_core") as core, patch.object(
            step, "preflight"
        ), patch.object(
            step, "compute_response", return_value={"stress_applied": False}
        ) as response:
            for _ in range(3):
                self.assertEqual(
                    step.execute("sample", "/unused", phase="defer")["status"],
                    "deferred_until_final",
                )
            response.assert_not_called()
            self.assertIn("Response deferred until the final", Path("sample.outputdmf2_response").read_text())
            result = step.execute("sample", "/unused", phase="final")
            self.assertEqual(result["status"], "complete")
            response.assert_called_once()
            self.assertEqual(core.call_count, 4)
            # Read the controls anew at execution time; never use a cached flag.
            Path("sample.stress_response.json").write_text('{"enabled": false}')
            self.assertEqual(
                step.execute("sample", "/unused", phase="final")["status"], "disabled"
            )
            response.assert_called_once()
            self.assertEqual(
                json.loads(Path("sample.stress_response.status.json").read_text())[
                    "status"
                ],
                "disabled",
            )

    def test_native_force_modes_follow_live_input(self):
        # Include both ordinary and screened input modes. Only the native
        # charge calculation runs on intermediate iterations.
        for requested, expected in (("FOR", "FOR"), ("FORT", "FOR"),
                                    ("STRF", "STRF"), ("STRFT", "STRF"),
                                    ("STRT", "STR"), ("TOT", "TOT")):
            with self.subTest(requested), case_directory(requested):
                modes = []
                def inspect(root, mpi, definition, output):
                    modes.append(step.input_mode(step.entries(Path(definition).read_text())[5]))
                with patch.object(step, "run_core", side_effect=inspect):
                    step.execute("sample", ".", phase="defer")
                    step.execute("sample", ".", phase="defer", final_force=True)
                    Path("sample.in2").write_text("FOR\n")
                    step.execute("sample", ".", phase="defer")
                self.assertEqual(modes, [expected, expected, "FOR"])

    def test_final_force_uses_private_input_only_for_charge(self):
        for mode, phase, expected in (("c", "final", "FOR"), ("c", "standalone", "TOT"),
                                      ("m", "final", "TOT")):
            with self.subTest(mode=mode, phase=phase), case_directory("TOT", mode):
                inp, definition = Path("sample.in2c"), Path("dmft2.def")
                Path("sample.in2").rename(inp)
                definition.write_text(definition.read_text().replace("sample.in2", inp.name))
                original = inp.read_bytes(), definition.read_bytes()
                def inspect(root, mpi, temp_definition, output):
                    current = step.entries(Path(temp_definition).read_text())[5]
                    self.assertEqual(step.input_mode(current), expected)
                    self.assertEqual(Path(current).read_text().splitlines()[1:],
                                     inp.read_text().splitlines()[1:])
                    raise RuntimeError("stop")
                with patch.object(step, "run_core", side_effect=inspect), self.assertRaisesRegex(RuntimeError, "stop"):
                    step.execute("sample", ".", phase=phase, final_force=True)
                self.assertEqual((inp.read_bytes(), definition.read_bytes()), original)
                self.assertFalse(list(Path(".").glob(".dmft2-response-input-*")))

    def test_named_final_measurements_from_tot(self):
        for value, native_mode, screened in [
            (True, "FOR", False), ("FOR", "FOR", False), ("FORT", "FOR", True),
            ("STRF", "STRF", False), ("STRFT", "STRF", True),
        ]:
            with self.subTest(value=value), case_directory("TOT"):
                native = []
                def inspect(root, mpi, definition, output):
                    native.append(step.input_mode(step.entries(Path(definition).read_text())[5]))
                with patch.object(step, "preflight"), patch.object(step, "run_core", side_effect=inspect), \
                     patch.object(step, "compute_response", return_value={}) as response:
                    step.execute("sample", ".", phase="defer", final_force=value)
                    result = step.execute("sample", ".", phase="final", final_force=value)
                self.assertEqual(native, ["TOT", native_mode])
                logical = "FOR" if value is True else value
                self.assertEqual(result["requested_mode"], logical)
                self.assertEqual(response.call_count, int(screened))
                if screened:
                    self.assertEqual(response.call_args.args[1]["mode"], logical)
                self.assertEqual(step.input_mode("sample.in2"), "TOT")

    def test_final_setting_validation_and_precedence(self):
        with self.assertRaisesRegex(ValueError, "force_last"):
            step.final_force_mode("STRTFT")
        self.assertEqual(step.final_force_mode(" fort "), "FORT")
        with case_directory("FOR"):
            self.assertEqual(step.request("sample", final_force="STRFT")["mode"], "FOR")
        with case_directory("TOT", charge="m"):
            self.assertEqual(step.request("sample", final_force="FORT")["mode"], "TOT")
        with case_directory("TOT"):
            Path("sample.stress_response.json").write_text('{"enabled":false}')
            req = step.request("sample", final_force="STRFT")
            self.assertEqual(req["mode"], "STRFT")
            self.assertFalse(req["enabled"])

    def test_energy_only_never_responds(self):
        with case_directory(charge="m"), patch.object(step, "run_core"), patch.object(
            step, "compute_response"
        ) as response:
            result = step.execute("sample", "/unused", phase="final")
            self.assertEqual(result["status"], "skipped_noncharge")
            response.assert_not_called()

    def test_final_native_force_for_scalar_and_magnetic_runs(self):
        import run_dmft as driver
        from contextlib import ExitStack
        from types import SimpleNamespace
        with case_directory("TOT"):
            Path("sample.inm").write_text("MSR1\n")
            for ext in ("", "dn"):
                with self.subTest(ext), ExitStack() as stack:
                    calls = []
                    def record(name, **kwargs):
                        calls.append(name)
                        if name in ("dmft2", "lattice_charge_step"):
                            self.assertEqual(kwargs["final_force"], expected_mode)
                    for name in ("dmft2", "lattice_charge_step", "lcore", "scf", "scf1", "mixer", "scfm"):
                        stack.enter_context(patch.object(driver, name,
                            side_effect=lambda *args, n=name, **kwargs: record(n, **kwargs)))
                    params = {"force_last": False}
                    args = (io.StringIO(), SimpleNamespace(case="sample", WIENROOT="."),
                            io.StringIO(), None, params, "", {"m_extn": ext}, True)
                    self.assertFalse(driver.final_response_pass(*args))
                    for value, expected_mode in ((True, "FOR"), ("FORT", "FORT"), ("STRFT", "STRFT")):
                        calls.clear()
                        params["force_last"] = value
                        self.assertTrue(driver.final_response_pass(*args))
                        # A screened magnetic measurement owns both channels.
                        first = "lattice_charge_step" if ext and value is True else "dmft2"
                        self.assertEqual(calls, [first, "lcore", "scf", "scf1", "mixer", "scfm"])
            Path("sample.inm").write_text("MSR1a\n")
            self.assertFalse(driver.final_response_pass(*args))
            Path("sample.inm").write_text("MSR1\n")
            Path("sample.in2").write_text("FOR\n")
            self.assertFalse(driver.final_response_pass(*args))

    def test_scheduled_lapw0_stress_tracks_active_in2_and_restores_input(self):
        import run_dmft as driver
        for suffix in ('in2', 'in2c'):
            for mode in ('FOR', 'FORT', 'STRF', 'STRFT', True, False):
                with self.subTest(suffix=suffix, mode=mode), case_directory('TOT'):
                    if suffix == 'in2c':
                        Path('sample.in2').rename('sample.in2c')
                        d = Path('dmft2.def')
                        d.write_text(d.read_text().replace("sample.in2'", "sample.in2c'"))
                        Path('sample.in2').write_text('FOR\n')
                    Path('sample.inm').write_text('MSR1\n')
                    path = Path('sample.in0')
                    original = b'TOT  XC_PBESOL\r\n  NR2V\tIFFT\r\n'
                    path.write_bytes(original)
                    with driver.final_lapw0_input('sample', mode, io.StringIO()):
                        expected = b'STR' if mode in ('STRF', 'STRFT') else b'TOT'
                        self.assertEqual(path.read_bytes(), expected + original[3:])
                        self.assertEqual(step.input_mode('sample.' + suffix), 'TOT')
                    self.assertEqual(path.read_bytes(), original)

    def test_final_lapw0_failure_and_live_edits(self):
        import run_dmft as driver
        with case_directory('TOT'):
            Path('sample.inm').write_text('MSR1\n')
            path, shared = Path('sample.in0'), Path('shared.in0')
            shared.write_text('TOT XC_PBE\n')
            path.symlink_to(shared)
            with self.assertRaisesRegex(RuntimeError, 'failure'):
                with driver.final_lapw0_input('sample', 'STRFT', io.StringIO()):
                    self.assertEqual(shared.read_text(), 'TOT XC_PBE\n')
                    self.assertEqual(path.read_text(), 'STR XC_PBE\n')
                    raise RuntimeError('failure')
            self.assertTrue(path.is_symlink())
            self.assertEqual(path.readlink(), shared)
            path.unlink()
            path.write_text('TOT XC_PBE\n')
            with driver.final_lapw0_input('sample', 'STRF', io.StringIO()):
                path.write_text('STR XC_PBESOL\n')
            self.assertEqual(path.read_text(), 'TOT XC_PBESOL\n')

    def test_final_lapw0_respects_explicit_inputs(self):
        import run_dmft as driver
        for in0, in2, inm, charge in [('STR', 'TOT', 'MSR1', 'c'), ('TOT', 'FOR', 'MSR1', 'c'),
                                       ('TOT', 'TOT', 'MSR1a', 'c'), ('TOT', 'TOT', 'MSR1', 'm')]:
            with self.subTest(in0=in0, in2=in2, inm=inm, charge=charge), case_directory(in2, charge):
                Path('sample.inm').write_text(inm + '\n')
                path = Path('sample.in0')
                path.write_text(in0 + ' XC_PBESOL\n')
                with driver.final_lapw0_input('sample', 'STRFT', io.StringIO()):
                    self.assertEqual(step.input_mode(path), in0)
                self.assertEqual(step.input_mode(path), in0)

    def test_mode_aliases_preserve_prepared_files(self):
        for alias, native in [("STRT", "STR"), ("STRFT", "STRF")]:
            with self.subTest(alias=alias), case_directory(alias):
                original = {
                    p: p.read_bytes()
                    for p in [
                        Path("sample.in2"),
                        Path("sample.indmf2"),
                        Path("dmft2.def"),
                    ]
                }
                with step.execution_definition("dmft2.def") as definition:
                    self.assertNotEqual(definition, "dmft2.def")
                    name = step.entries(Path(definition).read_text())[5]
                    self.assertEqual(step.input_mode(name), native)
                    self.assertEqual(
                        Path(name).read_text().splitlines()[1], " -10 20 0.0"
                    )
                self.assertFalse(Path(definition).exists())
                for p, b in original.items():
                    self.assertEqual(p.read_bytes(), b)

    def test_temporary_inputs_removed_on_failure(self):
        with case_directory(), patch.object(
            step, "run_core", side_effect=RuntimeError("core failed")
        ), patch.object(step, "preflight"), patch.object(
            step, "compute_response"
        ) as response:
            with self.assertRaisesRegex(RuntimeError, "core failed"):
                step.execute("sample", "/unused")
            response.assert_not_called()
            self.assertFalse(list(Path(".").glob(".dmft2-response-input-*")))
            self.assertEqual(step.input_mode("sample.in2"), "STRT")
            self.assertEqual(
                json.loads(Path("sample.stress_response.status.json").read_text())[
                    "status"
                ],
                "failed",
            )

    def test_response_failure_invalidates_previous_success(self):
        with case_directory(), patch.object(step, "run_core"), patch.object(
            step, "preflight"
        ), patch.object(
            step, "compute_response", side_effect=RuntimeError("no convergence")
        ):
            Path("sample.stress_response.status.json").write_text(
                '{"status":"complete","outputs":{"old":"old"}}'
            )
            with self.assertRaisesRegex(RuntimeError, "no convergence"):
                step.execute("sample", "/unused")
            result = json.loads(Path("sample.stress_response.status.json").read_text())
            self.assertEqual(result["status"], "failed")
            self.assertNotIn("outputs", result)
            summary = Path("sample.outputdmf2_response").read_text()
            self.assertIn("FAILED", summary)
            self.assertIn("no convergence", summary)

    def test_concurrent_run_cannot_overwrite_status(self):
        with case_directory(), patch.object(step, "run_core") as core:
            with open(".sample.dmft2-response.lock", "a+") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                Path("sample.stress_response.status.json").write_text(
                    '{"status":"running_response"}'
                )
                Path("sample.outputdmf2_response").write_text("active job")
                with self.assertRaises(BlockingIOError):
                    step.execute("sample", "/unused")
                self.assertEqual(
                    Path("sample.stress_response.status.json").read_text(),
                    '{"status":"running_response"}',
                )
                self.assertEqual(Path("sample.outputdmf2_response").read_text(), "active job")
            core.assert_not_called()

    def test_config_errors_are_not_silent(self):
        with case_directory():
            with self.assertRaisesRegex(ValueError, "does not exist"):
                step.request("sample", config_path="missing.json")
            Path("sample.stress_response.json").write_text('{"every":3}')
            with self.assertRaisesRegex(ValueError, "unknown response"):
                step.request("sample")
            Path("sample.stress_response.json").write_text('{"enabled":"false"}')
            with self.assertRaisesRegex(ValueError, "true or false"):
                step.request("sample")

    def test_def_localization_preserves_case_x_name(self):
        src = " 9, '/scratch/my_x_case.vector_x', 'old', 'unformatted', 9000\n 20, '/source/my_x_case.struct', 'old', 'formatted', 0\n"
        result = step.entries(step.local_definition(src))
        self.assertEqual(result[9], "my_x_case.vector")
        self.assertEqual(result[20], "my_x_case.struct")

    def test_actual_child_exit_checks(self):
        with case_directory("STR") as d:
            exe = d / "dmft2"
            for body, success in [
                ('echo "STOP DMFT2 END"', True),
                ('echo "STOP NUMERICAL ERROR"', False),
                ('echo "STOP DMFT2 END"; exit 4', False),
                ('echo error > dmft2.error; echo "STOP DMFT2 END"', False),
            ]:
                exe.write_text("#!/bin/sh\n" + body + "\n")
                exe.chmod(0o755)
                if success:
                    step.run_core(str(d), "", "dmft2.def", io.StringIO())
                else:
                    with self.assertRaises(RuntimeError):
                        step.run_core(str(d), "", "dmft2.def", io.StringIO())

    def test_run_dmft_has_one_final_call_outside_charge_loop(self):
        tree = ast.parse(Path(step.__file__).with_name("run_dmft.py").read_text())
        found = []
        passes = []

        def visit(node, loops):
            loops = loops + int(isinstance(node, (ast.For, ast.While)))
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == "dmft2"
            ):
                for kw in node.keywords:
                    if (
                        kw.arg == "response_phase"
                        and isinstance(kw.value, ast.Constant)
                        and kw.value.value == "final"
                    ):
                        found.append(loops)
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == "final_response_pass"
            ):
                passes.append(loops)
            for child in ast.iter_child_nodes(node):
                visit(child, loops)

        visit(tree, 0)
        self.assertEqual(found, [0])
        self.assertEqual(passes, [0])

    def test_preflight_still_rejects_spin_polarized_KS(self):
        with case_directory():
            req=step.request("sample")
            Path("sample.clmup").write_text("spin polarized KS")
            with self.assertRaisesRegex(ValueError,"unpolarized KS"):
                step.preflight("sample",req,"/unused","")

    def test_charge_rule_is_checked_before_starting_response(self):
        with case_directory():
            req = step.request("sample")
            req["config"]["solver"] = {"charge_integration": "incorrect"}
            with self.assertRaisesRegex(ValueError, "charge_integration"):
                step.preflight("sample", req, "/unused", "")
            for rule in ["native", "legacy_spline"]:
                req["config"]["solver"]["charge_integration"] = rule
                # The valid choice reaches the next missing input check.
                with self.assertRaisesRegex(ValueError, "lapw1.def"):
                    step.preflight("sample", req, "/unused", "")

    def test_final_sequence_and_failure_stop(self):
        from types import SimpleNamespace

        tree = ast.parse(Path(step.__file__).with_name("run_dmft.py").read_text())
        node = next(
            n
            for n in tree.body
            if isinstance(n, ast.FunctionDef) and n.name == "final_response_pass"
        )
        for enabled, fail in [(False, False), (True, False), (True, True)]:
            events = []

            def core(*args, **kw):
                events.append("dmft2:" + kw["response_phase"])
                if fail:
                    raise RuntimeError("response failed")

            scope = {
                "dmft2_step": SimpleNamespace(request=lambda case: {"enabled": enabled, "mode":"STRT"}),
                "dmft2": core,
            }
            for name in ["lcore", "scf", "scf1", "mixer", "scfm"]:
                scope[name] = lambda *args, n=name: events.append(n)
            exec(
                compile(
                    ast.Module(body=[node], type_ignores=[]), "run_dmft.py", "exec"
                ),
                scope,
            )
            arguments = [io.StringIO(), SimpleNamespace(case="sample", WIENROOT="."), io.StringIO(),
                         None, {"force_last": False}, "", {"m_extn": ""}, True]
            if fail:
                with self.assertRaisesRegex(RuntimeError, "response failed"):
                    scope["final_response_pass"](*arguments)
                self.assertEqual(events, ["dmft2:final"])
            else:
                self.assertEqual(scope["final_response_pass"](*arguments), enabled)
                self.assertEqual(
                    events,
                    (
                        ["dmft2:final", "lcore", "scf", "scf1", "mixer", "scfm"]
                        if enabled
                        else []
                    ),
                )


if __name__ == "__main__":
    unittest.main()
