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

    def test_energy_only_never_responds(self):
        with case_directory(charge="m"), patch.object(step, "run_core"), patch.object(
            step, "compute_response"
        ) as response:
            result = step.execute("sample", "/unused", phase="final")
            self.assertEqual(result["status"], "skipped_noncharge")
            response.assert_not_called()

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

    def test_concurrent_run_cannot_overwrite_status(self):
        with case_directory(), patch.object(step, "run_core") as core:
            with open(".sample.dmft2-response.lock", "a+") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                Path("sample.stress_response.status.json").write_text(
                    '{"status":"running_response"}'
                )
                with self.assertRaises(BlockingIOError):
                    step.execute("sample", "/unused")
                self.assertEqual(
                    Path("sample.stress_response.status.json").read_text(),
                    '{"status":"running_response"}',
                )
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

    def test_preflight_rejects_fixed_mu_and_magnetism_before_core(self):
        with case_directory():
            req = step.request("sample")
            Path("sample.indmf2").write_text("c 1.25 0 1.0 .T. 0.02\n")
            with self.assertRaisesRegex(ValueError, "recomputeEF=1"):
                step.preflight("sample", req, "/unused", "")
            Path("sample.indmfldn").write_text("magnetic")
            with self.assertRaisesRegex(ValueError, "nonmagnetic"):
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
                "dmft2_step": SimpleNamespace(requested=lambda case: enabled),
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
            arguments = [io.StringIO(), "sample", io.StringIO(), None, None, None]
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
