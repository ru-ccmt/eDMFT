"""Scheduling and SCF assembly for the explicit and screened spinor forces."""

from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import io
import json
import unittest
from unittest.mock import patch
import numpy as np
import dmft2_step as step
from lapw_response_prepare import working_directory
from stress_response_publish import force_replacement_rows


class SOCForceIntegrationTests(unittest.TestCase):
    def test_direct_force_follows_native_force_only(self):
        for mode, phase, final, option, expected in [
            ("FOR", "defer", False, None, True),
            ("FOR", "standalone", False, False, False),
            ("TOT", "defer", True, None, False),
            ("TOT", "final", True, None, True),
            ("FORT", "defer", False, None, False),
            ("FORT", "final", False, None, True),
            ("PRES", "defer", False, None, False),
            ("PRES", "final", False, None, True),
        ]:
            with self.subTest(
                mode=mode, phase=phase, final=final, option=option
            ), ExitStack() as ctx:
                work = ctx.enter_context(TemporaryDirectory())
                ctx.enter_context(working_directory(work))
                Path("test.in2").write_text(mode + "\n")
                Path("test.indmf2").write_text("c 0.0 1\n")
                Path("dmft2.def").write_text(
                    "2, 'test.indmf2', 'old', 'formatted', 0\n"
                    "5, 'test.in2', 'old', 'formatted', 0\n"
                    "9, 'test.vectorso', 'old', 'unformatted', 0\n"
                )
                if option is not None:
                    Path("test.stress_response.json").write_text(
                        json.dumps({"solver": {"soc_operator_force": option}})
                    )
                ctx.enter_context(
                    patch.object(step, "soc_operator_scope", return_value=None)
                )
                ctx.enter_context(
                    patch.object(
                        step,
                        "preflight",
                        return_value={
                            "all_forces_fixed": False,
                            "physical_policy": {"magnetic": False},
                        },
                    )
                )
                core = ctx.enter_context(patch.object(step, "run_core"))
                calc = ctx.enter_context(
                    patch.object(step, "compute_response", return_value={})
                )
                ctx.enter_context(patch("stress_response_summary.finish"))
                step.execute(
                    "test", ".", phase=phase, final_force=final, output=io.StringIO()
                )
                self.assertEqual(calc.called, expected)
                self.assertEqual(
                    "extra_env" in core.call_args.kwargs, expected and mode not in {"FORT", "PRES"}
                )
                if expected:
                    self.assertEqual(
                        "captured_soc" in calc.call_args.kwargs, mode not in {"FORT", "PRES"}
                    )
                self.assertEqual(Path("test.in2").read_text().strip(), mode)

    def test_operator_addition_and_total_force_replacement_in_site_frame(self):
        # Global z maps to local x. Native force components total 4 along x.
        structure = SimpleNamespace(
            nat=1, rotloc=np.array([[[0, 0, 1], [1, 0, 0], [0, 1, 0]]])
        )
        with patch("wstruct.Struct", return_value=structure), patch(
            "force_response_modes.check_force_frame"
        ):
            structure.ReadStruct = lambda fn: None
            native = "".join(
                f":{tag}001:   1.ATOM{1.:18.9f}{1.:18.9f}{0.:18.9f}{0.:18.9f}\n"
                for tag in ("FHF", "FVA", "FSU", "FCO")
            )
            report = dict(
                structure="test.struct",
                components=[{"atom": 1}],
                selected_atoms=[1],
                force_mRy_Bohr=[[0.0, 0.0, -0.25]],
            )
            for additive, delta in [(False, -4.25), (True, -0.25)]:
                if additive:
                    report["application"] = "add_soc_operator"
                rows = force_replacement_rows(native, report, reference_id="test")
                x = [r for r in rows.splitlines() if r.startswith(":FVA")][0].split()
                np.testing.assert_allclose([float(t) for t in x[3:6]], [delta, 0, 0])
                with self.assertRaisesRegex(ValueError, "already applied"):
                    force_replacement_rows(native + rows, report, reference_id="test")
