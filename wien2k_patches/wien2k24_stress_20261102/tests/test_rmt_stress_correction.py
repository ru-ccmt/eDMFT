#!/usr/bin/env python3
"""Small, calculation-free regression tests for rmt_stress_correction.py."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).parents[2] / "src" / "python" / "rmt_stress_correction.py"
if not SCRIPT.is_file():
    # In a source checkout the patch bundle is below the repository root.
    SCRIPT = Path(__file__).parents[3] / "src" / "python" / "rmt_stress_correction.py"

spec = importlib.util.spec_from_file_location("rmt_stress_correction", SCRIPT)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot load {SCRIPT}")
rmt = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = rmt
spec.loader.exec_module(rmt)


class RMTRadiusSensitivityTests(unittest.TestCase):
    def test_scale_rkmax_preserves_rest_of_line(self) -> None:
        source = "WFFIL EF=0\n  7.00     10   4   ELPA comment\nrest\n"
        self.assertEqual(
            rmt.scaled_in1(source, 1.001),
            "WFFIL EF=0\n  7.00700     10   4   ELPA comment\nrest\n",
        )

    def test_scale_all_rmt_fields_without_shifting_z(self) -> None:
        source = (
            "Fe NPT= 781 R0=.000050000 RMT=   2.20000 Z:  26.00000\n"
            "Se NPT= 781 R0=.000050000 RMT=   2.10000 Z:  34.00000\n"
        )
        result = rmt.scaled_struct(source, 0.999)
        self.assertIn("RMT=   2.19780 Z:  26.00000", result)
        self.assertIn("RMT=   2.09790 Z:  34.00000", result)

    def test_diagnostic_cannot_modify_a_stress_tensor(self) -> None:
        self.assertFalse(hasattr(rmt, "corrected_tensor"))
        self.assertIn("must not be added", rmt.__doc__)

    @staticmethod
    def side_result(
        *, returncode: int, energies: tuple[float, ...], matrix_sizes: tuple[int, ...]
    ) -> rmt.SideResult:
        return rmt.SideResult(
            side="plus",
            scale=1.001,
            directory="/temporary/plus",
            energy_ry=energies[-1],
            energy_history_ry=energies,
            matrix_sizes=matrix_sizes,
            energy_records=len(energies),
            command_returncode=returncode,
            converged=returncode == 0,
            charge_distance=2.0e-6,
        )

    def test_early_convergence_before_maximum_is_success(self) -> None:
        result = self.side_result(
            returncode=0,
            energies=(-10.0, -10.1, -10.1001),
            matrix_sizes=(292, 292, 292),
        )
        with tempfile.TemporaryDirectory() as directory:
            # Three completed cycles are valid even though twenty were allowed.
            rmt.validate_side_convergence(result, 20, Path(directory) / "run.log")

    def test_iteration_limit_without_convergence_is_error(self) -> None:
        result = self.side_result(
            returncode=1,
            energies=(-10.0, -10.1, -10.1001),
            matrix_sizes=(292, 292, 292),
        )
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(rmt.RMTStressError, "did not converge"):
                rmt.validate_side_convergence(result, 3, Path(directory) / "run.log")

    def test_equal_fixed_basis_accepts_unequal_cycle_counts(self) -> None:
        self.assertEqual(
            rmt.common_matrix_size((292, 292, 292), (292, 292, 292, 292)),
            292,
        )

    def test_basis_change_is_rejected(self) -> None:
        self.assertIsNone(rmt.common_matrix_size((292, 293), (292, 292)))
        self.assertIsNone(rmt.common_matrix_size((292,), (293,)))

    def test_sensitivity_history_is_aligned_at_final_cycle(self) -> None:
        history = rmt.sensitivity_history_gpa(
            (-10.0, -10.1, -10.2, -10.3),
            (-9.0, -9.1),
            delta=0.001,
            volume=100.0,
        )
        expected = tuple(
            -((plus - minus) / 0.002) / 300.0 * rmt.RY_PER_BOHR3_TO_GPA
            for minus, plus in zip((-10.2, -10.3), (-9.0, -9.1))
        )
        self.assertEqual(history, expected)

    def test_lapw_command_requires_energy_and_charge_convergence(self) -> None:
        self.assertEqual(
            rmt.make_lapw_command(
                "/opt/wien/run_lapw", "SrVO3", 20, 1.0e-6, 1.0e-5
            ),
            [
                "/opt/wien/run_lapw",
                "-f",
                "SrVO3",
                "-i",
                "20",
                "-ec",
                "1e-06",
                "-cc",
                "1e-05",
            ],
        )

    def test_temporary_inputs_disable_expensive_stress_only(self) -> None:
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            (directory / "case.in0").write_text("STR  XC_PBE comment\nrest\n")
            (directory / "case.in2c").write_text("STRF   comment\nrest\n")
            (directory / "case.in2").write_text("TOT    comment\nrest\n")
            changed = rmt.select_energy_only_inputs(directory, "case")
            self.assertEqual(changed, ("case.in0", "case.in2c"))
            self.assertEqual(
                (directory / "case.in0").read_text(), "TOT  XC_PBE comment\nrest\n"
            )
            self.assertEqual(
                (directory / "case.in2c").read_text(), "TOT   comment\nrest\n"
            )
            self.assertEqual(
                (directory / "case.in2").read_text(), "TOT    comment\nrest\n"
            )


if __name__ == "__main__":
    unittest.main()
