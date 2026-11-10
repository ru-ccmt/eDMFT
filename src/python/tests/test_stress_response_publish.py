"""Complete-total replacement: signs, native fixed columns and stale inputs."""

from pathlib import Path
import hashlib, json, tempfile, unittest
import numpy as np
from stress_response_publish import (
    replacement_rows,
    response_rows,
    mixer_accumulator,
    mark_published,
)


def fixture():
    text = ""
    for name, values in [
        ("VAL", [0.3, 0.4, 0.5]),
        ("CORE", [17.0, 18.0, 19.0]),
        ("DENS", [17.4, 18.6, 19.8]),
    ]:
        for i, row in enumerate(np.diag(values)):
            text += (
                f":STR_{name}{i+1:03d}:  " + "".join(f"{x:20.10f}" for x in row) + "\n"
            )
    return text


class PublicationTests(unittest.TestCase):
    def test_complete_tensor_replaces_old_result_with_native_sign(self):
        old = fixture()
        sigma = np.diag([-0.001, 0.002, 0.003])
        rows = replacement_rows(old, sigma, reference_id="test")
        np.testing.assert_allclose(-mixer_accumulator(old + rows), sigma, atol=3e-15)
        # Fortran SCFANA reads VAL using column 15, three F20.10 values.
        lines = [l for l in rows.splitlines() if l.startswith(":STR_VAL")]
        native = np.array(
            [
                [float(line[14 + 20 * j : 14 + 20 * (j + 1)]) for j in range(3)]
                for line in lines
            ]
        )
        np.testing.assert_allclose(
            -(mixer_accumulator(old) + native), sigma, atol=3e-15
        )

    def test_missing_core_or_duplicate_publication_rejected(self):
        text = fixture()
        sigma = np.eye(3) * 0.001
        with self.assertRaisesRegex(ValueError, "complete DENS"):
            replacement_rows(
                "\n".join(l for l in text.splitlines() if "CORE" not in l),
                sigma,
                reference_id="x",
            )
        with self.assertRaisesRegex(ValueError, "already applied"):
            replacement_rows(
                text + replacement_rows(text, sigma, reference_id="x"),
                sigma,
                reference_id="x",
            )

    def test_snapshot_and_report_must_match_then_publication_is_recorded(self):
        with tempfile.TemporaryDirectory() as folder:
            d = Path(folder)
            source = d / "sample.clmsum"
            source.write_text("central charge")
            report = d / "tensor.json"
            report.write_text(
                json.dumps(
                    dict(
                        converged=True,
                        full_tensor_within_symmetry=True,
                        stress_Ry_Bohr3=(np.eye(3) * 0.001).tolist(),
                    )
                )
            )
            sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
            status = dict(
                backend="analytic",
                stress_ready=True,
                status="complete",
                converged=True,
                reference_id="r",
                report=str(report),
                report_sha256=sha(report),
                parent_fingerprints={source.name: sha(source)},
            )
            path = d / "sample.stress_response.status.json"
            path.write_text(json.dumps(status))
            rows = response_rows("sample", fixture(), directory=d)
            self.assertIn(":STR_VAL001:", rows)
            mark_published("sample", rows, d / "sample.scf", directory=d)
            self.assertTrue(json.loads(path.read_text())["stress_applied"])
            source.write_text("different charge")
            with self.assertRaisesRegex(ValueError, "stale"):
                response_rows("sample", fixture(), directory=d)
            source.write_text("central charge")
            report.write_text("{}")
            with self.assertRaisesRegex(ValueError, "report changed"):
                response_rows("sample", fixture(), directory=d)

    def test_default_omission_can_be_published_with_valid_report(self):
        # The tested production prescription is not an explicit sensitivity
        # diagnostic. Publication still requires valid fingerprints/convergence.
        with tempfile.TemporaryDirectory() as folder:
            d = Path(folder)
            report = d / 'tensor.json'
            sigma = np.eye(3)*.001
            report.write_text(json.dumps(dict(
                converged=True, full_tensor_within_symmetry=True,
                diagnostic_only=False, stress_ready=True,
                charge_fixed_point=dict(residual_norm=.00503141,
                    warning_threshold=.001, within_threshold=False,
                    warning='Central charge-map residual exceeds warning threshold.'),
                window_gap_treatment=[dict(policy='omit_singular_pairs',
                    cutoff_Ry=1e-10, omitted_pairs=9, affected_kpoints=8)],
                stress_Ry_Bohr3=sigma.tolist())))
            status = dict(backend='analytic', stress_ready=True,
                status='complete', converged=True, reference_id='default-omission',
                report=str(report), parent_fingerprints={},
                report_sha256=hashlib.sha256(report.read_bytes()).hexdigest())
            (d / 'sample.stress_response.status.json').write_text(json.dumps(status))
            rows = response_rows('sample', fixture(), directory=d)
            np.testing.assert_allclose(-mixer_accumulator(fixture()+rows), sigma, atol=1e-14)

    def test_window_omission_report_cannot_reach_mixer(self):
        # Even a converged diagnostic carrying a full tensor and a valid hash
        # must not bypass the explicit window-regularity publication policy.
        with tempfile.TemporaryDirectory() as folder:
            d = Path(folder)
            report = d / 'tensor.json'
            report.write_text(json.dumps(dict(
                converged=True, full_tensor_within_symmetry=True,
                diagnostic_only=True, stress_Ry_Bohr3=np.eye(3).tolist())))
            status = dict(backend='analytic', stress_ready=True,
                status='complete', converged=True, reference_id='diagnostic',
                report=str(report), parent_fingerprints={},
                report_sha256=hashlib.sha256(report.read_bytes()).hexdigest())
            (d / 'sample.stress_response.status.json').write_text(json.dumps(status))
            with self.assertRaisesRegex(ValueError, 'complete supported tensor'):
                response_rows('sample', fixture(), directory=d)

    def test_no_tensor_for_ordinary_or_deferred_iteration(self):
        with tempfile.TemporaryDirectory() as folder:
            d = Path(folder)
            self.assertEqual(response_rows("sample", "", directory=d), "")
            (d / "sample.stress_response.status.json").write_text(
                json.dumps(dict(status="deferred_until_final"))
            )
            self.assertEqual(response_rows("sample", "", directory=d), "")


if __name__ == "__main__":
    unittest.main()
