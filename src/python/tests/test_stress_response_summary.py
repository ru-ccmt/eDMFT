"""Readable reports must match native SCF publication and preserve history."""
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch
import hashlib
import io
import json
import os
import tempfile
import unittest
import numpy as np
import merge_scf
import scf_merge
import stress_response_summary as summary
from stress_response_publish import mixer_accumulator, response_rows, mark_published
from test_stress_response_publish import fixture as tensors


@contextmanager
def case():
    old = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        os.chdir(d)
        try:
            root = Path(d)
            work = root/'sample.response_runs/test'
            work.mkdir(parents=True)
            report = dict(
                converged=True, full_tensor_within_symmetry=True,
                stress_Ry_Bohr3=np.diag([.001,.001,.002]).tolist(),
                strains=[np.eye(3).tolist()], explicit_Ry=[3.0],
                induced_Ry=[-2.5], derivatives_Ry=[.5],
                fixed_energy_explicit_parts=dict(lattice_geometry_Ry=[1.0],
                    core_explicit_Ry=[-5.0],density_explicit_Ry=[6.0]),
            )
            reportpath = work/'report.json'
            reportpath.write_text(json.dumps(report))
            status = dict(run_id='test',reference_id='reference',status='complete',
                requested_mode='STRFT',phase='final',backend='analytic',
                stress_ready=True,converged=True,stress_applied=False,
                workspace=str(work),report=str(reportpath),parent_fingerprints={},
                report_sha256=hashlib.sha256(reportpath.read_bytes()).hexdigest())
            Path('sample.stress_response.status.json').write_text(json.dumps(status))
            for suffix, label in [('0','DENS'),('2','VAL'),('c','CORE')]:
                Path('sample.scf'+suffix).write_text('\n'.join(l for l in tensors().splitlines()
                                                            if ':STR_'+label in l)+'\n')
            yield root,status,report
        finally:
            os.chdir(old)


class SummaryTests(unittest.TestCase):
    def test_both_mergers_publish_same_tensor_and_report_correct_physical_signs(self):
        reports=[]
        for standalone in (True,False):
            with self.subTest(standalone=standalone),case() as (root,status,report):
                files=[Path('sample.scf'+s) for s in ['0','2','c']]
                saved=[p.read_bytes() for p in files]
                summary.begin('sample',status)
                summary.finish('sample',status)
                if standalone:
                    merge_scf.merge('sample',files,replace=True)
                else:
                    scf_merge.scf('sample',dict(updn=False,m_extn=''))
                np.testing.assert_allclose(-mixer_accumulator(Path('sample.scf').read_text()),
                                           report['stress_Ry_Bohr3'],atol=1e-14)
                self.assertEqual(saved,[p.read_bytes() for p in files])
                text=Path('sample.outputdmf2_response').read_text()
                rows={line[:29].strip():np.array(list(map(float,line[29:].split())))
                      for line in text.splitlines() if line.startswith(
                          ('Original total','Net response correction','FINAL TOTAL'))}
                expected=np.array([.1,.2,.3,.2])*summary.RY_BOHR3_TO_GPA
                np.testing.assert_allclose(rows['Original total'],expected,atol=1e-6)
                np.testing.assert_allclose(rows['Original total']+rows['Net response correction'],
                                           rows['FINAL TOTAL'],atol=1.1e-6)
                self.assertIn('SCF publication complete',text)
                self.assertEqual((Path(status['workspace'])/'outputdmf2_response.txt').read_text(),text)
                reports.append(rows)
        for key in reports[0]:np.testing.assert_array_equal(reports[0][key],reports[1][key])

    def test_rebuild_uses_last_matching_iteration_not_whole_history(self):
        with case() as (root,status,report):
            files=[Path('sample.scf'+s) for s in ['0','2','c']]
            merge_scf.merge('sample',files,replace=True)
            current=Path('sample.scf').read_text()
            Path('sample.scf').write_text(current.replace('reference','old-reference')+current)
            before=Path('sample.scf').read_bytes()
            Path('sample.outputdmf2_response').unlink()
            path=summary.rebuild('sample')
            text=path.read_text()
            self.assertIn('Reconstructed from saved results',text)
            original=next(l for l in text.splitlines() if l.startswith('Original total'))
            self.assertAlmostEqual(float(original.split()[-1]),.2*summary.RY_BOHR3_TO_GPA,places=6)
            self.assertEqual(before,Path('sample.scf').read_bytes())

    def test_saved_report_tampering_is_rejected_before_overwriting_summary(self):
        with case() as (root,status,report):
            path=Path('sample.outputdmf2_response');path.write_text('keep me')
            Path(status['report']).write_text('{}')
            with self.assertRaisesRegex(ValueError,'hash mismatch'):summary.rebuild('sample')
            self.assertEqual(path.read_text(),'keep me')

    def test_searched_energy_completion_and_induced_term_are_separate(self):
        with case() as (root,status,report):
            text=summary._report_text(report)
            line=next(l for l in text.splitlines() if l.startswith('Radial-energy-search completion'))
            self.assertEqual(float(line.split()[-1]),1.0)
            self.assertIn('NOT the old-stress/net-correction split',text)
            self.assertIn('Compression-positive pressure: -',text)

    def test_window_omission_warning_is_in_readable_summary(self):
        with case() as (root, status, report):
            report['window_gap_treatment'] = [dict(spin_channel='',
                omitted_pairs=9, affected_kpoints=8,
                affected_kpoint_weight=.029513888888888888, cutoff_Ry=1e-10)]
            text = summary._report_text(report)
            self.assertIn('WARNING [primary]: Window cutoff', text)
            self.assertIn('9 singular cross-window band pairs on 8 k points', text)
            self.assertIn('1e-10 Ry', text)
            self.assertIn('0.029513889', text)

    def test_charge_warning_distinguishes_solver_and_central_convergence(self):
        with case() as (root, status, report):
            report['charge_fixed_point'] = dict(residual_norm=.00503141,
                warning_threshold=.001, within_threshold=False,
                warning='Central charge-map residual 0.00503141 exceeds threshold 0.001.')
            text = summary._report_text(report)
            self.assertIn('WARNING: Central charge-map residual 0.00503141', text)
            self.assertIn('Response solver converged: True', text)
            self.assertIn('0.001 / False', text)

    def test_shear_values_are_printed_when_present(self):
        tensor=np.array([[1.,2.,3.],[2.,4.,5.],[3.,5.,6.]])
        text=summary._stress_table([('Complete stress',tensor)])
        self.assertIn('Shear components',text)
        row=text.split('Shear components')[1].splitlines()[1]
        self.assertEqual(list(map(float,row[29:].split())),[5.,3.,2.])

    def test_reporting_failure_does_not_invalidate_successful_publication(self):
        with case() as (root,status,report):
            correction=response_rows('sample',tensors())
            with patch.object(summary,'_append',side_effect=OSError('disk full')), \
                 patch('sys.stderr',new_callable=io.StringIO) as warning:
                mark_published('sample',correction,'sample.scf',current_iteration=tensors())
            self.assertTrue(json.loads(Path('sample.stress_response.status.json').read_text())['stress_applied'])
            self.assertIn('disk full',warning.getvalue())


if __name__ == '__main__':
    unittest.main()
