"""Publication, retention, opt-out and scope tests using small local fixtures."""
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest
import response_cleanup as cleanup


def js(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,sort_keys=True)+'\n')


def fixture(root,case='FeSe'):
    work=root/(case+'.response_runs')/'current-run';central=work/'central';central.mkdir(parents=True)
    report=dict(case=case,converged=True,force_mRy_Bohr=[[.1,.2,.3],[-.1,-.2,-.3]],adjoint=dict(converged=True,full_frequency_check=True,residual_norm=1e-10,threshold=1e-8))
    for folder in ['total_force','total_stress']:js(work/folder/'report.json',report)
    (root/(case+'.scf')).write_text(':FORRESP_BEGIN: ref-current complete response\n:FORRESP_END: response\n')
    (root/(case+'.vector')).write_bytes(b'PARENT VECTOR')
    (root/(case+'.clmsum')).write_bytes(b'OLD density')
    status=dict(status='complete',converged=True,backend='analytic',requested_mode='FORT',run_id=work.name,workspace=str(work),reference_id='ref-current',force_applied=True,stress_applied=False,published_scf=str(root/(case+'.scf')),keep_temporary_files=False,parent_fingerprints={case+'.clmsum':cleanup._sha(root/(case+'.clmsum'))})
    for field,folder in [('report','total_stress'),('force_report','total_force')]:
        path=work/folder/'report.json';status[field]=str(path);status[field+'_sha256']=cleanup._sha(path)
    js(root/(case+'.stress_response.status.json'),status)
    (root/(case+'.scfm')).write_text(':DIS  1.e-7\n:FGL001: 1.ATOM .1 .2 .3 total forces\n')
    generated=['native.k00000001.bin','native.ham.k00000001.bin','response.k00000001.bin','response_dn.k00000001.bin','response.operator.k00000001.bin','response.spin',case+'.vector',case+'.vectorso',case+'.vectorsodn_1.gz']
    for name in generated:(central/name).write_bytes(('regenerate '+name).encode())
    for name in ['response.meta','soc.geometry','potential_core.npz','unknown.bin',case+'.vector_notes.txt','dmft2.def','reference.json']:(central/name).write_bytes(('keep '+name).encode())
    return work,central,status


class CleanupTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name).resolve()/'case';self.root.mkdir();self.work,self.central,self.status=fixture(self.root)
    def tearDown(self):self.tmp.cleanup()
    def call(self,status=None):return cleanup.cleanup_response('FeSe',directory=self.root,status=status)
    def save_status(self):
        js(self.root/'FeSe.stress_response.status.json',self.status)
        # Fixtures using this helper represent a new successful publication
        # followed by mixer, unless a test explicitly makes scfm stale/missing.
        (self.root/'FeSe.scfm').write_text(':DIS 1.e-7\n:FGL001: 1.ATOM .1 .2 .3 total forces\n')
    def rewrite_reports(self,change):
        for field in ['report','force_report']:
            p=Path(self.status[field]);d=json.loads(p.read_text());change(d);js(p,d);self.status[field+'_sha256']=cleanup._sha(p)
        self.save_status()

    def test_success_generic_case_preserves_parent_inputs_and_metadata(self):
        kept=[self.root/'FeSe.vector',self.root/'FeSe.clmsum',self.central/'response.meta',self.central/'soc.geometry',self.central/'potential_core.npz',self.central/'unknown.bin',self.central/'FeSe.vector_notes.txt',self.central/'dmft2.def',self.central/'reference.json']
        before={p:p.read_bytes() for p in kept};result=self.call()
        self.assertEqual(result['status'],'cleaned');self.assertEqual(result['removed_files'],9)
        for p,data in before.items():self.assertEqual(p.read_bytes(),data)
        manifest=json.loads((self.work/cleanup.MANIFEST).read_text());self.assertEqual(len(manifest['files']),9)
        self.assertTrue(all(len(x['sha256'])==64 and x['bytes']>0 for x in manifest['files']))
        self.assertTrue(all(not (self.work/x['path']).exists() for x in manifest['files']))

    def test_no_status_or_unapplied_native_returns_skip(self):
        p=self.root/'FeSe.stress_response.status.json';p.unlink();self.assertEqual(self.call()['reason'],'no_response_status')
        self.status['force_applied']=False;self.save_status();self.assertEqual(self.call()['reason'],'response_not_applied')
        self.assertTrue((self.central/'FeSe.vector').exists())

    def test_failed_applied_status_is_refused_without_deletion(self):
        self.status['status']='failed';self.save_status()
        with self.assertRaises(cleanup.CleanupError):self.call()
        self.assertTrue((self.central/'FeSe.vector').exists())

    def test_keep_option_preserves_everything(self):
        self.status['keep_temporary_files']=True;self.save_status();self.assertEqual(self.call()['reason'],'keep_temporary_files')
        self.assertTrue((self.central/'response.spin').exists());self.assertFalse((self.work/cleanup.MANIFEST).exists())

    def test_keep_option_requires_real_boolean(self):
        self.status['keep_temporary_files']='false';self.save_status()
        with self.assertRaises(cleanup.CleanupError):self.call()

    def test_parent_density_may_change_after_publication(self):
        (self.root/'FeSe.clmsum').write_bytes(b'NEW density after successful mixer')
        self.assertEqual(self.call()['status'],'cleaned')
        self.assertEqual((self.root/'FeSe.clmsum').read_bytes(),b'NEW density after successful mixer')

    def test_report_hash_mismatch_refuses_before_deletion(self):
        Path(self.status['force_report']).write_text('{}')
        with self.assertRaises(cleanup.CleanupError):self.call()
        self.assertTrue((self.central/'native.k00000001.bin').exists())

    def test_missing_adjoint_is_refused_for_screened_FORT(self):
        self.rewrite_reports(lambda d:d.pop('adjoint'))
        with self.assertRaisesRegex(cleanup.CleanupError,'missing its adjoint'):self.call()
        self.assertTrue((self.central/'FeSe.vector').exists())

    def test_full_frequency_failure_refuses_with_valid_hashes(self):
        self.rewrite_reports(lambda d:d['adjoint'].update(full_frequency_check=False))
        with self.assertRaisesRegex(cleanup.CleanupError,'full-frequency'):self.call()
        self.assertTrue((self.central/'FeSe.vector').exists())

    def test_FOR_explicit_operator_can_legitimately_lack_adjoint(self):
        self.status['requested_mode']='FOR';self.save_status();self.rewrite_reports(lambda d:d.pop('adjoint'))
        self.assertEqual(self.call()['status'],'cleaned')

    def test_publication_marker_or_case_mismatch_refuses(self):
        (self.root/'FeSe.scf').write_text(':FORRESP_BEGIN: old-ref result\n')
        with self.assertRaises(cleanup.CleanupError):self.call()
        (self.root/'FeSe.scf').write_text(':FORRESP_BEGIN: ref-current result\n')
        self.rewrite_reports(lambda d:d.update(case='Ce'))
        with self.assertRaises(cleanup.CleanupError):self.call()

    def test_older_run_and_EOS_archives_remain_untouched(self):
        paths=[self.root/'FeSe.response_runs/previous-failed/central/FeSe.vector',self.root/'EOS_archive/current-run/FeSe.vector',self.work/'eos_archive/FeSe.vector']
        for p in paths:p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'archive evidence')
        self.call();self.assertTrue(all(p.read_bytes()==b'archive evidence' for p in paths))

    def test_external_symlink_and_hardlink_are_preserved(self):
        outside=self.root.parent/'outside';outside.mkdir();original=outside/'source';original.write_bytes(b'external data')
        (self.central/'native.k00000099.bin').symlink_to(original)
        (self.central/'external_directory').symlink_to(outside,target_is_directory=True)
        os.link(original,self.central/'native.k00000098.bin')
        self.call();self.assertEqual(original.read_bytes(),b'external data')
        self.assertTrue((self.central/'native.k00000099.bin').is_symlink());self.assertTrue((self.central/'external_directory').is_symlink())
        self.assertTrue((self.central/'native.k00000098.bin').exists())

    def test_another_run_id_or_external_report_is_refused(self):
        self.status['run_id']='previous-failed';self.save_status()
        with self.assertRaises(cleanup.CleanupError):self.call()
        self.status['run_id']=self.work.name
        outside=self.root.parent/'report.json';js(outside,dict(converged=True))
        self.status['force_report']=str(outside);self.status['force_report_sha256']=cleanup._sha(outside);self.save_status()
        with self.assertRaises(cleanup.CleanupError):self.call()

    def test_archived_case_status_cannot_cleanup_original_calculation(self):
        archive=self.root.parent/'EOS_archive';archive.mkdir()
        js(archive/'FeSe.stress_response.status.json',self.status);(archive/'FeSe.scf').write_text((self.root/'FeSe.scf').read_text())
        with self.assertRaises(cleanup.CleanupError):cleanup.cleanup_response('FeSe',directory=archive)
        self.assertTrue((self.central/'FeSe.vector').exists())

    def test_supplemental_status_must_match_persisted_publication(self):
        supplied=dict(self.status,reference_id='different')
        with self.assertRaises(cleanup.CleanupError):self.call(supplied)

    def test_idempotence_preserves_first_hash_inventory(self):
        self.call();p=self.work/cleanup.MANIFEST;old=p.read_bytes()
        self.assertEqual(self.call()['status'],'already_cleaned');self.assertEqual(p.read_bytes(),old)

    def test_stress_only_publication_is_supported(self):
        self.status.update(force_applied=False,stress_applied=True,requested_mode='STRT');self.save_status()
        (self.root/'FeSe.scf').write_text(':STRRESP_BEGIN: ref-current complete response\n')
        self.assertEqual(self.call()['status'],'cleaned')

    def test_nonempty_mixer_error_preserves_every_regenerable_file(self):
        (self.root/'mixer.error').write_text('mixer failed')
        with self.assertRaisesRegex(cleanup.CleanupError,'error file'):self.call()
        self.assertTrue((self.central/'FeSe.vector').exists())
        self.assertFalse((self.work/cleanup.MANIFEST).exists())

    def test_stale_or_missing_scfm_preserves_every_regenerable_file(self):
        p=self.root/'FeSe.scfm';publication=(self.root/'FeSe.stress_response.status.json').stat().st_mtime_ns
        os.utime(p,ns=(publication-1000000000,publication-1000000000))
        with self.assertRaisesRegex(cleanup.CleanupError,'predates'):self.call()
        self.assertTrue((self.central/'FeSe.vector').exists())
        p.unlink()
        with self.assertRaises(cleanup.CleanupError):self.call()
        self.assertTrue((self.central/'FeSe.vector').exists())

    def test_fresh_but_empty_mixer_output_is_not_success_evidence(self):
        (self.root/'FeSe.scfm').write_text('old header only\n')
        with self.assertRaisesRegex(cleanup.CleanupError,'no charge/force/stress'):self.call()
        self.assertTrue((self.central/'FeSe.vector').exists())

    def test_empty_native_error_files_do_not_prevent_cleanup(self):
        (self.root/'mixer.error').write_bytes(b'')
        self.assertEqual(self.call()['status'],'cleaned')

if __name__=='__main__':unittest.main(verbosity=2)
