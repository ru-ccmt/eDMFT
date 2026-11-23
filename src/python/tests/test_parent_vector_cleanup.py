import json,os,unittest
from pathlib import Path
import response_cleanup as cleanup
import test_response_cleanup as fixtures

class ParentVectorTests(fixtures.CleanupTests):
    # This class reuses the success-publication fixture, not the workspace tests.
    def call_parent(self):return cleanup.cleanup_parent_vectors('FeSe',directory=self.root)
    def test_parent_scalar_SOC_split_and_compressed_vectors_removed(self):
        names=['FeSe.vector','FeSe.vector_12','FeSe.vectorso','FeSe.vectorso_16','FeSe.vectorsodn_63','FeSe.vectorsodn_2.gz','FeSe.vectorup_3','FeSe.vectordn','FeSe.vectorsoup_7']
        for n in names:(self.root/n).write_bytes(n.encode())
        keep=['FeSe.energyso','FeSe.clmsum','sig.inp','FeSe.struct','FeSe.vectorso_notes.txt','Other.vectorso_1','fort_inputs.tgz']
        for n in keep:(self.root/n).write_bytes(('keep '+n).encode())
        before={n:(self.root/n).read_bytes() for n in keep}
        out=self.call_parent();self.assertEqual(out['removed_files'],len(names))
        self.assertTrue(all(not (self.root/n).exists() for n in names))
        for n,b in before.items():self.assertEqual((self.root/n).read_bytes(),b)
        self.assertTrue((self.central/'FeSe.vector').exists())
        m=json.loads((self.work/'parent_vector_cleanup.json').read_text());self.assertEqual({x['path'] for x in m['files']},set(names))
        self.assertTrue(all(len(x['sha256'])==64 for x in m['files']))
    def test_parent_keep_option_preserves_vectors(self):
        self.status['keep_temporary_files']=True;self.save_status()
        self.assertEqual(self.call_parent()['reason'],'keep_temporary_files');self.assertTrue((self.root/'FeSe.vector').exists())
    def test_parent_cleanup_requires_successful_current_publication(self):
        for kind in ['failed','mixer_error','bad_report']:
            with self.subTest(kind=kind):
                if kind=='failed':self.status['status']='failed';self.save_status()
                elif kind=='mixer_error':self.status['status']='complete';self.save_status();(self.root/'mixer.error').write_text('failed')
                else:(self.root/'mixer.error').write_text('');Path(self.status['force_report']).write_text('{}')
                with self.assertRaises(cleanup.CleanupError):self.call_parent()
                self.assertTrue((self.root/'FeSe.vector').exists())
    def test_parent_symlinks_hardlinks_and_other_directories_preserved(self):
        external=self.root.parent/'external';external.write_bytes(b'keep')
        (self.root/'FeSe.vectorso_1').symlink_to(external)
        os.link(external,self.root/'FeSe.vectorso_2')
        sub=self.root/'unrelated';sub.mkdir();(sub/'FeSe.vectorso_3').write_bytes(b'keep nested')
        self.call_parent();self.assertEqual(external.read_bytes(),b'keep');self.assertTrue((self.root/'FeSe.vectorso_1').is_symlink());self.assertTrue((self.root/'FeSe.vectorso_2').exists());self.assertTrue((sub/'FeSe.vectorso_3').exists())
    def test_parent_preserved_alias_does_not_become_dangling(self):
        target=self.root/'FeSe.vectorso_0';target.write_bytes(b'keep target of alias')
        alias=self.root/'FeSe.vectorso';alias.symlink_to(target.name)
        self.call_parent();self.assertTrue(alias.exists());self.assertTrue(target.exists())
        self.assertEqual(alias.read_bytes(),b'keep target of alias')
    def test_parent_after_workspace_cleanup_and_idempotence(self):
        cleanup.cleanup_response('FeSe',directory=self.root)
        old=(self.work/'response_cleanup.json').read_bytes()
        self.call_parent();saved=(self.work/'parent_vector_cleanup.json').read_bytes()
        self.assertEqual(self.call_parent()['status'],'already_cleaned')
        self.assertEqual((self.work/'parent_vector_cleanup.json').read_bytes(),saved)
        self.assertEqual((self.work/'response_cleanup.json').read_bytes(),old)

# Do not inherit all workspace tests again in the reported test count.
for name in list(fixtures.CleanupTests.__dict__):
    if name.startswith('test_'):setattr(ParentVectorTests,name,None)
if __name__=='__main__':unittest.main(verbosity=2)
