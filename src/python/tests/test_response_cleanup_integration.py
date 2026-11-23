import ast
from contextlib import contextmanager
import io,json,os,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import dmft2_step as step

class IntegrationTests(unittest.TestCase):
    def function(self, events, fail=None):
        source=Path(step.__file__).with_name('run_dmft.py')
        tree=ast.parse(source.read_text())
        node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='final_response_pass')
        scope={'dmft2_step':SimpleNamespace(request=lambda case:dict(enabled=True,mode='FORT'))}
        for name in ['dmft2','lcore','scf','scf1','mixer','scfm']:
            def call(*args,n=name,**kwargs):
                events.append(n)
                if fail==n:raise RuntimeError('failure in '+n)
            scope[name]=call
        exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),scope)
        return scope['final_response_pass']
    def args(self):
        return [io.StringIO(),SimpleNamespace(case='sample',WIENROOT='.'),io.StringIO(),None,{'force_last':False},'',{'m_extn':''},True]
    def test_cleanup_occurs_only_after_successful_mixer_and_scfm(self):
        events=[];args=self.args()
        def clean(case):events.append('cleanup');return dict(removed_files=2,removed_bytes=99)
        def parent(case):events.append('parent_cleanup');return dict(removed_files=1,removed_bytes=10)
        with patch.dict(sys.modules,{'response_cleanup':SimpleNamespace(cleanup_response=clean,cleanup_parent_vectors=parent)}):
            self.assertTrue(self.function(events)(*args))
        self.assertEqual(events,['dmft2','lcore','scf','scf1','mixer','scfm','cleanup','parent_cleanup'])
        self.assertIn('109 bytes',args[2].getvalue())
    def test_response_or_mixer_failure_retains_temporaries(self):
        for fail in ['dmft2','lcore','mixer','scfm']:
            events=[]
            with patch.dict(sys.modules,{'response_cleanup':SimpleNamespace(cleanup_response=lambda case:events.append('cleanup'),cleanup_parent_vectors=lambda case:events.append('parent_cleanup'))}):
                with self.assertRaises(RuntimeError):self.function(events,fail)(*self.args())
            self.assertNotIn('cleanup',events);self.assertNotIn('parent_cleanup',events)
    def test_cleanup_error_is_visible_and_does_not_invalidate_forces(self):
        events=[];args=self.args()
        def clean(case):raise OSError('storage error')
        with patch.dict(sys.modules,{'response_cleanup':SimpleNamespace(cleanup_response=clean,cleanup_parent_vectors=lambda case:dict(removed_files=0,removed_bytes=0))}):
            self.assertTrue(self.function(events)(*args))
        self.assertIn('WARNING',args[2].getvalue());self.assertIn('storage error',args[2].getvalue())
    def test_retention_setting_is_validated_and_captured_in_status(self):
        old=Path.cwd()
        with tempfile.TemporaryDirectory() as tmp:
            os.chdir(tmp)
            try:
                Path('dmft2.def').write_text("2, 'sample.indmf2', 'old', 'formatted', 0\n5, 'sample.in2', 'old', 'formatted', 0\n")
                Path('sample.indmf2').write_text('c 1.0 1\n');Path('sample.in2').write_text('FORT\n')
                cfg=Path('sample.stress_response.json')
                for value in [True,False]:
                    cfg.write_text(json.dumps({'keep_temporary_files':value}))
                    self.assertIs(step.request('sample')['config']['keep_temporary_files'],value)
                    with patch.object(step,'run_core'),patch.object(step,'preflight',return_value={}):
                        result=step.execute('sample','/unused',phase='defer')
                    self.assertIs(result['keep_temporary_files'],value)
                for value in ['false',0,None]:
                    cfg.write_text(json.dumps({'keep_temporary_files':value}))
                    with self.assertRaisesRegex(ValueError,'keep_temporary_files'):step.request('sample')
            finally:os.chdir(old)
if __name__=='__main__':unittest.main(verbosity=2)
