import glob,io,json,os,unittest,textwrap
from pathlib import Path
from types import SimpleNamespace
import convert_processes
import dmft2_step as step
import test_response_cleanup as fixtures

class FooterTests(fixtures.CleanupTests):
    def run_footer(self,keep=False,mode='FOR'):
        root=self.root
        (root/'dmft2.def').write_text("2, 'FeSe.indmf2', 'old', 'formatted', 0\n5, 'FeSe.in2', 'old', 'formatted', 0\n")
        (root/'FeSe.indmf2').write_text('c 1 0\n');(root/'FeSe.in2').write_text(mode+'\n')
        (root/'FeSe.stress_response.json').write_text(json.dumps(dict(keep_temporary_files=keep)))
        code=Path(step.__file__).with_name('run_dmft.py').read_text()
        marker='    # Final cleanup must obey response retention and preserve linked vector targets.\n'
        tail=textwrap.dedent(code[code.index(marker):]);out=io.StringIO()
        scope=dict(os=os,glob=glob,dmft2_step=step,w2k=SimpleNamespace(case='FeSe'),istep0=1,block=None,fh_info=out)
        old=Path.cwd();os.chdir(root)
        try:exec(compile(tail,'actual_run_dmft_footer','exec'),scope)
        finally:os.chdir(old)
        return out.getvalue()
    def test_footer_preserves_split_vector_maps(self):
        maps = {'_processes_': 'native assignments', '_processes_0': 'rank assignments', '_processes_1': ''}
        for name, content in maps.items():
            (self.root/name).write_text(content)
        self.run_footer()
        for name, content in maps.items():
            self.assertEqual((self.root/name).read_text(), content)

    def test_wien_conversion_publishes_full_count_before_retiring_old_maps(self):
        previous = Path.cwd()
        os.chdir(self.root)
        try:
            Path('_processes_').write_text('old dftKS layout')
            Path('_processes_9').write_text('old WIEN2k rank')
            Path('.processes').write_text('1 : node : 2 : 1 : 1 : 0\n2 : node : 3 : 1 : 1 : 0\n')
            for i in (1,2):
                Path(f'lapw1_{i}.def').write_text(f"10, 'FeSe.vector_{i}',\n11, 'FeSe.energy_{i}',\n")
            convert_processes.convert(False, io.StringIO())
            rank_map = Path('_processes_0')
            saved = rank_map.read_text()
            self.assertEqual([line.split()[-2:] for line in saved.splitlines()], [['nk_total','5']]*2)
            self.assertFalse(Path('_processes_').exists())
            self.assertFalse(Path('_processes_9').exists())
            Path('lapw1_2.def').write_text('')
            with self.assertRaises(ValueError):
                convert_processes.convert(False, io.StringIO())
            self.assertEqual(rank_map.read_text(), saved)
        finally:
            os.chdir(previous)

    def test_footer_honors_keep_for_scalar_SOC_and_zero_files(self):
        for name in ['FeSe.vector_0','FeSe.vectorso_0','debug.empty']:(self.root/name).write_bytes(b'')
        self.run_footer(keep=True)
        for name in ['FeSe.vector','FeSe.vector_0','FeSe.vectorso_0','debug.empty']:self.assertTrue((self.root/name).exists())
        self.assertTrue((self.central/'response.spin').exists())
    def test_footer_cleans_ordinary_FOR_and_preserves_zero_alias_target(self):
        target=self.root/'FeSe.vector_0';target.write_bytes(b'')
        alias=self.root/'FeSe.vector_alias';alias.symlink_to(target.name)
        (self.root/'FeSe.vectorso_0').write_bytes(b'regenerate')
        (self.root/'dangling').symlink_to('missing')
        self.run_footer(mode='FOR')
        self.assertFalse((self.root/'FeSe.vector').exists());self.assertFalse((self.root/'FeSe.vectorso_0').exists())
        self.assertTrue(alias.is_symlink());self.assertTrue(alias.exists());self.assertTrue(target.exists())
        self.assertTrue((self.root/'dangling').is_symlink());self.assertTrue((self.root/'FeSe.clmsum').exists())
    def test_footer_retains_native_error_and_vectors_on_failed_mixer(self):
        error=self.root/'dmft2.error.0';error.write_text('failure')
        output=self.run_footer()
        self.assertTrue(error.exists());self.assertEqual(error.read_text(),'failure');self.assertTrue((self.root/'FeSe.vector').exists())
        self.assertIn('WARNING',output)
for name in list(fixtures.CleanupTests.__dict__):
    if name.startswith('test_'):setattr(FooterTests,name,None)
if __name__=='__main__':unittest.main(verbosity=2)
