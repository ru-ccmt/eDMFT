"""Magnetic detection, common-mu sequencing and preservation of live inputs."""
from pathlib import Path
from contextlib import contextmanager
import io,os,tempfile,unittest
from unittest.mock import patch
from lapw_response_policy import detect_policy,control_fields
from dmft2_spin import run_channels
import dmft2_step as step


@contextmanager
def fixture(ferro=False, fixed=True, magnetic=True):
    old=os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        os.chdir(tmp)
        try:
            Path('sample.inm').write_text('MSR1 0.0 '+('YES' if fixed else 'NO')+'\n')
            Path('params.dat').write_text(f'Ferro={ferro!r}\nrecomputeEF={int(fixed)}\n')
            Path('sample.in2').write_text('STRFT\n')
            Path('EF.dat').write_text('1.25\n')
            for suffix in ['', 'dn'] if magnetic else ['']:
                Path('sample.indmfl'+suffix).write_text('1 4 1 5\n')
                Path('sample.indmf2'+suffix).write_text(f'# keep header\nc 1.25 {int(fixed and not ferro and not suffix)} 0.7 .T. 0.02\nkeep trailing controls\n')
                Path('dmft2'+suffix+'.def').write_text(''.join(f"{u}, '{name}', 'old', 'formatted', 0\n" for u,name in [(2,'sample.indmf2'+suffix),(5,'sample.in2'),(7,'sample.indmfl'+suffix),(15,'sample.bnds'+suffix)]))
            yield
        finally:os.chdir(old)


class SpinWorkflowTests(unittest.TestCase):
    def test_detects_ferro_fixed_N_despite_temporary_zero(self):
        with fixture(ferro=True):
            p=detect_policy('.', 'sample')
            self.assertTrue(p.magnetic and p.ferro and p.fixed_number)
            self.assertEqual(p.channels,('', 'dn'))
        with fixture(fixed=False):
            p=detect_policy('.', 'sample')
            self.assertFalse(p.fixed_number or p.normalization_enabled)
        with fixture(magnetic=False):self.assertFalse(detect_policy('.','sample').magnetic)

    def test_ordinary_pair_shares_first_mu_and_preserves_inputs(self):
        with fixture():
            old={p:p.read_bytes() for p in Path('.').glob('*') if p.name!='EF.dat'}
            seen=[]
            def run(label,df):
                _,_,words=control_fields(step.entries(Path(df).read_text())[2]);seen.append((label,words[:3]))
                if label=='dmft2':Path('EF.dat').write_text('1.47\n')
            run_channels(['dmft2.def','dmft2dn.def'],detect_policy('.','sample'),run,lambda _:self.fail('unexpected fermif'))
            self.assertEqual(seen,[('dmft2',['c','1.25','1']),('dmft2dn',['c','1.47','0'])])
            for p,b in old.items():self.assertEqual(p.read_bytes(),b)
            self.assertFalse(list(Path('.').glob('.dmft2-*-input-*')))

    def test_ferro_spectra_then_common_mu_then_both_charges(self):
        with fixture(ferro=True):
            seen=[]
            def run(label,df):
                _,_,w=control_fields(step.entries(Path(df).read_text())[2]);seen.append((label,w[0],int(w[2]),float(w[1])))
            def fermi(bands):
                self.assertEqual(bands,['sample.bnds','sample.bndsdn']);seen.append('fermif');Path('EF.dat').write_text('1.61\n')
            run_channels(['dmft2.def','dmft2dn.def'],detect_policy('.','sample'),run,fermi)
            self.assertEqual([x[0] if isinstance(x,tuple) else x for x in seen],['dmft2_eigen','dmft2_eigendn','fermif','dmft2','dmft2dn'])
            self.assertEqual(seen[-2][1:],('c',0,1.61));self.assertEqual(seen[-1][1:],('c',0,1.61))

    def test_export_preserves_physical_ensemble_but_never_refits_mu(self):
        with fixture(ferro=True):
            p=detect_policy('.','sample');seen=[]
            def run(label,df):seen.append(control_fields(step.entries(Path(df).read_text())[2])[2][:3])
            run_channels(['dmft2.def','dmft2dn.def'],p,run,lambda _:self.fail('export refit'),export=True)
            self.assertTrue(p.fixed_number);self.assertEqual(seen,[['c','1.25','0']]*2)

    def test_executor_averages_once_before_shared_response(self):
        with fixture():
            events=[]
            with patch.object(step,'preflight'),patch.object(step,'run_core',side_effect=lambda *a:events.append('core')),patch('dmft2_spin.average_charge',side_effect=lambda *a:events.append('average')),patch.object(step,'compute_response',side_effect=lambda *a:events.append('response') or {}):
                step.execute('sample','/unused',phase='final',output=io.StringIO())
            self.assertEqual(events,['core','core','average','response'])

    def test_no_extra_charge_solve_during_native_search_phase(self):
        for ferro,expected in [(False,['dmft2']),(True,['dmft2_eigen','dmft2_eigendn','fermif'])]:
            with fixture(ferro=ferro):
                seen=[]
                run_channels(['dmft2.def','dmft2dn.def'],detect_policy('.','sample'),lambda label,df:seen.append(label),lambda _:seen.append('fermif'),search_only=True)
                self.assertEqual(seen,expected)

    def test_ferro_rejects_ambiguous_physical_ensemble(self):
        with fixture(ferro=True):
            for value in [2, '1', -1]:
                Path('params.dat').write_text(f'Ferro=True\nrecomputeEF={value!r}\n')
                with self.assertRaisesRegex(ValueError, 'physical recomputeEF'):
                    detect_policy('.', 'sample')
            with self.assertRaisesRegex(ValueError, 'booleans'):
                detect_policy('.', 'sample', override={'ferro': True, 'fixed_number': 1})

    def test_parameter_read_does_not_execute_file(self):
        with fixture():
            with Path('params.dat').open('a') as f:f.write("raise RuntimeError('must not execute')\n")
            self.assertTrue(detect_policy('.','sample').fixed_number)
