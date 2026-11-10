from pathlib import Path
root=Path('/private/tmp/edmft_shear_magnetic_20261109')
p=root/'stage/lapw_response_policy.py'
s=p.read_text().replace("    fixed = override.get('fixed_number', bool(physical_mode))", "    if 'fixed_number' not in override and (type(physical_mode) is not int or physical_mode not in (0, 1)):\n        raise ValueError('physical recomputeEF must be 0 or 1 for response')\n    fixed = override.get('fixed_number', bool(physical_mode))")
s=s.replace("    params = {} if", "    if any(not isinstance(value, bool) for value in override.values()):\n        raise ValueError('physical policy flags must be booleans')\n    params = {} if")
p.write_text(s)
p=root/'stage/tests/test_dmft2_spin.py'
s=p.read_text().replace('    def test_parameter_read_does_not_execute_file(self):', '''    def test_ferro_rejects_ambiguous_physical_ensemble(self):
        with fixture(ferro=True):
            for value in [2, '1', -1]:
                Path('params.dat').write_text(f'Ferro=True\\nrecomputeEF={value!r}\\n')
                with self.assertRaisesRegex(ValueError, 'physical recomputeEF'):
                    detect_policy('.', 'sample')
            with self.assertRaisesRegex(ValueError, 'booleans'):
                detect_policy('.', 'sample', override={'ferro': True, 'fixed_number': 1})

    def test_parameter_read_does_not_execute_file(self):''')
p.write_text(s)
# Controlled asymmetric-spin native smoke test. Modify only this isolated
# reference, never the user's saved MnTe calculations or the running snapshot.
import json,shutil,numpy as np
cfg=json.loads((root/'prepare_MnTe.json').read_text())
ref=root/'ferro_reference';ref.mkdir(exist_ok=False)
work=root/'native_MnTe'
for p in work.iterdir():
    if (p.name.startswith('0_800_MnTe.') and p.suffix in {'.struct','.in0','.in1','.in1c','.in2','.inc','.inm','.klist','.kgen','.rotlm','.indmfl','.indmfldn','.indmfi','.indmf2','.indmf2dn','.clmsum'}) or p.name in {'EF.dat','projectorw.dat','Eorb_imp.dat','params.dat','lapw1.def','dmft2.def','dmft2dn.def','sig.inp','sig.inp1','sig.inp2','sig.inp1dn','sig.inp2dn'}:
        shutil.copy2(p,ref/p.name)
for p in [ref/'sig.inp1dn',ref/'sig.inp2dn']:
    data=np.loadtxt(p);data[:,1::2]+=.2;np.savetxt(p,data,fmt='%.18e')
with (ref/'params.dat').open('a') as f:f.write('\n# Isolated asymmetric-spin sequencing test only.\nFerro=True\nrecomputeEF=1\n')
cfg.update(reference=str(ref),workspace=str(root/'native_Ferro'),lapw1_definition=str(ref/'lapw1.def'),dmft2_definition=str(ref/'dmft2.def'),dmft2dn_definition=str(ref/'dmft2dn.def'))
(root/'prepare_Ferro.json').write_text(json.dumps(cfg,indent=2)+'\n')
