from pathlib import Path
import sys,re,importlib,hashlib,json,unittest
root=Path('/private/tmp/edmft_charge_closure_20261108')
repo=Path('/Users/haule/projects/eDMFT_new')
source=Path(sys.argv[1]);sys.path[:0]=[str(source),str(repo/'src/python/tests')]
for name in ['lapw_charge_response','lapw_potential_core_state','dmft2_step','wien_response_operators','wien_frozen_sigma_response','lapw_hartree_periodic']:
    mod=importlib.import_module(name)
    assert Path(mod.__file__).parent==source,(name,mod.__file__)
    print('Installed import:',mod.__file__,flush=True)
modules=sorted(set(re.findall(r'\((test_\w+)\.',(root/'full_tests.log').read_text()))|{'test_lapw_core_response','test_lapw_xc_response'})
result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromNames(modules))
if not result.wasSuccessful() or result.skipped or result.testsRun!=161:raise SystemExit(1)
for item in json.loads((repo/'src/doc/stress/charge_closure_20261108/installation_manifest.json').read_text()):
    assert hashlib.sha256(Path(item['destination']).read_bytes()).hexdigest()==item['sha256'],item['destination']
print('161 tests passed, zero skips; all installed hashes verified',flush=True)
