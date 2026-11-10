from pathlib import Path
import sys,re,importlib,hashlib,json,unittest
root=Path('/private/tmp/edmft_shared_gamma_response_20261108')
repo=Path('/Users/haule/projects/eDMFT_new')
source=Path(sys.argv[1]);sys.path[:0]=[str(source),str(repo/'src/python/tests')]
for name in ['lapw_response_chain','frozen_sigma_tangent']:
    mod=importlib.import_module(name)
    assert Path(mod.__file__).parent==source,(name,mod.__file__)
    print('Installed import:',mod.__file__,flush=True)
modules=sorted(set(re.findall(r'\((test_\w+)\.',(root/'full_tests.log').read_text()))|{'test_lapw_core_response','test_lapw_xc_response'})
result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromNames(modules))
if not result.wasSuccessful() or result.skipped or result.testsRun!=169:raise SystemExit(1)
for item in json.loads((repo/'src/doc/stress/shared_gamma_response_20261108/installation_manifest.json').read_text()):
    assert hashlib.sha256(Path(item['destination']).read_bytes()).hexdigest()==item['sha256'],item['destination']
print('169 tests passed, zero skips; all installed hashes verified',flush=True)
