"""Post-installation tests with installed imports only; no staged modules."""
from pathlib import Path
import sys,re,unittest,hashlib,importlib,json
root=Path(__file__).resolve().parent
repo=Path('/Users/haule/projects/eDMFT_new')
module_dir=repo/'src/python' if sys.argv[1]=='repository' else Path('/Users/haule/dbin')
sys.path[:0]=[str(module_dir),str(repo/'src/python/tests')]
record=json.loads((repo/'src/doc/stress/xc_core_response_20261108/INSTALLATION.json').read_text())
for row in record['files']:
 p=Path(row['target']);assert hashlib.sha256(p.read_bytes()).hexdigest()==row['new_sha256'],str(p)
for name in ['lapw_hartree_periodic','lapw_xc_response','lapw_core_response']:
 m=importlib.import_module(name)
 assert Path(m.__file__).parent==module_dir,(name,m.__file__)
 print('Installed import:',m.__file__,flush=True)
modules=sorted(set(re.findall(r'\((test_\w+)\.',(root/'full_tests.log').read_text())))
for name in ['test_lapw_core_response','test_lapw_xc_response']:
 if name not in modules:modules.append(name)
suite=unittest.defaultTestLoader.loadTestsFromNames(modules)
result=unittest.TextTestRunner(verbosity=2).run(suite)
print('Verified all installed hashes and',result.testsRun,'tests; mode',sys.argv[1],flush=True)
raise SystemExit(not result.wasSuccessful() or bool(result.skipped))
