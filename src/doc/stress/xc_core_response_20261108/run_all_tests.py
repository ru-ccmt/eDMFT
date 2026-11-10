from pathlib import Path
import sys,re,unittest
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),str(root/'stage/tests'),'/Users/haule/projects/eDMFT_new/src/python','/Users/haule/projects/eDMFT_new/src/python/tests']
text=Path('/private/tmp/edmft_hartree_periodic_20261108/full_tests.log').read_text()
modules=sorted(set(re.findall(r'\((test_\w+)\.',text)))
modules+=['test_lapw_core_response','test_lapw_xc_response']
print('Modules:',', '.join(modules),flush=True)
suite=unittest.defaultTestLoader.loadTestsFromNames(modules)
result=unittest.TextTestRunner(verbosity=2).run(suite)
raise SystemExit(not result.wasSuccessful())
