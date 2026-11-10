from pathlib import Path
import hashlib,json,os,sys,unittest
r=Path(__file__).resolve().parent
record=Path('/Users/haule/projects/eDMFT_new/src/doc/stress/shear_magnetic_20261109/INSTALL.json')
for f in json.loads(record.read_text())['files']:
    assert hashlib.sha256(Path(f['path']).read_bytes()).hexdigest()==f['sha256'], f['path']
os.environ['EDMFT_RADIAL_LIBRARY']='/Users/haule/dbin/liblapw_response.so'
os.environ['EDMFT_HARTREE_REFERENCE_LIBRARY']='/private/tmp/edmft_hartree_radial_20261108/native/libhartree_reference.dylib'
os.environ['EDMFT_HARTREE_CHARGE_LIBRARY']='/private/tmp/edmft_hartree_periodic_20261108/native/libcharge_reference.dylib'
os.environ['EDMFT_MIXER_NORMALIZATION_LIBRARY']='/private/tmp/edmft_charge_closure_20261108/native/libmixer_normalization.dylib'
tests=Path('/Users/haule/projects/eDMFT_new/src/python/tests')
sys.path[:0]=['/Users/haule/dbin',str(tests)]
import dmft2_step,dmft2_spin,lapw_response_policy,lapw_total_stress,lapw_total_stress_driver
for m in [dmft2_step,dmft2_spin,lapw_response_policy,lapw_total_stress,lapw_total_stress_driver]:
    assert Path(m.__file__).parent==Path('/Users/haule/dbin'),m.__file__
suite=unittest.defaultTestLoader.discover(str(tests))
result=unittest.TextTestRunner(verbosity=1).run(suite)
raise SystemExit(not result.wasSuccessful())
