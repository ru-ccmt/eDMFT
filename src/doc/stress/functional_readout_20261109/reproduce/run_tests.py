"""Run the suite using only installed source and installed test files."""
from pathlib import Path
import os
import re
import sys
import unittest

# Independent compiled native controls; override paths for another checkout.
for key,value in {
    'EDMFT_RADIAL_LIBRARY':'/Users/haule/dbin/liblapw_response.so',
    'EDMFT_MIXER_NORMALIZATION_LIBRARY':'/private/tmp/edmft_charge_closure_20261108/native/libmixer_normalization.dylib',
    'EDMFT_HARTREE_REFERENCE_LIBRARY':'/private/tmp/edmft_hartree_radial_20261108/native/libhartree_reference.dylib',
    'EDMFT_HARTREE_CHARGE_LIBRARY':'/private/tmp/edmft_hartree_periodic_20261108/native/libcharge_reference.dylib',
}.items():os.environ.setdefault(key,value)
checkpoint_root=Path(__file__).resolve().parent
if checkpoint_root.name=='reproduce': checkpoint_root=checkpoint_root.parent
repo = Path('/Users/haule/projects/eDMFT_new')
source = Path(os.environ['EDMFT_TEST_SOURCE'])
sys.path[:0] = [str(checkpoint_root / "stage"), str(checkpoint_root / "stage/tests"), str(source), str(repo / "src/python/tests")]
text = (repo / 'src/doc/stress/charge_closure_20261108/post_install_repository.log').read_text()
modules = sorted(set(re.findall(r'\((test_\w+)\.', text)) | {
    'test_lapw_compact_response', 'test_lapw_valence_state', 'test_lapw_radial_quadrature',
    'test_lapw_charge_response', 'test_lapw_potential_core_state',
    'test_lapw_core_response', 'test_lapw_xc_response', 'test_lapw_screened_state',
    'test_lapw_projector_functional', 'test_lapw_strain_potential', 'test_lapw_geometry_response', 'test_lapw_strain_charge', 'test_lapw_strain_field', 'test_lapw_strain_functional', 'test_lapw_physical_potential',
})
print('Installed source:', source, flush=True)
suite = unittest.defaultTestLoader.loadTestsFromNames(modules)
result = unittest.TextTestRunner(verbosity=2).run(suite)
assert result.testsRun == 225
assert result.wasSuccessful() and not result.skipped
helpers = ['lapw_strain_charge','lapw_strain_field','lapw_geometry_response','lapw_strain_potential','lapw_hydrostatic_potential','lapw_projector_functional', 'lapw_screened_state', 'lapw_compact_blocks', 'lapw_compact_response', 'lapw_native_snapshot',
           'lapw_valence_state', 'lapw_radial_native', 'lapw_spherical_response',
           'lapw_nonspherical_response', 'lapw_potential_response']
for name in helpers:
    actual = Path(sys.modules[name].__file__).resolve()
    expected = (source / (name + '.py')).resolve()
    assert actual == expected, (actual, expected)
for name in ['test_lapw_strain_charge','test_lapw_strain_field','test_lapw_geometry_response','test_lapw_strain_potential','test_lapw_projector_functional', 'test_lapw_screened_state', 'test_lapw_compact_response', 'test_lapw_valence_state', 'test_lapw_radial_quadrature']:
    assert Path(sys.modules[name].__file__).resolve() == (repo / 'src/python/tests' / (name + '.py')).resolve()
for name,sub in [('lapw_strain_functional','stage'),('test_lapw_strain_functional','stage/tests'),('lapw_physical_potential','stage'),('test_lapw_physical_potential','stage/tests')]:
    assert Path(sys.modules[name].__file__).resolve() == checkpoint_root / sub / (name+'.py')
print('Verified staged readout helper/test and installed dependencies.', flush=True)
