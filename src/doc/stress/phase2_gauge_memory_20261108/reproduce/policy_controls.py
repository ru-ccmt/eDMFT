from pathlib import Path
import sys,json
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_valence_state import load_native_valence_state
from lapw_radial_native import NativeRadialContext
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python/tests')
from test_lapw_radial_native import inputs
import numpy as np
report={}
for case in ['FeSe','SrVO3']:
    source=Path('/private/tmp/edmft_full_k_response_20261108')/case/'charge/0'
    try:
        with load_native_valence_state(source,case,complex_vectors=case=='SrVO3',energy_policy='fixed',matrix_quadrature='atparn',library=str(root/'build_final/liblapw_response.so')):
            raise AssertionError('inconsistent policy was accepted')
    except ValueError as e:
        assert 'assembled overlap' in str(e)
        report[case]=str(e)
v,kw,r=inputs();kw['library']='/Users/haule/dbin/liblapw_response.so'
with NativeRadialContext(v,**kw) as old, NativeRadialContext(v,**(kw|dict(library=str(root/'build_final/liblapw_response.so')))) as new:
    report['legacy_primitive_max']=float(np.max(abs(old.primitive-new.primitive)))
    report['legacy_ABC_max']=float(np.max(abs(old.lo_coefficients-new.lo_coefficients)))
    report['legacy_ABC_relative']=float(np.linalg.norm(old.lo_coefficients-new.lo_coefficients)/max(1,np.linalg.norm(old.lo_coefficients)))
    assert report['legacy_primitive_max']<1e-13 and report['legacy_ABC_relative']<2e-13
try:
    NativeRadialContext(v,quadrature='wien2k24',**kw)
    raise AssertionError('old library silently accepted new policy')
except RuntimeError as e:
    assert 'updated response library' in str(e)
    report['old_library_new_policy']=str(e)
(root/'policy_controls.json').write_text(json.dumps(report,indent=2)+'\n')
print(report,flush=True)
