"""Material smoke test using staged helpers plus dbin dependencies only."""
from pathlib import Path
import sys,json
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
out={}
for case in ['FeSe','SrVO3']:
 with load_native_valence_state(root/case/'charge/0',case,complex_vectors=case=='SrVO3',energy_policy='fixed',library='/Users/haule/dbin/liblapw_response.so') as f:
  pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
  select,coords=pc.hartree.native_potential_output(tuple(f.bands.spherical.slices),{a:op.layout.harmonics for a,(op,ks) in f.bands.ns.items()},f.bands.ir[0].representatives,complex_interstitial=case=='SrVO3')
  action=(f.valence_operator(fixed_number=True)@select+pc.core.operator)@pc.potential_operator
  direction=np.load('/private/tmp/edmft_hartree_periodic_20261108/'+case+'_direction.npz')['density']
  actual=action.matvec(direction)
  expected=np.load(root/(case+'_fixed_n_compact.npz'))['density']
  error=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected))
  assert error<2e-11
  out[case]=dict(relative_to_repository_dependencies=error)
  print(case,error,flush=True)
(root/'dbin_material.json').write_text(json.dumps(out,indent=2)+'\n')
