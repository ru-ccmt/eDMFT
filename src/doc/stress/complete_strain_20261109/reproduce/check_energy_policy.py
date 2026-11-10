from pathlib import Path
import sys,json,numpy as np
root=Path(__file__).resolve().parent;sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_radial_energy_policy import searched_energy_derivatives
from lapw_potential_core_state import load_potential_core_state
case=sys.argv[1];pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
source=np.load(root/case/'current_readout/source.npz');sol=np.load(root/case/'current_readout/solution.npz')
v=source['potential_source']+sol['potential'];A,detail=searched_energy_derivatives(root/case/'current',case,complex_vectors=case=='SrVO3',library='/Users/haule/dbin/liblapw_response.so')
fd=json.loads((root/case/'energy_parameters/REPORT.json').read_text())['derivative'];out={}
for (a,l),mat in A.items():
 sl=next(sl for key,sl in zip(pc.hartree.keys,pc.hartree.slices) if key==(a,0,0))
 de=mat@v[sl];out[f'{a},{l}']=dict(predicted_fixed_energy_response=de.tolist(),native_full_response=fd[f'{a},{l}'])
print(json.dumps(dict(root_checks=detail,derivatives=out),indent=2),flush=True)
(root/case/'energy_policy.json').write_text(json.dumps(dict(root_checks=detail,derivatives=out),indent=2)+'\n')
np.savez_compressed(root/case/'energy_policy.npz',**{f'{a},{l}':v for (a,l),v in A.items()})
