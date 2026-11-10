from pathlib import Path
import sys,json,time,numpy as np,traceback
root=Path(__file__).parent;sys.path.insert(0,str(root/'stage'))
from lapw_native_potential_export import prepare_native_potential_core
from lapw_potential_core_state import load_potential_core_state
from lapw_density_functional import PreparedDensityFunctional
from frozen_sigma_klein import last_tag_value
from lapw_hydrostatic_potential import hydrostatic_explicit_potential
results={}
for mat,d,c in [('MnTe','dmft','0_800_MnTe'),('MnO','MnO_center','MnO')]:
 w=root/'native'/mat/d;t=time.perf_counter()
 try:
  prepare_native_potential_core(w,c,w/'potential_core.npz')
  pc=load_potential_core_state(w/'potential_core.npz')
  scalar=PreparedDensityFunctional(pc)
  ref=last_tag_value(w/(c+'.scf0'),':DEN')
  row=dict(value_ry=scalar.value_ry,native_ry=ref,error_ry=scalar.value_ry-ref,transpose_error=pc.potential_operator.check_transpose(samples=1),volume=pc.hartree.volume)
 except Exception as exc:
  row=dict(error=type(exc).__name__+': '+str(exc));traceback.print_exc()
 row['seconds']=time.perf_counter()-t;results[mat]=row;print(mat,json.dumps(row),flush=True)
 (root/'potential_probe.json').write_text(json.dumps(results,indent=2)+'\n')
