from pathlib import Path
from contextlib import ExitStack
import sys,time,json,traceback,numpy as np
root=Path(__file__).parent;sys.path.insert(0,str(root/'stage'))
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from lapw_magnetic_response import combine_valence
from wien_response_density import Density
results={}
for mat,d,c in [('MnTe','dmft','0_800_MnTe'),('MnO','MnO_center','MnO')]:
 t=time.perf_counter();w=root/'native'/mat/d
 try:
  pc=load_potential_core_state(w/'potential_core.npz');core=Density(w/(c+'.clmcor'),w/(c+'.struct'))
  with ExitStack() as stack:
   branches=[]
   for channel in (['','dn'] if mat=='MnTe' else ['']):
    v=stack.enter_context(load_native_valence_state(w,c,complex_vectors=False,energy_policy='fixed',matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so',diagnostic_gap_cutoff=1e-10,spin_channel=channel,cache_mb=128))
    branches.append(v);print(mat,channel,'valence prepared',v.preparation_seconds,flush=True)
   v=combine_valence(*branches) if len(branches)==2 else branches[0]
   state=NativeScreenedState(v,pc,core,target_charge=np.dot(pc.structure.mult,pc.structure.Znuc),fixed_number=False)
   row=dict(metadata=state.metadata,transpose=state.jacobian.check_transpose(samples=1,tolerance=2e-10))
 except Exception as exc: row=dict(error=type(exc).__name__+': '+str(exc));traceback.print_exc()
 row['seconds']=time.perf_counter()-t;results[mat]=row;print(mat,json.dumps(row),flush=True);(root/'valence_probe.json').write_text(json.dumps(results,indent=2)+'\n')
