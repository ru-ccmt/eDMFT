from pathlib import Path
import sys,os,json,shutil,subprocess,time,numpy as np
root=Path(__file__).parent;sys.path.insert(0,str(root/'stage'))
from lapw_potential_core_state import load_potential_core_state
from lapw_density_functional import PreparedDensityFunctional
from lapw_charge_response import native_charge_covectors
from wien_response_density import Density
from frozen_sigma_klein import last_tag_value

env=dict(os.environ,WIENROOT='/Users/haule/wien2k_24',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
results={}
for mat,d,c in [('MnTe','dmft','0_800_MnTe'),('MnO','MnO_center','MnO')]:
 w=root/'native'/mat/d; pc=load_potential_core_state(w/'potential_core.npz'); f=PreparedDensityFunctional(pc)
 den=Density(w/(c+'.clmsum'),w/(c+'.struct')); n=den.vector(); delta=den.copy()
 for (a,l,m),v in delta.rad.items(): v*=.1*np.exp(-np.exp(delta.grids[a]))*(1+.2*abs(l))
 for g in delta.pw: delta.pw[g]*=.08
 direction=delta.vector(); mt,ir=native_charge_covectors(den,pc.structure,inversion=pc.metadata['inversion']); charge=mt+ir
 direction-=n*(charge@direction)/(charge@n)
 analytic=f.action(direction); values={}
 for h in [.001,.0005]:
  endpoints=[]
  for sign in [-1,1]:
   dest=root/'native_density_probe'/mat/(str(h)+'_'+str(sign));dest.mkdir(parents=True,exist_ok=True)
   for suf in ['struct','in0','inm','clmsum']:
    shutil.copy2(w/(c+'.'+suf),dest/(c+'.'+suf))
   shutil.copy2(w/'lapw0.def',dest/'lapw0.def')
   den.from_vector(n+sign*h*direction).write(dest/(c+'.clmsum'))
   with (dest/'native.log').open('w') as log:
    done=subprocess.run(['/Users/haule/dbin/lapw0_response','lapw0.def'],cwd=dest,env=env,stdout=log,stderr=subprocess.STDOUT)
   if done.returncode or 'LAPW0 END' not in ' '.join((dest/'native.log').read_text().split()): raise RuntimeError(str(dest))
   endpoints.append(last_tag_value(dest/(c+'.scf0'),':DEN'))
  values[str(h)]={'finite_difference_Ry':(endpoints[1]-endpoints[0])/(2*h),'error_Ry':(endpoints[1]-endpoints[0])/(2*h)-analytic}
 results[mat]=dict(analytic_Ry=analytic,charge_direction=float(charge@direction),endpoints=values)
 print(mat,json.dumps(results[mat]),flush=True);(root/'native_density_probe.json').write_text(json.dumps(results,indent=2)+'\n')
