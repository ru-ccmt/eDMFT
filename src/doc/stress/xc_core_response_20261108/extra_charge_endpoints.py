from pathlib import Path
import os,sys,shutil,subprocess
import numpy as np
root=Path(__file__).resolve().parent
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from wien_response_density import Density
case=sys.argv[1];base=root/case/'xc/0';den=Density(base/(case+'.clmsum'),base/(case+'.struct'))
d=np.load('/private/tmp/edmft_hartree_periodic_20261108/'+case+'_direction.npz')['density']
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
for step in [.001,-.001,.0005,-.0005]:
 folder=root/case/'xc'/str(step);folder.mkdir(parents=True,exist_ok=True)
 for name in ['lapw0.def','xc_funcs.h',case+'.struct',case+'.in0',case+'.inm']:shutil.copy2(base/name,folder/name)
 den.from_vector(den.vector()+step*d).write(folder/(case+'.clmsum'))
 with (folder/'run.log').open('w') as log:subprocess.run([str(root/'native/SRC_lapw0/lapw0'),'lapw0.def'],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 core=root/case/'hxc_core'/str(step);core.mkdir(parents=True,exist_ok=True)
 for name in [case+'.struct',case+'.vsp',case+'.vns']:shutil.copy2(folder/name,core/name)
 for name in [case+'.inc','lcore.def']:shutil.copy2(root/case/'core/0'/name,core/name)
 (core/'core_native.bin').unlink(missing_ok=True)
 with (core/'run.log').open('w') as log:subprocess.run([str(root/'native/SRC_lcore/lcore'),'lcore.def'],cwd=core,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 print(case,step,'potential/core complete',flush=True)
