import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys,json
sys.dont_write_bytecode=True
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from check_projector_sternheimer_fd import read_dump,read_normalization
r=Path(__file__).parent/'FeSe';m=np.load(r/'matrices.npz');D=m['total'];R=m['response']+m['norm_response'];target=float(np.einsum('kij,kji->',R,D).real);fd={}
for h in [1e-5,5e-6]:
 val=0.
 for k,V in enumerate(D,1):
  d=read_dump(r/'center'/f'dmft_sternheimer.k{k:08d}.bin');n=read_normalization(r/'center'/f'dmft_sternheimer_norm.k{k:08d}.bin');W=np.arange(d.dm_nemin-d.nemin,d.dm_nemax-d.nemin+1);contraction=[]
  for sign in [1,-1]:
   e,v=eigh(np.diag(d.energies)+sign*h*V);vw=v[:,W];u=vw@(vw.conj().T@d.projector.reshape(d.nstored,-1));p=vw@(vw.conj().T@n.projector.reshape(d.nstored,-1))
   contraction.append(2*np.real(np.sum(d.gradient*u[W].reshape(d.gradient.shape))+np.sum(n.gradient*p[W].reshape(n.gradient.shape))))
  val+=(contraction[0]-contraction[1])/(2*h)
 fd[str(h)]={'derivative_Ry':float(val),'difference_Ry':float(val-target)}
result={'analytic_Ry':target,'finite_differences':fd};print(json.dumps(result,indent=2));(r/'projector_fd_checks.json').write_text(json.dumps(result,indent=2))
