"""Actual reconstruction routine on a rank-deficient U_W with exact static moments."""
import ctypes as ct,json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
import os
p=Path(os.environ['EDMFT_STRESS_TEST_DIR']);lib=ct.CDLL(str(p/'audit_routines.dylib'))
F=lambda a:np.asfortranarray(a);I=lambda x:F(np.array([x],np.int32));Z=lambda s:np.zeros(s,complex,order='F')
nb,nw,nc=6,3,2;energies=np.array([-1.4,-.9,-.25,.4,1.1,1.8])
rng=np.random.default_rng(214);u=(rng.normal(size=(nb,nc))+1j*rng.normal(size=(nb,nc)))*.25
u[:nw,1]=0;sig=np.array([[.3,.2+.1j],[.2-.1j,-.1]])
eff=np.diag(energies[:nw])+u[:nw]@sig@u[:nw].conj().T
val,v=eigh(eff);temp=.23;mu=-.5;rho=(v/(1+np.exp((val-mu)/temp)))@v.conj().T;rhoe=eff@rho
cdirect=sig@u[:nw].conj().T@rho
rr=Z((nb,nb));er=Z((nb,nb));g=Z((nw,nc,1))
args=[rr,er,g,F(u[:,:,None]),F(rho),F(rhoe),F(energies),*[I(x) for x in [1,nb,1,nw,nb,nw,nc,1,nc,1]],I(1),I(nc),I(nc),F(np.arange(1,nc+1,dtype=np.int32)[:,None]),F(np.ones((nc,nc,1),np.int32)),F(cdirect.T[:,:,None]),I(1)]
lib.buildwindowresponsematrices_(*[a.ctypes.data_as(ct.c_void_p) for a in args])
cfit=g[:,:,0].T;dh=Z((nb,nb));dh[nw:,:nw]=(rng.normal(size=(nb-nw,nw))+1j*rng.normal(size=(nb-nw,nw)))*.2;dh[:nw,nw:]=dh[nw:,:nw].conj().T
rot=dh[nw:,:nw]/(energies[:nw][None,:]-energies[nw:,None]);du=rot.conj().T@u[nw:]
def fun(t):
 en,vv=eigh(np.diag(energies)+t*dh);pp=vv[:,:nw].conj().T@u
 vals=eigh(np.diag(en[:nw])+pp@sig@pp.conj().T,eigvals_only=True)
 return -temp*np.logaddexp(0,-(vals-mu)/temp).sum()
r={'rank_U_window':int(np.linalg.matrix_rank(u[:nw])),'orbital_dimension':nc,'B_reconstruction_residual':float(np.max(abs(u[:nw]@cfit-(rhoe-np.diag(energies[:nw])@rho)))),'C_error':float(np.max(abs(cfit-cdirect))),'functional_fd':float((fun(1e-5)-fun(-1e-5))/2e-5),'direct_C_contraction':float(2*np.trace(cdirect@du).real),'production_C_contraction':float(2*np.trace(cfit@du).real)}
(p/'rank_deficient_C.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

assert r['C_error']<1e-12,r
assert abs(r['production_C_contraction']-r['functional_fd'])<1e-9,r
print('PASS: direct C retains the rank-deficient null-space component')
