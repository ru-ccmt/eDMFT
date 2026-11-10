"""APW boundary matrix versus independent angular quadrature, complex density."""
from pathlib import Path
import re,subprocess,ctypes as ct,json
import numpy as np
from scipy.special import sph_harm_y
ROOT=Path(__file__).resolve().parent;DFT=Path('/Users/haule/wien2k_24/SRC_lapw2');DMFT=Path('/Users/haule/projects/eDMFT_new/src/dmft2')
L=3;NH=16;NB=4;R=1.7
stub='''module param
integer,parameter::lmax2=3,nume=4,nloat=2,nrad=41
end module
module atspdt
real*8::p(0:3),pe(0:3),dp(0:3),dpe(0:3),el(0:3)=0
end module
module lo
integer::ilo(0:3)=2
logical::loor(2,0:3)=.true.
real*8::plo(2,0:3),dplo(2,0:3),a1lo(41,2,0:3)=0
end module
module xa
real*8::weight(4)=1
end module
module radial_functions
real*8::VR(41)=0
end module
module struk
integer::jri(1)=41
real*8::r0(1)=0.1d0,dx(1)=0.07d0
end module
subroutine setup(val,der)
use atspdt
use lo
real*8::val(4,0:3),der(4,0:3)
p=val(1,:);pe=val(2,:);dp=der(1,:);dpe=der(2,:)
plo=val(3:4,:);dplo=der(3:4,:)
end subroutine
'''
rng=np.random.default_rng(8182)
val=rng.normal(size=(4,L+1));der=rng.normal(size=(4,L+1))
a=rng.normal(size=(NH,NB))+1j*rng.normal(size=(NH,NB));b=rng.normal(size=(NH,NB))+1j*rng.normal(size=(NH,NB));c=rng.normal(size=(NH,NB,2))+1j*rng.normal(size=(NH,NB,2));ds=rng.normal(size=(NH,NB))+1j*rng.normal(size=(NH,NB))
rho=rng.normal(size=(NB,NB))+1j*rng.normal(size=(NB,NB));rho=(rho+rho.conj().T)/2
w,u=np.linalg.eigh(rho)
pp=a@u;bb=b@u;cc=np.einsum('ink,nm->imk',c,u);ss=ds@u
states=[(l,m) for l in range(L+1) for m in range(-l,l+1)]
x,wt=np.polynomial.legendre.leggauss(22);phi=np.arange(44)*2*np.pi/44
x,phi=np.meshgrid(x,phi,indexing='ij');qw=np.repeat(wt,44)*2*np.pi/44
nv=np.stack([np.sqrt(1-x*x)*np.cos(phi),np.sqrt(1-x*x)*np.sin(phi),x],axis=-1).reshape(-1,3)
y=np.array([sph_harm_y(l,m,np.arccos(nv[:,2]),np.arctan2(nv[:,1],nv[:,0])) for l,m in states]).T
psi=np.zeros_like(a);dm=np.zeros_like(a)
for i,(l,m) in enumerate(states):
 psi[i]=val[0,l]*a[i]+val[1,l]*b[i]+c[i]@val[2:,l]
 dm[i]=der[0,l]*a[i]+der[1,l]*b[i]+c[i]@der[2:,l]
psi=y@psi;jump=y@(dm-ds);ref=np.zeros((3,3),complex)
for i in range(3):
 for j in range(3):
  mat=R**2*psi.conj().T@((qw*nv[:,i]*nv[:,j])[:,None]*jump)
  ref[i,j]=np.trace(rho@mat)
results={}
for name in ['lapw2','dmft2']:
 folder=ROOT/('surface_'+name);folder.mkdir(exist_ok=True)
 if name=='lapw2':body=(DFT/'stress_APW_surface.f').read_text()
 else:body=re.search(r'(?ims)^subroutine stress_APW_surface\(.*?^end subroutine\s*$',(DMFT/'stress_kernels.f90').read_text()).group()
 src=folder/'probe.f90';src.write_text(stub+body)
 more=['c_alpha_m.f','gaunt.f','t3j.f','t3j0.f','notri.f','ph.f']+(['spline.f','dergl.f'] if name=='lapw2' else [])
 subprocess.run(['/opt/homebrew/bin/gfortran','-O0','-g','-fcheck=all','-ffree-form','-ffree-line-length-none','-shared','-fPIC',str(src),*[str(DFT/f) for f in more],'-o',str(folder/'probe.dylib')],cwd=folder,check=True)
 lib=ct.CDLL(str(folder/'probe.dylib'))
 def call(fun,*args):getattr(lib,fun)(*[np.asfortranarray(arg).ctypes.data_as(ct.c_void_p) for arg in args])
 call('setup_',val,der)
 out=np.zeros((3,3,NB),complex,order='F');d1=out.copy(order='F');d2=out.copy(order='F')
 args=[np.array([1],np.int32),np.array([NB],np.int32),np.array([1],np.int32),np.array([R]),pp,bb,cc,ss]
 if name=='lapw2':args+=[np.zeros_like(ss),out,d1,d2]
 else:args+=[out]
 call('stress_apw_surface_',*args)
 got=np.einsum('n,ijn->ij',w,out);error=float(abs(got-ref).max());assert error<2e-10,(name,error)
 results[name]={'complex_matrix_contraction_error':error}
(ROOT/'surface_results.json').write_text(json.dumps(results,indent=2)+'\n');print(results)
