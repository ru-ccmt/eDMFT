from pathlib import Path
import re,subprocess,ctypes as ct,json
import numpy as np
from scipy.special import spherical_jn
root=Path(__file__).resolve().parent;dft=Path('/Users/haule/wien2k_24/SRC_lapw2');dmft=Path('/Users/haule/projects/eDMFT_new/src/dmft2')
ng=7;nb=4
stub=f'''module param
integer,parameter::nume={nb}
end module
module defs
real*8,parameter::pi=3.141592653589793238462643d0
complex*16,parameter::zeroc=(0d0,0d0)
end module
module structure
integer::mult(1)=1
real*8::pos(3,1),rmt(1),br1(3,3),vol
end module
module struk
use structure
end module
module com
integer::nat=1
end module
module xa
real*8::weight({nb})
end module
module xa3
complex*16,allocatable::a(:,:)
integer::k3(3,{ng}),kx({ng}),ky({ng}),kz({ng})
end module
module stress_pw
real*8::f1p({ng},3)
end module
'''
wrapper='''subroutine probe(g,c,q,w,p,r,b,v,out)
use param
use structure,only:pos,rmt,br1,vol
use xa,only:weight
use xa3,only:a,k3,kx,ky,kz
use stress_pw,only:f1p
implicit none
integer::g(3,7)
complex*16::c(7,4)
real*8::q(7,3),w(4),p(3),r,b(3,3),v,out(3,3)
k3=g;kx=g(1,:);ky=g(2,:);kz=g(3,:)
if(.not.allocated(a))allocate(a(7,4))
a=c;f1p=q;weight=w
pos(:,1)=p;rmt(1)=r;br1=b;vol=v
out=0
CALL_KERNEL
end subroutine
'''
def cx(s):return '\n'.join(line.replace('!_COMPLEX','') for line in s.splitlines() if not line.startswith('!_REAL'))
libs={}
for key in ['lapw2','dmft2']:
 folder=root/('interstitial_'+key);folder.mkdir(exist_ok=True)
 if key=='lapw2':
  bodies=cx((dft/'step_fun_kb.F').read_text())+cx((dft/'val_kin_interstitial.F').read_text())
  call='call val_kin_interstitial(7,out,1,4)'
 else:
  s=(dmft/'stress_interstitial.f90').read_text()
  bodies=re.search(r'(?ims)^SUBROUTINE step_fun\(.*?^END\s*$',s).group()+'\n'+re.search(r'(?ims)^subroutine val_kin_interstitial\(.*?^end subroutine\s*$',s).group()
  call='call val_kin_interstitial(7,out,1,4,w)'
 source=folder/'probe.f90';source.write_text(stub+bodies+'\n'+wrapper.replace('CALL_KERNEL',call))
 subprocess.run(['/opt/homebrew/bin/gfortran','-O0','-g','-fcheck=all','-ffree-form','-ffree-line-length-none','-shared','-fPIC',str(source),'-framework','Accelerate','-o',str(folder/'probe.dylib')],cwd=folder,check=True)
 libs[key]=ct.CDLL(str(folder/'probe.dylib'))
g=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[-1,1,0],[1,1,-1],[-1,0,1]],np.int32)
b=np.diag([1.1,.9,1.3]);v=(2*np.pi)**3/np.linalg.det(b);radius=.7;pos=np.array([.13,.27,.39]);k=np.array([.12,.07,.15]);q=(g+k)@b.T
rng=np.random.default_rng(1783);c=.2*(rng.normal(size=(ng,nb))+1j*rng.normal(size=(ng,nb)))
rho=rng.normal(size=(nb,nb))+1j*rng.normal(size=(nb,nb));rho=(rho+rho.conj().T)/2
w,u=np.linalg.eigh(rho);pseudo=c@u
dg=g[:,None,:]-g[None,:,:];norm=np.linalg.norm(dg@b.T,axis=2)*radius
j=np.ones_like(norm)/3;mask=norm>0;j[mask]=spherical_jn(1,norm[mask])/norm[mask]
theta=np.eye(ng)-(4*np.pi*radius**3/v)*j*np.exp(-2j*np.pi*np.einsum('ijx,x->ij',dg,pos))
reference=np.zeros((3,3))
for a in range(3):
 for bidx in range(3):
  ket=-c.conj().T@theta@(q[:,a,None]*q[:,bidx,None]*c)
  reference[a,bidx]=-2*np.trace(rho@ket).real
out={}
for key,lib in libs.items():
 result=np.zeros((3,3),order='F')
 args=[np.asfortranarray(g.T),np.asfortranarray(pseudo),np.asfortranarray(q),w,pos,np.array([radius]),np.asfortranarray(b),np.array([v]),result]
 lib.probe_(*[a.ctypes.data_as(ct.c_void_p) for a in args])
 err=float(abs(result-reference).max());assert err<1e-11,(key,err)
 out[key]={'max_error':err,'tensor':result.tolist()}
(root/'interstitial_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
