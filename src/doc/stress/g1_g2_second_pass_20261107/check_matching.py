"""FD of actual matching derivatives against direct endpoint Bessel/Ylm matching.
Tests APW, LAPW, two LOs per l, complex phases, general local rotations,
rectangular pseudo-band transformation, and occupied bands below DMFT window.
"""
from pathlib import Path
import re,ctypes as ct,subprocess,json,hashlib
import numpy as np
from scipy.special import spherical_jn,sph_harm_y
from scipy.spatial.transform import Rotation
ROOT=Path(__file__).resolve().parent;SRC=Path('/Users/haule/projects/eDMFT_new/src/dmft2')
L=2;NH=9;NPW=10;NLO=18;NB=5

def routine(path,name):return re.search(rf'(?ims)^\s*subroutine {name}\(.*?^\s*end subroutine {name}\b',path.read_text()).group()+'\n'
stubs='''module param
integer,parameter::lmax2=2,lomax=2,nume=5,nloat=2,iblock=4,nmat=10
end module
module defs
real*8,parameter::pi=3.1415926535897932384626433832795d0,two=2d0
complex*16,parameter::imag=(0d0,1d0)
end module
module xa
real*8::fj(0:2,10,1),dfj(0:2,10,1)
end module
module xa3
real*8::bk3(3,10),bk3lo(3,18)
complex*16::As(10,5,1),As_lo(18,5,1)
end module
module structure
real*8::rmt(1)=2d0,vol=200d0,rotij(3,3,1),br1(3,3),pos(3,1),tauij(3,1),rotloc(3,3,1)
integer::mult(1)=1
end module
module atspdt
real*8::p(0:2),dp(0:2),pe(0:2),dpe(0:2)
end module
module lo
integer::nlo=18,nlov=0,ilo(0:2)=2
logical::lapw(0:2)=[.true.,.false.,.true.]
real*8::alo(0:2,2),blo(0:2,2),clo(0:2,2)
end module
module dmfts
real*8::shft(1,3)
end module
'''
setter='''subroutine setup(q,qlo,c,cl,br,rt,rs,ps,ta,sh,boundary,locoeff,j,dj)
use structure
use xa
use xa3
use atspdt
use lo
use dmfts
implicit none
real*8::q(3,10),qlo(3,18),br(3,3),rt(3,3),rs(3,3),ps(3),ta(3),sh(3),boundary(3,4),locoeff(3,2,3),j(3,10),dj(3,10)
complex*16::c(10,5),cl(18,5)
bk3=q;bk3lo=qlo;As(:,:,1)=c;As_lo(:,:,1)=cl
br1=br;rotloc(:,:,1)=rt;rotij(:,:,1)=rs;pos(:,1)=ps;tauij(:,1)=ta;shft(1,:)=sh
p=boundary(:,1);dp=boundary(:,2);pe=boundary(:,3);dpe=boundary(:,4)
alo=locoeff(:,:,1);blo=locoeff(:,:,2);clo=locoeff(:,:,3);fj(:,:,1)=j;dfj(:,:,1)=dj
end subroutine
subroutine test(a,b,c,w,sa,sb,sc)
use alm_blm_clm
implicit none
complex*16::a(9,5),b(9,5),c(9,5,2),w(3,3),sa(9,5,6),sb(9,5,6),sc(9,5,2,6)
call cmp_alm_blm_clm_strain(sa,sb,sc,a,b,c,w,1,1,1,28,18,1,3,5,4,3,2,1)
end subroutine
'''
source=stubs+'module mod_lo\ncontains\n'+routine(SRC/'lomain.f90','lomain_strain')+'end module\n'
source+='module alm_blm_clm\ncontains\n'+routine(SRC/'alm_blm_clm.f90','cmp_alm_blm_clm_strain')+'end module\n'
source+=routine(SRC/'stress_kernels.f90','c_factor')+setter
(ROOT/'matching_fixture.f90').write_text(source)
subprocess.run(['/opt/homebrew/bin/gfortran','-O0','-g','-fcheck=all','-ffree-line-length-none','-shared','-fPIC',str(ROOT/'matching_fixture.f90'),str(SRC/'ylm.f90'),'-framework','Accelerate','-o',str(ROOT/'matching.dylib')],cwd=ROOT,check=True)
lib=ct.CDLL(str(ROOT/'matching.dylib'))
def call(name,*args):getattr(lib,name)(*[a.ctypes.data_as(ct.c_void_p) for a in args])
def arr(x):return np.asfortranarray(x)
rng=np.random.default_rng(93992)
q=rng.normal(size=(3,NPW));q[:,0]=0;qlo=rng.normal(size=(3,NLO))
c=rng.normal(size=(NPW,NB))+1j*rng.normal(size=(NPW,NB));cl=rng.normal(size=(NLO,NB))+1j*rng.normal(size=(NLO,NB))
br=np.array([[.65,.17,-.08],[.05,.82,.21],[.09,-.03,.74]])
rt=Rotation.from_rotvec([.43,-.31,.28]).as_matrix();rs=np.array([[0.,-1,0],[1,0,0],[0,0,1.]])
ps=np.array([.13,.28,.44]);ta=np.array([.5,0,.25]);sh=np.array([.04,.02,-.01])
boundary=rng.normal(size=(3,4));boundary[:,0]+=[2,3,4];locoeff=rng.normal(size=(3,2,3))
w=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
R=2.;vol=200.

def harmon(v,l):
 norm=np.linalg.norm(v,axis=0);z=np.divide(v[2],norm,out=np.ones_like(norm),where=norm>0)
 return np.array([sph_harm_y(l,m,np.arccos(np.clip(z,-1,1)),np.arctan2(v[1],v[0])) for m in range(-l,l+1)])
def coefficients(F):
 B=np.linalg.solve(F.T,br);J=np.linalg.det(F);fac=4*np.pi*R**2/np.sqrt(vol*J)
 v=B@rs@q;vloc=rt@v;mag=np.linalg.norm(v,axis=0)
 ph=np.exp(2j*np.pi*(ps@rs@q+ta@q+sh@q))
 a=np.zeros((NH,NB),complex);b=a.copy();lo=np.zeros((NH,NB,2),complex)
 for l in range(L+1):
  js=spherical_jn(l,mag*R);djs=mag*spherical_jn(l,mag*R,derivative=True)
  if l==1:ha=js/(R*R*boundary[l,0]);hb=np.zeros_like(ha)
  else:ha=djs*boundary[l,2]-js*boundary[l,3];hb=js*boundary[l,1]-djs*boundary[l,0]
  yp=harmon(vloc,l).conj()*ph
  a[l*l:(l+1)**2]=(yp*ha)@c;b[l*l:(l+1)**2]=(yp*hb)@c
  # LO indices follow l, jlo, m1, exactly the vector-file enumeration.
 offset=0
 for l in range(L+1):
  for j in range(2):
   ids=slice(offset,offset+2*l+1);offset+=2*l+1
   v=rt@B@rs@qlo[:,ids]
   ph=np.exp(2j*np.pi*(ps@rs@qlo[:,ids]+ta@qlo[:,ids]+sh@qlo[:,ids]))
   z=(harmon(v,l).conj()*ph)@cl[ids]
   a[l*l:(l+1)**2]+=locoeff[l,j,0]*z;b[l*l:(l+1)**2]+=locoeff[l,j,1]*z
   lo[l*l:(l+1)**2,:,j]+=locoeff[l,j,2]*z
 for l in range(L+1):
  a[l*l:(l+1)**2]*=fac*1j**l;b[l*l:(l+1)**2]*=fac*1j**l;lo[l*l:(l+1)**2]*=fac*1j**l
 for t in [a,b]:t[:,2:4]=t[:,2:5]@w[:,:2];t[:,4]=0
 for j in range(2):lo[:,2:4,j]=lo[:,2:5,j]@w[:,:2];lo[:,4,j]=0
 return a,b,lo

mag=np.linalg.norm(br@rs@q,axis=0)
j=np.array([spherical_jn(l,mag*R) for l in range(3)])
dj=np.array([mag*spherical_jn(l,mag*R,derivative=True) for l in range(3)])
call('setup_',*[arr(t) for t in [q,qlo,c,cl,br,rt,rs,ps,ta,sh,boundary,locoeff,j,dj]])
a,b,lo=map(arr,coefficients(np.eye(3)));sa=np.zeros((NH,NB,6),complex,order='F');sb=sa.copy(order='F');sc=np.zeros((NH,NB,2,6),complex,order='F')
call('test_',a,b,lo,arr(w),sa,sb,sc)
axes=[(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)];results={}
for h in [1e-5,3e-6]:
 dif=[np.zeros_like(t) for t in [sa,sb,sc]]
 for iq,(aa,bb) in enumerate(axes):
  E=np.zeros((3,3));E[aa,bb]+=.5;E[bb,aa]+=.5
  pp=coefficients(np.eye(3)+h*E);mm=coefficients(np.eye(3)-h*E)
  for dst,p,m in zip(dif,pp,mm):dst[...,iq]=(p-m)/(2*h)
 errors=[float(np.max(abs(t-d))) for t,d in zip([sa,sb,sc],dif)]
 results[str(h)]=dict(zip(['salm','sblm','sclm'],errors))
 print(h,errors)
 assert max(errors)<2e-7,errors
(ROOT/'matching_results.json').write_text(json.dumps(results,indent=2)+'\n')
print('PASS actual APW/LAPW/LO matching derivatives versus direct endpoints')
