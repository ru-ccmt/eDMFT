"""Independent finite-H/O test of the actual G1 and complete G2 kernels.
Regular radial solutions all solve the SAME harmonic-oscillator radial ODE.
Compact C1 combinations remove surface terms without on-shell band assumptions.
"""
from pathlib import Path
import os,re,subprocess,ctypes as ct,json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import null_space
from scipy.special import sph_harm_y
ROOT=Path(__file__).resolve().parent
SRC=Path('/Users/haule/projects/eDMFT_new/src/dmft2')
N=1601;L=2;NB=4;R=2.6;dx=.006;r0=R*np.exp(-(N-1)*dx);mesh=r0*np.exp(dx*np.arange(N));NH=(L+1)**2;NP=4*NH
# Frobenius expansion for u=r^l v, v''+2(l+1)v'/r+(E-r^2)v=0.
# w=du/dE_Ry solves the same ODE with the inhomogeneity -u.
def radial(l,E,rs):
 def rhs(r,y):
  v,dv,w,dw=y
  return [dv,-2*(l+1)*dv/r+(r*r-E)*v,dw,-2*(l+1)*dw/r+(r*r-E)*w-v]
 eps=1e-7;den=2*(2*l+3)
 y0=[1-E*eps**2/den,-2*E*eps/den,-eps**2/den,-2*eps/den]
 sol=solve_ivp(rhs,(eps,R),y0,method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
 assert sol.success
 v,dv,w,dw=sol.sol(rs)
 u=rs**l*v;du=l*rs**(l-1)*v+rs**l*dv
 ud=rs**l*w;dud=l*rs**(l-1)*w+rs**l*dw
 ddu=-2*du/rs+(l*(l+1)/rs**2+rs**2-E)*u
 ddud=-2*dud/rs+(l*(l+1)/rs**2+rs**2-E)*ud-u
 return (u,du,ddu),(ud,dud,ddud)
energies=np.zeros((4,L+1),order='F')
rg,wg=np.polynomial.legendre.leggauss(180);rg=(rg+1)*R/2;wg=wg*R/2
f=np.zeros((N,4,L+1),order='F');df=f.copy(order='F')
g=np.zeros((len(rg),4,L+1));dg=g.copy();ddg=g.copy();end=np.zeros((2,4,L+1))
for l in range(L+1):
 E=2*l+3
 for p,e in [(0,E),(2,E+4),(3,E+8)]:
  energies[p,l]=e
  a,b=radial(l,e,mesh);q,qd=radial(l,e,rg);z,zd=radial(l,e,np.array([R]))
  for pp,aa,qq,zz in [(p,a,q,z)]+([(1,b,qd,zd)] if p==0 else []):
   f[:,pp,l]=aa[0];df[:,pp,l]=aa[1];g[:,pp,l]=qq[0];dg[:,pp,l]=qq[1];ddg[:,pp,l]=qq[2];end[:,pp,l]=[zz[0][0],zz[1][0]]
 energies[1,l]=E
rad=np.zeros((4,L+1,4,L+1,3),order='F');srad=np.zeros((4,4,L+1),order='F')
for i in range(4):
 for l in range(L+1):
  for j in range(4):
   for lp in range(L+1):
    rad[i,l,j,lp]=[np.dot(wg*rg**2*g[:,i,l],ddg[:,j,lp]),np.dot(wg*rg*g[:,i,l],dg[:,j,lp]),np.dot(wg*g[:,i,l],g[:,j,lp])]
   srad[i,j,l]=np.dot(wg*rg**2*g[:,i,l],g[:,j,l])
# Compile unchanged production bodies with small standalone module fixtures.
fixture=(SRC/'tests/g2_fixture.f90').read_text().replace('nrad=801','nrad='+str(N)).replace('ncom=2','ncom=2,nume=4,lmax2n=lmax2')
fixture=fixture.replace('real(8) :: el(0:lmax2)','real(8) :: el(0:lmax2),p(0:lmax2)=0,dp(0:lmax2)=0,pe(0:lmax2)=0,dpe(0:lmax2)=0')
extra='''module defs
complex*16,parameter::zeroc=(0d0,0d0)
end module
module stress_kb
use param
'''
for q in range(1,4):
 for pair in ['00','01','10','11']:extra+=f'real*8::tc{q}s{pair}(0:lmax2,0:lmax2)\n'
for pre in ['juddu','judu','juu']:
 for pair in ['02','12','20','21']:extra+=f'real*8::{pre}_{pair}(nloat,0:lmax2,0:lmax2)\n'
 extra+=f'real*8::{pre}_22(nloat,nloat,0:lmax2,0:lmax2)\n'
extra+='end module\n'
s=(SRC/'stress_kernels.f90').read_text()
def routine(name):return re.search(rf'(?ims)^\s*subroutine {name}\(.*?^\s*end subroutine {name}\b',s).group()+'\n'
prod=(SRC/'stress_g2_matrices.f90').read_text()
wrap='''subroutine initialize(f,df,en,rad,ov,rzero,step)
use param
use structure
use lo
use atspdt,only:el
use w_atpar
use stress_kb
use stress_g2_matrices
implicit none
real*8::f(nrad,4,0:lmax2),df(nrad,4,0:lmax2),en(4,0:lmax2),rad(4,0:lmax2,4,0:lmax2,3),ov(4,4,0:lmax2),rzero,step,r(nrad)
integer::i,j,l,lp,q
r0=rzero;dx=step
rotloc=0;do i=1,3;rotloc(i,i,1)=1;enddo
do i=1,nrad;r(i)=rzero*exp((i-1)*step);enddo
w_lm(:,1,1)=[0,0];w_lm(:,2,1)=[2,0]
w_rrad2=0;w_rade2=0;w_b1lo=0;w_radp2=0;w_radpe2=0;w_bp1lo=0
do l=0,lmax2
el(l)=en(1,l);elo(l,:)=en(3:4,l)
w_rrad1(:,l,1)=r*f(:,1,l);w_rade1(:,l,1)=r*f(:,2,l)
w_radp1(:,l,1)=r*df(:,1,l);w_radpe1(:,l,1)=r*df(:,2,l)
do i=1,2
w_a1lo(:,i,l,1)=r*f(:,i+2,l);w_ap1lo(:,i,l,1)=r*df(:,i+2,l)
enddo
ri_mat(:,:,l,1)=ov(:,:,l)
do lp=0,lmax2
'''
for q in range(1,4):
 for i,j in [(0,0),(0,1),(1,0),(1,1)]:wrap+=f'tc{q}s{i}{j}(l,lp)=rad({i+1},l,{j+1},lp,{q})\n'
for q,pre in enumerate(['juddu','judu','juu'],1):
 wrap+=f'''do i=1,nloat
{pre}_20(i,l,lp)=rad(i+2,l,1,lp,{q})
{pre}_21(i,l,lp)=rad(i+2,l,2,lp,{q})
{pre}_02(i,lp,l)=rad(1,l,i+2,lp,{q})
{pre}_12(i,lp,l)=rad(2,l,i+2,lp,{q})
do j=1,nloat
{pre}_22(i,j,l,lp)=rad(i+2,l,j+2,lp,{q})
enddo
enddo
'''
wrap+='''enddo
enddo
call prepare_g2_matrices()
end subroutine
subroutine contract(coeff,deriv,weight,energy_part,out,parts)
use param
use stress_g2_matrices
implicit none
integer,parameter::nh=(lmax2+1)**2,np=4*nh
complex*16::coeff(np,nume),deriv(np,nume,6)
real*8::weight(nume),out(3,3),parts(3,3,3)
integer::energy_part
complex*16::a(nh,nume),b(nh,nume),c(nh,nume,nloat),sa(nh,nume,6),sb(nh,nume,6),sc(nh,nume,nloat,6)
complex*16::sph(3,3,nume),ns0(3,3,nume),nsr(3,3,nume),nsa(3,3,nume),sp0(3,3,nume),spr(3,3,nume),spa(3,3,nume),kin(3,3,nume),g20(3,3,nume),g21(3,3,nume)
real*8::eb(nume)
integer::l,p,m,i,ly,iq,k
i=0
do l=0,lmax2
do p=1,4
do m=-l,l
i=i+1;ly=l*(l+1)+m+1
select case(p)
case(1)
a(ly,:)=coeff(i,:);sa(ly,:,:)=deriv(i,:,:)
case(2)
b(ly,:)=coeff(i,:);sb(ly,:,:)=deriv(i,:,:)
case default
c(ly,:,p-2)=coeff(i,:);sc(ly,:,p-2,:)=deriv(i,:,:)
end select
enddo
enddo
enddo
do iq=1,3
sa(:,:,iq)=sa(:,:,iq)-0.5d0*a;sb(:,:,iq)=sb(:,:,iq)-0.5d0*b;sc(:,:,:,iq)=sc(:,:,:,iq)-0.5d0*c
enddo
eb=0
call evaluate_g2(1,1,nume,a,b,c,sa,sb,sc,eb,sph,ns0,nsr,nsa,sp0,spr,spa,.true.)
g20=sph+ns0+nsr+nsa
parts=0
if(energy_part==0) then
call val_kinstrmt(1,nume,a,b,c,kin,2.6d0)
do k=1,nume
parts(:,:,1)=parts(:,:,1)+2*weight(k)*real(kin(:,:,k))
parts(:,:,2)=parts(:,:,2)+2*weight(k)*real(g20(:,:,k))
enddo
else
eb=1
call evaluate_g2(1,1,nume,a,b,c,sa,sb,sc,eb,sph,ns0,nsr,nsa,sp0,spr,spa,.true.)
g21=sph+ns0+nsr+nsa
do k=1,nume
parts(:,:,3)=parts(:,:,3)-2*weight(k)*real(g20(:,:,k)-g21(:,:,k))
enddo
endif
out=sum(parts,dim=3)
end subroutine
'''
(ROOT/'combined_fixture.f90').write_text(fixture+'\n'+extra+prod+'\n'+routine('val_kinstrmt')+wrap)
subprocess.run(['/opt/homebrew/bin/gfortran','-O0','-g','-fcheck=all','-ffree-form','-ffree-line-length-none','-shared','-fPIC',str(ROOT/'combined_fixture.f90'),str(SRC/'ylm.f90'),str(SRC/'atpar/rint13.f90'),*[str(Path('/Users/haule/wien2k_24/SRC_lapw2')/f) for f in ['c_alpha_m.f','gaunt.f','t3j.f','t3j0.f','notri.f','ph.f']],'-framework','Accelerate','-o',str(ROOT/'combined.dylib')],cwd=ROOT,check=True)
lib=ct.CDLL(str(ROOT/'combined.dylib'))
def call(fun,*args):getattr(lib,fun)(*[a.ctypes.data_as(ct.c_void_p) for a in args])
call('initialize_',f,df,energies,rad,srad,np.array([r0]),np.array([dx]))
# Independent primitive S,T, and Cartesian quadratic potential matrices.
states=[(l,m,p) for l in range(L+1) for p in range(4) for m in range(-l,l+1)]
S=np.zeros((NP,NP));T=S.copy();Q=np.zeros((NP,NP,3,3),complex)
x,ww=np.polynomial.legendre.leggauss(24);phi=np.arange(48)*2*np.pi/48
xx,pp=np.meshgrid(x,phi,indexing='ij');aw=np.repeat(ww,48)*2*np.pi/48
nv=np.stack([np.sqrt(1-xx**2)*np.cos(pp),np.sqrt(1-xx**2)*np.sin(pp),xx],axis=-1).reshape(-1,3)
y=np.array([sph_harm_y(l,m,np.arccos(nv[:,2]),np.arctan2(nv[:,1],nv[:,0])) for l,m,p in states]).T
for i,(l,m,p) in enumerate(states):
 for j,(lp,mp,q) in enumerate(states):
  if (l,m)==(lp,mp):
   S[i,j]=np.dot(wg*rg**2*g[:,p,l],g[:,q,lp])
   T[i,j]=np.dot(wg*rg**2*dg[:,p,l],dg[:,q,lp])+l*(l+1)*np.dot(wg*g[:,p,l],g[:,q,lp])
  rr=np.dot(wg*rg**4*g[:,p,l],g[:,q,lp])
  for a in range(3):
   for b in range(3):Q[i,j,a,b]=rr*np.dot(aw*nv[:,a]*nv[:,b]*y[:,i].conj(),y[:,j])
W=np.eye(3)+.31*np.sqrt(5/(16*np.pi))*np.diag([-1,-1,2])
V=np.einsum('ijab,ab->ij',Q,W);H=T+V
# C1 compact combinations in the nullspace of value AND slope at R.
Z=np.zeros((NP,2*NH))
for l in range(L+1):
 z=null_space(end[:,:,l])
 for m in range(-l,l+1):
  ids=[states.index((l,m,p)) for p in range(4)];ly=l*(l+1)+m
  Z[np.ix_(ids,[2*ly,2*ly+1])]=z
rng=np.random.default_rng(318431)
C=Z@(rng.normal(size=(2*NH,NB))+1j*rng.normal(size=(2*NH,NB)))
e,u=np.linalg.eigh(C.conj().T@S@C);C=C@(u/np.sqrt(e))@u.conj().T
ks,u=np.linalg.eigh(C.conj().T@H@C);C=C@u
axes=[(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]
D=np.stack([Z@(rng.normal(size=(2*NH,NB))+1j*rng.normal(size=(2*NH,NB)))*.07 for _ in axes],axis=-1)
rho=rng.normal(size=(NB,NB))+1j*rng.normal(size=(NB,NB));rho=.15*(rho+rho.conj().T)+np.eye(NB)
rhoE=rng.normal(size=(NB,NB))+1j*rng.normal(size=(NB,NB));rhoE=.3*(rhoE+rhoE.conj().T)

def analytic(rho,rhoE,D):
 total=np.zeros((3,3));part=np.zeros((3,3,3))
 for i,M in enumerate([rho,rhoE]):
  w,u=np.linalg.eigh(M);cc=np.asfortranarray(C@u);dd=np.asfortranarray(np.einsum('inb,nm->imb',D,u))
  out=np.zeros((3,3),order='F');parts=np.zeros((3,3,3),order='F')
  call('contract_',cc,dd,w,np.array([i],np.int32),out,parts);total+=out;part+=parts
 return total,part

def value(t,E,Dq,rho,rhoE):
 F=np.eye(3)+t*E;Fi=np.linalg.inv(F);J=np.linalg.det(F);A=C+t*Dq
 VV=np.einsum('ijab,ab->ij',Q,Fi.T@W@Fi)
 return (np.trace(rho@(A.conj().T@(T+VV)@A))-np.trace(rhoE@(A.conj().T@S@A))).real/J
results={}
for name,r,rE,der in [('DMFT_noncommuting',rho,rhoE,D),('DFT_diagonal',np.diag([.9,.7,.3,.1]),np.diag(ks*np.array([.9,.7,.3,.1])),D),('geometric_cancellation',rho,rhoE,np.zeros_like(D))]:
 an,parts=analytic(r,rE,der);fd=np.zeros((3,3));h=2e-5
 for q,(a,b) in enumerate(axes):
  E=np.zeros((3,3));E[a,b]+=.5;E[b,a]+=.5
  fd[a,b]=fd[b,a]=(value(h,E,der[:,:,q],r,rE)-value(-h,E,der[:,:,q],r,rE))/(2*h)
 error=float(abs(an-fd).max());results[name]={'max_error':error,'analytic':an.tolist(),'finite_difference':fd.tolist(),'contributions':parts.tolist()}
 print(name,error,flush=True)
 assert error<3e-5,(name,error,an,fd)
(ROOT/'combined_results.json').write_text(json.dumps(results,indent=2)+'\n')
print('PASS combined G1+G2 full-H/O finite difference')
