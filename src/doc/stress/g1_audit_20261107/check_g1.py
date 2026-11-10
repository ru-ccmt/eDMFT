"""Read-only audit: compile unchanged G1 routines, supply exact radial integrals.
The oracle differentiates Cartesian polynomials; it does not reuse Gaunt algebra.
"""
from pathlib import Path
import re,subprocess,ctypes as ct,json,os,hashlib
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parent
DFT=Path('/Users/haule/wien2k_24/SRC_lapw2')
DMFT=Path('/Users/haule/projects/eDMFT_new/src/dmft2')
L=3;NH=(L+1)**2;NP=4;NB=1
def routine(p,name):
    s=p.read_text()
    return re.search(rf'^\s*subroutine {name}\(.*?^\s*end subroutine {name}\b',s,re.I|re.S|re.M).group(0)+'\n'

stub=f'''module param
integer,parameter::lmax2={L},lmax2n={L},nume=1,nloat=2,lomax={L}
end module
module atspdt
use param
real*8::p(0:lmax2)=0,dp(0:lmax2)=0,pe(0:lmax2)=0,dpe(0:lmax2)=0
end module
module defs
complex*16,parameter::zeroc=(0d0,0d0)
end module
module lo
use param
integer::ilo(0:lmax2)=2
logical::loor(nloat,0:lmax2)=.true.
end module
module com
logical::rel=.false.
end module
module stress_kb
use param
'''
for q in range(1,4):
    for pair in ['00','01','10','11']:stub+=f'real*8::tc{q}s{pair}(0:lmax2,0:lmax2)\n'
for pre in ['juddu','judu','juu']:
    for pair in ['02','12','20','21']:stub+=f'real*8::{pre}_{pair}(nloat,0:lmax2,0:lmax2)\n'
    stub+=f'real*8::{pre}_22(nloat,nloat,0:lmax2,0:lmax2)\n'
stub+='end module\n'
setup='''subroutine fill_radial(rad)
use param
use stress_kb
implicit none
real*8,intent(in)::rad(4,0:lmax2,4,0:lmax2,3)
integer::l,lp,i,j
do l=0,lmax2
do lp=0,lmax2
'''
for q in range(1,4):
    for i,j in [(0,0),(0,1),(1,0),(1,1)]:setup+=f'tc{q}s{i}{j}(l,lp)=rad({i+1},l,{j+1},lp,{q})\n'
for q,pre in enumerate(['juddu','judu','juu'],1):
    setup+='do i=1,nloat\n'
    setup+=f'''{pre}_20(i,l,lp)=rad(i+2,l,1,lp,{q})
{pre}_21(i,l,lp)=rad(i+2,l,2,lp,{q})
{pre}_02(i,lp,l)=rad(1,l,i+2,lp,{q})
{pre}_12(i,lp,l)=rad(2,l,i+2,lp,{q})
do j=1,nloat
{pre}_22(i,j,l,lp)=rad(i+2,l,j+2,lp,{q})
enddo
enddo
'''
setup+='enddo\nenddo\nend subroutine\n'
libs={}
for key,path in [('lapw2',DFT/'val_kinstrmt.F'),('dmft2',DMFT/'stress_kernels.f90'),
                 ('lapw2_candidate',DFT/'val_kinstrmt.F'),('dmft2_candidate',DMFT/'stress_kernels.f90')]:
    folder=ROOT/key;folder.mkdir(exist_ok=True)
    body=routine(path,'val_kinstrmt')
    if key.endswith('_candidate'):
        # TEST COPY ONLY. Original source and installed binaries are untouched.
        # Restore the product-rule cross term and remove the LO-only doubling.
        body=body.replace('g1_deriv_scale = 0.5d0','g1_deriv_scale = 1.d0')
        if key.startswith('lapw2'):
            for n in [1,2]:
                old=f'radial_buf{n} = radial_buf{n} + 2.d0* &'
                assert body.count(old)==1
                body=body.replace(old,f'radial_buf{n} = radial_buf{n} + &')
        else:
            old='ohalf*g1_deriv_scale*gaunt1*2.d0* &'
            assert body.count(old)==1
            body=body.replace(old,'ohalf*g1_deriv_scale*gaunt1* &')
            old='2.d0*conjg( clm(index1,num,jlo) )*clm(index2,num,jlop)*&'
            assert body.count(old)==1
            body=body.replace(old,'conjg( clm(index1,num,jlo) )*clm(index2,num,jlop)*&')
    (folder/'actual_g1.f90').write_text(stub+body+setup)
    command=['/opt/homebrew/bin/gfortran','-O0','-g','-fcheck=all','-ffree-form','-ffree-line-length-none','-shared','-fPIC',str(folder/'actual_g1.f90'),*[str(DFT/f) for f in ['c_alpha_m.f','gaunt.f','t3j.f','t3j0.f','notri.f','ph.f']],'-o',str(folder/'g1.dylib')]
    subprocess.run(command,cwd=folder,check=True)
    libs[key]=ct.CDLL(str(folder/'g1.dylib'))

# Real primitive radial functions u_lp(r)=r^l(1+a_p*r^2+b_p*r^4).
# Their Cartesian products with real solid harmonics are polynomials, so
# both the radial integrals and Hessian oracle have no numerical derivative.
a=np.array([.17,-.13,.31,-.21]);b=np.array([.03,.07,-.02,.05]);R=1.3
rad=np.zeros((4,L+1,4,L+1,3),order='F')
for i in range(4):
 for l in range(L+1):
  for j in range(4):
   for lp in range(L+1):
    for ai,p in zip([1,a[i],b[i]],[l,l+2,l+4]):
     for aj,q in zip([1,a[j],b[j]],[lp,lp+2,lp+4]):
      base=ai*aj*R**(p+q+1)/(p+q+1)
      rad[i,l,j,lp]+=base*np.array([q*(q-1),q,1])
for lib in libs.values():lib.fill_radial_(rad.ctypes.data_as(ct.c_void_p))

x,y,z=sp.symbols('x y z',real=True);r2=x*x+y*y+z*z
# Choose harmonics with distinct m and complex mixing. Definitions agree
# with Condon--Shortley normalized complex spherical harmonics.
solid=[(0,0,1/sp.sqrt(4*sp.pi)),
       (1,0,sp.sqrt(3/(4*sp.pi))*z),
       (1,1,-sp.sqrt(3/(8*sp.pi))*(x+sp.I*y)),
       (2,0,sp.sqrt(5/(16*sp.pi))*(2*z*z-x*x-y*y)),
       (2,1,-sp.sqrt(15/(8*sp.pi))*z*(x+sp.I*y)),
       (2,2,sp.sqrt(15/(32*sp.pi))*(x+sp.I*y)**2),
       (3,0,sp.sqrt(7/(16*sp.pi))*z*(2*z*z-3*x*x-3*y*y))]
states=[(l,m,p,sp.expand(s*(1+float(a[p])*r2+float(b[p])*r2*r2))) for l,m,s in solid for p in range(4)]
N=len(states)
qr,wr=np.polynomial.legendre.leggauss(14);rr=(qr+1)*R/2;wr=wr*R/2*rr**2
qt,wt=np.polynomial.legendre.leggauss(18);ph=np.arange(40)*2*np.pi/40
rrr,tt,pp=np.meshgrid(rr,qt,ph,indexing='ij')
coords=[(rrr*np.sqrt(1-tt*tt)*np.cos(pp)).ravel(),(rrr*np.sqrt(1-tt*tt)*np.sin(pp)).ravel(),(rrr*tt).ravel()]
w=(wr[:,None,None]*wt[None,:,None]*np.ones((1,1,len(ph)))*2*np.pi/len(ph)).ravel()
def values(expr):return np.broadcast_to(np.asarray(sp.lambdify((x,y,z),expr,'numpy')(*coords),complex),w.shape)
f=np.array([values(s[3]) for s in states]).T
hess=np.empty((3,3,N,N),complex)
for aa in range(3):
 for bb in range(3):
  d=np.array([values(sp.diff(s[3],(x,y,z)[aa],(x,y,z)[bb])) for s in states]).T
  hess[aa,bb]=f.conj().T@(w[:,None]*d)

def code(which,c,legacy=False):
    if legacy:os.environ['WIEN_STRESS_G1_LEGACY']='1'
    else:os.environ.pop('WIEN_STRESS_G1_LEGACY',None)
    alm=np.zeros((NH,1),complex,order='F');blm=alm.copy(order='F');clm=np.zeros((NH,1,2),complex,order='F')
    for coeff,(l,m,p,_) in zip(c,states):
        lm=l*(l+1)+m
        if p==0:alm[lm,0]+=coeff
        elif p==1:blm[lm,0]+=coeff
        else:clm[lm,0,p-2]+=coeff
    ans=np.zeros((3,3,1),complex,order='F');parts=np.zeros((3,3,1,4),complex,order='F')
    args=[np.array([1],np.int32),np.array([1],np.int32),alm,blm,clm,ans]
    if which.startswith('lapw2'):args+=[parts]
    args+=[np.array([R])]
    libs[which].val_kinstrmt_(*[v.ctypes.data_as(ct.c_void_p) for v in args])
    return ans[:,:,0],parts[:,:,0,:]
def direct(c):return np.einsum('i,abij,j->ab',c.conj(),hess,c)
rng=np.random.default_rng(7382)
tests={}
for name,allowed in [('u_only',[0]),('u_udot',[0,1]),('LO_only',[2,3]),('all_primitives',[0,1,2,3])]:
    c=(rng.normal(size=N)+1j*rng.normal(size=N))*.13
    c=np.array([v if s[2] in allowed else 0 for v,s in zip(c,states)])
    oracle=direct(c)
    tests[name]={}
    for key in libs:
        for legacy in [False,True]:
            result,parts=code(key,c,legacy)
            tests[name][key+('_legacy' if legacy else '_default')]=dict(
                error=float(np.max(abs(result-oracle))),real_error=float(np.max(abs((result-oracle).real))),
                trace_error=float(abs(np.trace(result-oracle))),
                diagonal=result.diagonal().real.tolist(),reference_diagonal=oracle.diagonal().real.tolist())
    tests[name]['lapw2_dmft2_difference']=float(np.max(abs(code('lapw2',c)[0]-code('dmft2',c)[0])))

# A full complex density matrix, with arbitrary phases/off-diagonal entries:
# contractions by signed pseudo-bands versus direct matrix trace.
mix=rng.normal(size=(N,4))+1j*rng.normal(size=(N,4));mix*=.1
rho=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));rho=(rho+rho.conj().T)/2
e,v=np.linalg.eigh(rho);pseudo=mix@v
target=np.einsum('mn,im,abij,jn->ab',rho.T,mix.conj(),hess,mix)
tot=sum(e[i]*direct(pseudo[:,i]) for i in range(4))
tests['pseudo_band_identity_error']=float(np.max(abs(tot-target)))
assert tests['pseudo_band_identity_error']<1e-10
assert tests['u_only']['lapw2_legacy']['error']<1e-10
assert tests['u_udot']['dmft2_legacy']['error']<1e-10
assert tests['u_only']['lapw2_default']['real_error']>1e-3
assert tests['LO_only']['lapw2_legacy']['real_error']>1e-3
for name in ['u_only','u_udot','LO_only','all_primitives']:
    for key in ['lapw2_candidate_default','dmft2_candidate_default']:
        assert tests[name][key]['error']<1e-10,(name,key,tests[name][key])

# The simplest counterexample is psi=r Y_10=sqrt(3/(4*pi))*z.
# Its Cartesian Hessian vanishes identically, with no boundary assumptions.
plain=np.zeros_like(rad,order='F')
for l in range(L+1):
 for lp in range(L+1):
    plain[:,l,:,lp,:]=R**(l+lp+1)/(l+lp+1)*np.array([lp*(lp-1),lp,1])[None,None,:]
for lib in libs.values():lib.fill_radial_(plain.ctypes.data_as(ct.c_void_p))
c=np.zeros(N,complex);c[next(i for i,s in enumerate(states) if s[:3]==(1,0,0))]=1
tests['linear_cartesian_function']={key:code(key,c)[0].real.tolist() for key in libs}
tests['linear_cartesian_function']['expected']='zero Hessian'
tests['linear_cartesian_function']['default_error_prediction']=[R**3/15,R**3/15,-2*R**3/15]
(ROOT/'g1_results.json').write_text(json.dumps(tests,indent=2)+'\n')
(ROOT/'source_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [DFT/'val_kinstrmt.F',DMFT/'stress_kernels.f90',DFT/'kin_str_integrals.F',DMFT/'stress_integrals.f90']},indent=2)+'\n')
print(json.dumps(tests,indent=2))
