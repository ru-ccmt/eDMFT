"""Differentiate one finite functional containing every projector chain.

Calls the existing production Fortran routines, compiled separately with
bounds checks.  Endpoint eigensystems and inverse square roots are independent.
No production sources are changed.
"""
import ctypes as ct
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.special import expit

root = Path(__file__).parent
lib = ct.CDLL(str(root/'harness/audit_routines.dylib'))
F = np.asfortranarray
I = lambda a: np.asfortranarray(a, dtype=np.int32)
Z = lambda s: np.zeros(s, complex, order='F')
def call(name, *args):
    getattr(lib, name+'_')(*[a.ctypes.data_as(ct.c_void_p) for a in args])
def herm(a): return (a+a.conj().T)/2
def invsqrt(a):
    d,v=eigh(a)
    return (v/np.sqrt(d))@v.conj().T
rng=np.random.default_rng(71704)
def rand(*s): return rng.normal(size=s)+1j*rng.normal(size=s)
n,nc,nk=7,3,3
lo,hi=1,5
nw=hi-lo
W=slice(lo,hi)
weights=np.array([.2,.3,.5])
ee=np.array([-2.1,-1.2,-.6,.2,.8,1.7,2.5])
eps=[ee+.03*k for k in range(nk)]
p=[rand(n,nc)*.4 for k in range(nk)]
dh=[herm(rand(n,n))*.13 for k in range(nk)]
do=[herm(rand(n,n))*.045 for k in range(nk)]
dq=[rand(n,nc)*.06 for k in range(nk)]
sig=herm(rand(nc,nc))*.35
bath_coupling=rand(nc,nc)*.28
bath_energies=np.array([-.8,.1,1.3])
dynamic=False
temp,mu=.24,-.27
rot=np.diag([1,-1,1]).astype(complex)
def sym(a, enabled): return .5*(a+rot.conj().T@a@rot) if enabled else a
maps=[I([v]) for v in [1,n,lo+1,hi,n,nw,nc,1,nc,1]]
tail=[I([1]),I([nc]),I([nc]),I(np.arange(1,nc+1)[:,None])]
si=I(np.ones((nc,nc,1)))

def endpoint(t, sym_on):
    pp=[];en=[]
    for k in range(nk):
        e,v=eigh(np.diag(eps[k])+t*dh[k],np.eye(n)+t*do[k])
        en.append(e)
        pp.append(v[:,W].conj().T@(p[k]+t*dq[k]))
    s=sym(sum(w*a.conj().T@a for w,a in zip(weights,pp)),sym_on)
    x=invsqrt(s)
    free=0.
    for w,e,pw in zip(weights,en,pp):
        u=pw@x
        eff=np.diag(e).astype(complex)
        eff[W,W]+=u@sig@u.conj().T
        if dynamic:
            ub=np.zeros((n,nc),complex);ub[W]=u@bath_coupling
            eff=np.block([[eff,ub],[ub.conj().T,np.diag(bath_energies)]])
        ev=eigh(eff,eigvals_only=True)
        free-=w*temp*np.logaddexp(0,-(ev-mu)/temp).sum()
        if dynamic:
            free+=w*temp*np.logaddexp(0,-(bath_energies-mu)/temp).sum()
    return free

def analytic(sym_on, sym_adjoint):
    s=sym(sum(w*a[W].conj().T@a[W] for w,a in zip(weights,p)),sym_on)
    x=invsqrt(s)
    ds_exp=sym(sum(w*(b[W].conj().T@a[W]+a[W].conj().T@b[W])
                  for w,a,b in zip(weights,p,dq)),sym_on)
    dxrel=Z((nc,nc,1,1))
    call('computenormalizationderivative',dxrel,F(ds_exp[:,:,None,None]),
         F(x[:,:,None]),si,I([nc]),I([nc]),I([1]),I([0]),I([1]))
    dx=x@dxrel[:,:,0,0]
    moment=np.zeros((nc,nc),complex)
    pieces=dict(oneparticle=0.,explicit=0.,norm_explicit=0.,window=0.,norm_window=0.)
    for k,w in enumerate(weights):
        uall=p[k]@x;u=uall[W]
        eff=np.diag(eps[k]).astype(complex)
        eff[W,W]+=u@sig@u.conj().T
        if dynamic:
            ub=np.zeros((n,nc),complex);ub[W]=u@bath_coupling
            eff=np.block([[eff,ub],[ub.conj().T,np.diag(bath_energies)]])
        ev,v=eigh(eff)
        rhofull=w*(v*expit((mu-ev)/temp))@v.conj().T
        rho=rhofull[:n,:n]
        rhoe=(eff@rhofull)[:n,:n]
        c=sig@u.conj().T@rho[W,W]
        if dynamic:c+=bath_coupling@rhofull[n:,W]
        moment+=u.T@c.T
        pieces['oneparticle']+=np.trace(rho@dh[k]-rhoe@do[k]).real
        pieces['explicit']+=2*np.trace(c@dq[k][W]@x).real
        pieces['norm_explicit']+=2*np.trace(c@p[k][W]@dx).real
        rr,er,grad=Z((n,n)),Z((n,n)),Z((nw,nc,1))
        call('buildwindowresponsematrices',rr,er,grad,F(uall[:,:,None]),
             F(rho[W,W]),F(rhoe[W,W]),F(eps[k]),*maps,*tail,si,
             F(c.T[:,:,None]),I([1]))
        pieces['window']+=np.trace(rr@dh[k]-er@do[k]).real
    if sym_adjoint: moment=sym(moment.conj(),True).conj()
    qb=Z((nc,nc,1))
    call('buildnormalizationadjoint',qb,F(moment[:,:,None]),F(x[:,:,None]),
         si,I([nc]),I([nc]),I([1]))
    for k,w in enumerate(weights):
        rr,er,grad=Z((n,n)),Z((n,n)),Z((nw,nc,1))
        call('buildnormalizationwindowresponsematrices',rr,er,grad,
             F(p[k][:,:,None]),qb,F([w]),F(eps[k]),*maps,*tail)
        pieces['norm_window']+=np.trace(rr@dh[k]-er@do[k]).real
    return dict(pieces=pieces,total=float(sum(pieces.values())))

results={}
for dynamic in [False,True]:
 for sym_on in [False,True]:
    fd={str(h):(endpoint(h,sym_on)-endpoint(-h,sym_on))/(2*h)
        for h in [1e-4,1e-5,1e-6]}
    a=analytic(sym_on,False)
    b=analytic(sym_on,True) if sym_on else a
    results[f'dynamic={dynamic},sym={sym_on}']=dict(fd=fd,production_default=a,with_symmetry_adjoint=b)
    assert abs(b['total']-fd['1e-05'])<2e-9
print(json.dumps(results,indent=2))
(root/'joint_projector_results.json').write_text(json.dumps(results,indent=2)+'\n')

# The symmetry-adjoint failure also exists when Sigma itself respects the
# orbital group. It is not an artifact of a symmetry-breaking self-energy.
saved_sig,saved_bath=sig.copy(),bath_coupling.copy()
sig=sym(sig,True);bath_coupling=sym(bath_coupling,True)
h=1e-5;fd=(endpoint(h,True)-endpoint(-h,True))/(2*h)
a=analytic(True,False);b=analytic(True,True)
extra={'finite_difference':fd,'default_chain':a['total'],'adjoint_complete':b['total']}
assert abs(fd-b['total'])<2e-9
(root/'symmetry_preserving_sigma.json').write_text(json.dumps(extra,indent=2)+'\n')
sig,bath_coupling=saved_sig,saved_bath
