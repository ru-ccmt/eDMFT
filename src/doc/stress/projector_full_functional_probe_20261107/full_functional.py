"""Independent trace-log endpoints for fixed and moving physical embedding.

Uses SrVO3 central band matrices exported by production dmft2. Frozen
self-energy, chemical potential, basis, radial orbitals and geometry. No SCF.
All internal energies are Ry. Weights include spin degeneracy exactly once.
The trial potential is t Ry times the interstitial characteristic function.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.special import expit
import json,time
import re

root=Path(__file__).parent
prior=Path('/private/tmp/projector_stationarity_probe_20261107')
case=root/'inputs' if (root/'inputs').is_dir() else prior/'SrVO3'
evidence=root/'band_matrix_evidence.npz'
if not evidence.exists():evidence=prior/'band_matrix_evidence.npz'
d=np.load(evidence)
eps=d['energies'];q=d['raw_projector'][:,:,:3,0]
ufull=d['projector'][:,:,:3,0];pot=d['potential']
nk,nb=eps.shape
W=np.arange(11,31)
kw=np.array([float(s.split()[5]) for s in (case/'SrVO3.klist').read_text().splitlines()
             if s.strip() and not s.startswith('END')]);kw/=kw.sum()
RY=13.60569253
dat=np.loadtxt(case/'sig.inp1')/RY
om=dat[:-1,0];T=om[0]/np.pi
sigma=dat[:-1,1::2]+1j*dat[:-1,2::2]
sinf=dat[-1,1::2]
mu=float((case/'EF.dat').read_text())/RY
fw=np.load(root/'weights.npy')
z=1j*om
assert abs(fw.sum()-2387)<1e-8
# The production moment selector reported GOOD for all three components.
# Refit its two-moment model independently using numpy least squares.
fitx=1/z[-24:]
coef=np.linalg.lstsq(np.column_stack([np.ones(24),fitx]),
    z[-24:,None]*(sigma[-24:]-sinf),rcond=None)[0]
s1=coef[0].real

def sym(s):
    # P4mm active t2g representation: xz/yz form a two-dimensional irrep,
    # xy is a separate irrep. This equals the full forward group average.
    av=(s[0,0].real+s[1,1].real)/2
    return np.diag([av,av,s[2,2].real])

def invsqrt(s):
    e,v=eigh(s)
    assert e.min()>0
    return (v/np.sqrt(e))@v.conj().T

X0=np.linalg.lstsq(q.reshape(-1,3),ufull.reshape(-1,3),rcond=None)[0]
X0=(X0+X0.conj().T)/2
Sraw0=np.einsum('k,kna,knb->ab',kw,q[:,W].conj(),q[:,W])
Ssaved=np.linalg.inv(X0)@np.linalg.inv(X0)
Soffset=Ssaved-sym(Sraw0)
u0=np.zeros_like(q);u0[:,W]=ufull[:,W]

def embedding(t,mode,operator=pot):
    hs=np.array([np.diag(e) for e in eps],complex)+t*operator
    if mode=='fixed':return hs,u0
    raw=[]
    for k in range(nk):
        _,v=eigh(hs[k]);v=v[:,W]
        raw.append(v@(v.conj().T@q[k]))
    raw=np.array(raw)
    if mode=='moving_fixed_X':return hs,raw@X0
    s=sym(np.einsum('k,kna,knb->ab',kw,raw.conj(),raw))
    # Production uses a saved central X. Anchor the endpoint function to
    # that same X; the tiny offset is reported and a fresh-X control follows.
    if mode=='moving':s=s+Soffset
    elif mode!='moving_fresh_X':raise ValueError(mode)
    return hs,raw@invsqrt(s)

def one_fsp(h,u):
    hinf=h+(u*sinf)@u.conj().T
    e,v=eigh(hinf);xi=e-mu
    a=v.conj().T@u
    resolvent=1/(z[:,None]-xi)
    gl=np.einsum('ia,wi,ib->wab',a.conj(),resolvent,a,optimize=True)
    ratio=np.eye(3)[None,:,:]-gl*(sigma-sinf)[:,None,:]
    ld=np.linalg.slogdet(ratio)[1]
    static=-T*np.logaddexp(0,-xi/T).sum()
    dynamic=-2*T*np.dot(fw,ld)
    b=np.sum(np.abs(a)**2*s1[None,:],axis=1)
    c0=1.
    divided=(expit(-c0/T)-expit(-xi/T))/(c0-xi)
    near=np.abs(c0-xi)<1e-8
    divided[near]=-expit(-c0/T)*expit(c0/T)/T
    exact=np.dot(b,divided)
    grid=2*T*np.dot(fw,np.real((resolvent@b)/(z-c0)))
    return np.array([static,dynamic,exact-grid])

def fsp(h,u):
    return sum(2*w*one_fsp(hk,uk) for w,hk,uk in zip(kw,h,u))

if (root/'central_density.npy').exists():
    rho=np.load(root/'central_density.npy')
else:
    from check_oneparticle_stress_fd import read_one_particle
    rho=[]
    for k in range(nk):
        rr=read_one_particle(case/f'dmft_rho_rhoe.k{k+1:08d}.bin')
        r=np.zeros((nb,nb),complex);r[:rr.nfull,:rr.nfull]=rr.rho
        # Upper retained bands outside the DMFT window have negligible ordinary
        # Fermi occupation but include it for the finite-space functional.
        for a in range(rr.nfull,nb):r[a,a]=2*kw[k]*expit((mu-eps[k,a])/T)
        rho.append(r)
    rho=np.array(rho)
    np.save(root/'central_density.npy',rho)
bare=float(np.einsum('kij,kji->',rho,pot).real)
N=float(np.trace(rho,axis1=1,axis2=2).real.sum())
pred=float(np.einsum('kij,kji->',d['response']+d['norm_response'],pot).real)
start=time.monotonic()
base=fsp(*embedding(0,'fixed'))
out={'central_fsp_components_Ry':base.tolist(),'bare_density_slope_Ry':bare,
 'density_particle_number':N,'predicted_projector_slope_Ry':pred,
 'central_normalization_relative_residual':float(np.linalg.norm(X0@sym(Sraw0)@X0-np.eye(3))),
 'normalization_offset_norm':float(np.linalg.norm(Soffset)),
 'tail_moment1_Ry2':s1.tolist(),'mu_Ry':mu,'temperature_Ry':T,'steps':{}}
production_vsum=float(re.search(r'^:VSUM\s*:.*?=\s*([-+\d.Ee]+)',
    (case/'SrVO3.scf2').read_text(),re.M).group(1))
out['central_VSUM_check']={'production_Ry':production_vsum,
    'independent_Ry':float(base.sum()+mu*41),
    'difference_Ry':float(base.sum()+mu*41-production_vsum)}
assert abs(out['central_VSUM_check']['difference_Ry'])<1e-8
for step in [1e-3,3e-4,1e-4,3e-5]:
    row={}
    for mode in ['fixed','moving_fixed_X','moving']:
        fp=fsp(*embedding(step,mode));fm=fsp(*embedding(-step,mode))
        slope=float(np.sum(fp-fm)/(2*step))
        row[mode]={'fsp_slope':slope,'tangent_functional_slope':slope-bare,
            'component_slopes':((fp-fm)/(2*step)).tolist()}
    row['moving_minus_fixed']=row['moving']['fsp_slope']-row['fixed']['fsp_slope']
    row['projector_prediction_error']=row['moving_minus_fixed']-pred
    out['steps'][str(step)]=row
    print('step',step,json.dumps(row),flush=True)
    (root/'results.json').write_text(json.dumps(out,indent=2)+'\n')

# Fresh overlap rather than anchored saved normalization: central discrepancy
# is only ~1e-8; nevertheless report its effect explicitly.
h=1e-4
fresh=(fsp(*embedding(h,'moving_fresh_X')).sum()-fsp(*embedding(-h,'moving_fresh_X')).sum())/(2*h)
out['fresh_normalization_fsp_slope_Ry']=float(fresh)
identity=np.broadcast_to(np.eye(nb),pot.shape)
up=fsp(*embedding(h,'fixed',identity)).sum()
um=fsp(*embedding(-h,'fixed',identity)).sum()
out['uniform_shift_control']={'fixed_fsp_slope':float((up-um)/(2*h)),
    'tangent_slope':float((up-um)/(2*h)-N),
    'embedding_difference':float(np.max(np.abs(embedding(h,'moving',identity)[1]-u0)))}
uniform=[]
for uh in [1e-4,3e-5]:
    fp=fsp(*embedding(uh,'fixed',identity));fm=fsp(*embedding(-uh,'fixed',identity))
    bp=fsp(*embedding(uh,'moving',identity));bm=fsp(*embedding(-uh,'moving',identity))
    uniform.append({'h':uh,'tangent_slope':float((fp-fm).sum()/(2*uh)-N),
        'tail_slope':float((fp-fm)[2]/(2*uh)),
        'moving_minus_fixed_slope':float(((bp-bm)-(fp-fm)).sum()/(2*uh))})
out['uniform_shift_control']['step_check']=uniform

# Algebraic separation of the same matrix perturbation. These two controls
# are band-matrix diagnostics, not separately asserted to be local potentials.
cross=pot.copy();mask=np.zeros(nb,bool);mask[W]=True
cross[:,mask[:,None]==mask[None,:]]=0
block=pot-cross
out['band_block_controls']={}
for label,op in [('within_blocks_only',block),('cross_window_only',cross)]:
    vals=[]
    for sign in [1,-1]:
        ff=fsp(*embedding(sign*h,'fixed',op)).sum()
        fb=fsp(*embedding(sign*h,'moving',op)).sum()
        vals.append(float(fb-ff))
    out['band_block_controls'][label]={'endpoint_differences':vals,
        'moving_minus_fixed_slope':float((vals[0]-vals[1])/(2*h))}

# If the window spans every retained state, Pi_W is the identity at every
# endpoint. Check this with actual endpoint eigenvectors, without assuming it.
full_s=sym(np.einsum('k,kna,knb->ab',kw,q.conj(),q));full_x=invsqrt(full_s)
full_u=q@full_x
full_ctrl=[]
for t in [h,-h]:
    hh=np.array([np.diag(e) for e in eps],complex)+t*pot
    rr=[]
    for k in range(nk):
        _,v=eigh(hh[k]);rr.append(v@(v.conj().T@q[k]))
    rr=np.array(rr);xx=invsqrt(sym(np.einsum('k,kna,knb->ab',kw,rr.conj(),rr)))
    full_ctrl.append(float(fsp(hh,rr@xx).sum()-fsp(hh,full_u).sum()))
out['complete_stored_space_control']={'endpoint_functional_differences_Ry':full_ctrl,
    'slope_difference_Ry':float((full_ctrl[0]-full_ctrl[1])/(2*h))}
out['seconds']=time.monotonic()-start
(root/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print('FINAL',json.dumps(out,indent=2),flush=True)
