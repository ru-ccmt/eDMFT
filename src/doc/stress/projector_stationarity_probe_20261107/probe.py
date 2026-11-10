"""Local-potential directional residual in the central stored KS space.

The perturbation is +t Ry in the interstitial, zero inside all MT spheres.
Geometry, radial basis, physical local orbitals, effective self-energy and
chemical potential are fixed. This computes the projector part of dGamma/dt,
not a stress, and not a self-consistent density or full finite-Gamma test.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from pathlib import Path
import json,time,hashlib
import numpy as np
from scipy.linalg import eigh
from check_projector_sternheimer_fd import read_dump,read_normalization,read_response_matrix
from check_interstitial_overlap_stress import read_atoms,theta

root=Path(__file__).parent
case=root/'SrVO3'
lines=(case/'SrVO3.struct').read_text().splitlines()
cell=[float(lines[3][i:i+10]) for i in range(0,60,10)]
assert np.allclose(cell[3:],[90,90,90]) and lines[1].startswith('P')
br=np.diag(2*np.pi/np.array(cell[:3]));vol=np.prod(cell[:3])
atoms=read_atoms(case/'SrVO3.struct')

def response(p,g,e,window,outside):
    rr=np.zeros((len(e),len(e)),complex)
    for j,n in enumerate(window):
        d=e[n]-e[outside]
        assert np.min(np.abs(d))>1e-8
        value=np.einsum('aij,ij->a',p[outside],g[j])/d
        rr[outside,n]=value
        rr[n,outside]=value.conj()
    return rr

totals={'fixed_X':0.,'normalization':0.,'uniform_shift':0.}
steps=[1e-3,3e-4,1e-4,3e-5]
fd={str(h):0. for h in steps}
max_matrix_difference=0.
min_window_gap=np.inf
min_interstitial_eigenvalue=np.inf
max_interstitial_eigenvalue=-np.inf
points=[]
matrix_evidence=[]
start=time.monotonic()
paths=sorted(case.glob('dmft_sternheimer.k????????.bin'))
assert len(paths)==75,len(paths)
for path in paths:
    d=read_dump(path)
    norm=read_normalization(case/f'dmft_sternheimer_norm.k{d.kpoint:08d}.bin')
    saved=read_response_matrix(case/f'stress_response.k{d.kpoint:08d}.bin')
    window=np.arange(d.dm_nemin-d.nemin,d.dm_nemax-d.nemin+1)
    outside=np.setdiff1d(np.arange(d.nstored),window)
    min_window_gap=min(min_window_gap,np.min(np.abs(d.energies[window,None]-d.energies[outside])))
    npw=d.nbasis-d.nlo
    oi=theta(d.labels[:,:npw],br,atoms,vol,0.)
    dv=np.zeros((d.nstored,d.nstored),complex)
    for spinor in range(d.ncomponents):
        a=d.coefficients[:npw,:,spinor]
        dv+=a.conj().T@oi@a
    assert np.max(np.abs(dv-dv.conj().T))<1e-12
    ev=eigh(dv,eigvals_only=True)
    min_interstitial_eigenvalue=min(min_interstitial_eigenvalue,float(ev[0]))
    max_interstitial_eigenvalue=max(max_interstitial_eigenvalue,float(ev[-1]))
    rr=response(d.projector,d.gradient,d.energies,window,outside)
    rn=response(norm.projector,norm.gradient,d.energies,window,outside)
    max_matrix_difference=max(max_matrix_difference,float(np.max(np.abs(rr-saved.density))))
    sx=float(np.trace(rr@dv).real)
    sn=float(np.trace(rn@dv).real)
    totals['fixed_X']+=sx;totals['normalization']+=sn
    totals['uniform_shift']+=float(np.trace(rr+rn).real)
    # Independent endpoint eigensystems evaluate the spectral projection
    # Pi_W(t)q in the original fixed Hilbert space, avoiding phase alignment
    # or arbitrary internal rotations within degenerate window subspaces.
    # Contraction uses the central exact functional gradients; therefore this
    # checks the window variation, not a full nonlinear trace-log endpoint.
    pk=d.projector.reshape(d.nstored,-1)
    pn=norm.projector.reshape(d.nstored,-1)
    for h in steps:
        vals=[]
        for sign in [1,-1]:
            e,v=eigh(np.diag(d.energies)+sign*h*dv)
            v=v[:,window]
            up=(v@(v.conj().T@pk)).reshape(d.projector.shape)
            pp=(v@(v.conj().T@pn)).reshape(norm.projector.shape)
            vals.append(2*float(np.real(np.sum(d.gradient*up[window])+np.sum(norm.gradient*pp[window]))))
        fd[str(h)]+=(vals[0]-vals[1])/(2*h)
    points.append({'k':int(d.kpoint),'fixed_X':sx,'normalization':sn})
    matrix_evidence.append({'energies':d.energies,'potential':dv,
        'projector':d.projector,'gradient':d.gradient,
        'raw_projector':norm.projector,'norm_gradient':norm.gradient,
        'response':rr,'norm_response':rn})
    if d.kpoint%15==0:print('k',d.kpoint,'elapsed',round(time.monotonic()-start,2),flush=True)
total=totals['fixed_X']+totals['normalization']
result={'scope':__doc__,'material':'SrVO3','kpoints':len(paths),'nstored':43,
 'window':[12,31],'perturbation':'delta V=t Ry times interstitial characteristic function',
 'projector_residual_Ry_per_unit_t':total,'components':totals,
 'spectral_projector_fd_contracted_with_central_gradients':fd,
 'fd_errors':{s:x-total for s,x in fd.items()},
 'not_a_full_trace_log_fd':True,'fixed_X_dump_matrix_max_difference':max_matrix_difference,
 'minimum_window_boundary_gap_Ry':float(min_window_gap),
 'interstitial_matrix_eigenvalue_range':[min_interstitial_eigenvalue,max_interstitial_eigenvalue],
 'per_k':points,'seconds':time.monotonic()-start,
 'installed_executable_sha256':hashlib.sha256(Path('/Users/haule/dbin/dmft2').read_bytes()).hexdigest()}
(root/'results.json').write_text(json.dumps(result,indent=2)+'\n')
np.savez_compressed(root/'band_matrix_evidence.npz',
    **{name:np.stack([point[name] for point in matrix_evidence])
       for name in matrix_evidence[0]})
print(json.dumps({k:v for k,v in result.items() if k!='per_k'},indent=2))
