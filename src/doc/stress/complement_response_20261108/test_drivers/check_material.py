"""Compare projected complement with full spectral/native material reference."""
from pathlib import Path
import sys,time,json
import numpy as np
ROOT=Path(__file__).resolve().parent
OLD=Path('/private/tmp/edmft_coupled_valence_response_20261108')
sys.path[:0]=[str(OLD),'/Users/haule/projects/eDMFT_new/src/python']
from build_material import build
sys.path.insert(0,str(ROOT/'stage'))
from lapw_complement_response import LAPWComplementBasis,EmptyComplementResponse
from frozen_sigma_tangent import CachedFrozenSigmaResponse,real_inner
from check_oneparticle_stress_fd import read_basis_matrix

case=sys.argv[1];nr={'FeSe':66,'SrVO3':43}[case]
label={'FeSe':'FeSe_b509_one','SrVO3':'SrVO3_b703_one'}[case]
source=OLD/label/'central'
fixture=build(case,source);full=fixture['model'];records=fixture['records']
saved=np.load(OLD/(label+'_coupled.npz'));h,o,p=[saved[k] for k in ['h','o','p']]

def native(local):
    meta=full.meta
    result=np.zeros((*local.shape[:2],int(max(meta.nindo)),len(meta.cix)),complex)
    for orb,c in enumerate(meta.cix):
        if c>=0:
            for j in range(meta.nindo[orb]):result[:,:,j,orb]=local[:,:,meta.isx[j,orb],c]
    return result

start=time.perf_counter()
small=CachedFrozenSigmaResponse(full.meta,full.eps[:,:nr],native(full.raw[:,:nr]),
                  native(full.apply_x(full.raw[:,:nr],full.X)),full.kw,full.W)
bases=[];rotations=[];rawq=[];hqr=[];oqr=[];checks=[]
for k,record in enumerate(records):
    hh=read_basis_matrix(source/f'native.ham.k{k+1:08d}.bin',b'EDMFTHV1')
    oo=read_basis_matrix(source/f'native.k{k+1:08d}.bin',b'EDMFTOV1')
    np.testing.assert_array_equal(hh.labels,record.labels)
    c=record.coefficients[:,:,0]
    basis=LAPWComplementBasis.prepare(hh.matrix,oo.matrix,c[:,:nr],full.eps[k,:nr],
                  mu=full.meta.mu,temperature=full.meta.temperature)
    # Full spectral states are used here ONLY to compare gauges with the old
    # independent oracle. The solver above received H/O and retained C alone.
    rotation=c[:,nr:].conj().T@oo.matrix@basis.complement
    np.testing.assert_allclose(rotation.conj().T@rotation,np.eye(len(rotation)),atol=1e-9)
    rotations.append(rotation);bases.append(basis)
    rawq.append((rotation.conj().T@full.raw[k,nr:].reshape(len(rotation),-1)).reshape(-1,*full.raw.shape[2:]))
    hqr.append(rotation.conj().T@h[k,nr:,:nr]);oqr.append(rotation.conj().T@o[k,nr:,:nr])
    checks.append(basis.checks|basis.solver.info())
model=EmptyComplementResponse(small,[b.solver for b in bases],rawq)
prepare_seconds=time.perf_counter()-start
def forward():return model.response(h[:,:nr,:nr],o[:,:nr,:nr],p[:,:nr],hqr,oqr,fixed_number=False)
rr,qr,info=forward();actual=np.zeros_like(h);actual[:,:nr,:nr]=rr
for k,u in enumerate(rotations):actual[k,nr:,:nr]=u@qr[k];actual[k,:nr,nr:]=actual[k,nr:,:nr].conj().T
expected,_=full.generalized_response(h,o,p,fixed_number=False)
v=fixture['bands'].pack({a:.07*cs[0].radius*np.sin(1.3*cs[0].radius) for a,cs in fixture['contexts'].items()})
physical=fixture['spatial'].action(actual)+fixture['direct'].action(v)

rng=np.random.default_rng(159);y=rng.normal(size=rr.shape)+1j*rng.normal(size=rr.shape);y=(y+y.swapaxes(-1,-2).conj())/2
yq=[rng.normal(size=x.shape)+1j*rng.normal(size=x.shape) for x in qr]
def reverse():return model.adjoint(y,yq,fixed_number=False)
hb,ob,rb,hq,oq=reverse()
left=real_inner(y,rr)+2*sum(real_inner(a,b) for a,b in zip(yq,qr))
right=(real_inner(hb,h[:,:nr,:nr])+real_inner(ob,o[:,:nr,:nr])+real_inner(rb,p[:,:nr])+
       2*sum(real_inner(a,b) for a,b in zip(hq,hqr))+2*sum(real_inner(a,b) for a,b in zip(oq,oqr)))
residuals=[]
for k,b in enumerate(bases):
    rhs=hqr[k]-oqr[k]*b.solver.energies
    residuals.append(b.solver.residual(b.solver.action(rhs),rhs))
def timing(operation):
    times=[]
    for _ in range(5):
        start=time.perf_counter();operation();times.append(time.perf_counter()-start)
    return float(np.median(times))
full_y=np.zeros_like(h);full_y[:,:nr,:nr]=y
for k,u in enumerate(rotations):full_y[k,nr:,:nr]=u@yq[k];full_y[k,:nr,nr:]=full_y[k,nr:,:nr].conj().T
report=dict(case=case,retained_bands=nr,checks=checks,preparation_seconds=prepare_seconds,
     coefficient_response_relative=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected)),
     spatial_response_relative_to_full=float(np.linalg.norm(physical-saved['analytic'])/np.linalg.norm(saved['analytic'])),
     spatial_response_relative_to_native_fd=float(np.linalg.norm(physical-saved['last_fd'])/np.linalg.norm(saved['last_fd'])),
     transpose_error=float(abs(left-right)/max(1,abs(left),abs(right))),solve_relative_residuals=residuals,
     forward_seconds=timing(forward),reverse_seconds=timing(reverse),
     full_forward_seconds=timing(lambda:full.generalized_response(h,o,p,fixed_number=False)),
     full_reverse_seconds=timing(lambda:full.generalized_adjoint(full_y,fixed_number=False)),
     full_cache=full.cache_info(),retained_cache=small.cache_info())
(ROOT/(case+'_complement.json')).write_text(json.dumps(report,indent=2)+'\n')
np.savez(ROOT/(case+'_complement.npz'),coefficient=actual,spatial=physical)
print(json.dumps(report,indent=2),flush=True)
for cs in fixture['contexts'].values():
    for c in cs.values():c.close()
