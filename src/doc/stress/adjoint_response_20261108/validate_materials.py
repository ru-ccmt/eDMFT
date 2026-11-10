import json, time, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent/'src'))
from frozen_sigma_response import FrozenSigmaResponse, herm
from frozen_sigma_tangent import CachedFrozenSigmaResponse, real_inner
from stress_response_adjoint import PairedOperator, ScreenedResponse, electronic_operator

results={}
for case in ['SrVO3','FeSe']:
    prefix='/private/tmp/edmft_linear_response_20261108/'+case+'/central/response'
    old,_=FrozenSigmaResponse.from_export(prefix)
    new,_=CachedFrozenSigmaResponse.from_export(prefix, cache_mb=1024)
    rng=np.random.default_rng(112)
    shape=(*new.eps.shape,new.eps.shape[1])
    x=herm(rng.normal(size=shape)+1j*rng.normal(size=shape))*0.02
    y=herm(rng.normal(size=shape)+1j*rng.normal(size=shape))*0.02
    start=time.monotonic();a,ai=old.response(x);old_time=time.monotonic()-start
    start=time.monotonic();b,bi=new.response(x);cold=time.monotonic()-start
    start=time.monotonic();c,ci=new.response(x);warm=time.monotonic()-start
    start=time.monotonic();rev=new.adjoint(y);adjtime=time.monotonic()-start
    left,right=real_inner(y,c),real_inner(rev,x)
    grad=new.projector_potential_gradient()
    # Test the six-source adjoint using this material's full electronic map
    # and a SMALL SYNTHETIC spatial projection. These sources are algebraic
    # controls, not six physical strains and not new material stress values.
    chi,coords=electronic_operator(new)
    basis=rng.normal(size=(coords.size,8))
    basis=np.linalg.qr(basis,mode='reduced')[0]*5
    potential=PairedOperator(8,coords.size,lambda v:basis@v,lambda v:basis.T@v,'synthetic potential')
    density=PairedOperator(coords.size,8,lambda v:basis.T@v,lambda v:basis@v,'synthetic reconstruction')
    j=density@chi@potential
    csmall=potential.rmatvec(coords.pack(grad))
    sources=rng.normal(size=(8,6))
    solver=ScreenedResponse(j)
    start=time.monotonic();acon,lam,arep=solver.contractions(csmall,sources,rtol=1e-11,atol=1e-14);adj6=time.monotonic()-start
    start=time.monotonic();forward=np.column_stack([solver.solve(sources[:,i],rtol=1e-11,atol=1e-13)[0] for i in range(6)]);fwd6=time.monotonic()-start
    sixerror=float(np.max(abs(acon-csmall@forward)))
    row={'shape':shape, 'window':len(new.W), 'frequencies':len(new.meta.omega),
         'forward_max_abs_error':float(np.max(abs(a-b))),
         'forward_relative_error':float(np.linalg.norm(a-b)/np.linalg.norm(a)),
         'adjoint_left':left,'adjoint_right':right,
         'adjoint_relative_error':abs(left-right)/max(abs(left),abs(right),1e-30),
         'reference_seconds':old_time,'cold_seconds':cold,'warm_seconds':warm,
         'adjoint_seconds':adjtime,'kernel_speedup':old_time/warm,
         'cache':new.cache_info(),'fixed_N_dN':bi['dN'],
         'common_projector_gradient_norm':float(np.linalg.norm(grad)),
         'six_source_algebraic_control':{'max_contraction_error':sixerror,'six_forward_seconds':fwd6,'one_adjoint_seconds':adj6,'adjoint_report':arep,'physical_strain_sources':False},
         'scope':'saved central basis and stored bands; no radial/core/spatial-potential response'}
    results[case]=row
    np.save(Path(__file__).parent/(case+'_projector_gradient.npy'),grad)
    Path(__file__).with_name('material_validation.json').write_text(json.dumps(results,indent=2)+'\n')
    print(case,json.dumps(row),flush=True)
    assert row['forward_relative_error']<1e-11
    assert row['adjoint_relative_error']<1e-10
    assert sixerror<1e-11
