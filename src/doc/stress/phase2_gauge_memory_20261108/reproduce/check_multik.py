"""Three physical native k points, ragged LAPW bases, global X and mu."""
from pathlib import Path
import sys,time,json,resource,os
import numpy as np
root=Path(__file__).resolve().parent
output=Path(os.environ.get('EDMFT_CHECK_OUTPUT',root));output.mkdir(exist_ok=True)
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from wien_response_density import Density
from lapw_charge_response import NativeChargeNormalization
data_root=Path(os.environ.get('EDMFT_MATERIAL_CHECK_ROOT','/private/tmp/edmft_full_k_response_20261108'))
case=sys.argv[1];source=data_root/case/'charge/0'
old=Path('/private/tmp/edmft_xc_core_response_20261108')/case

def rel(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30))
with load_native_valence_state(source,case,complex_vectors=(case=='SrVO3'),energy_policy='fixed',matrix_quadrature='wien2k24',library='/private/tmp/edmft_phase2_gauge_memory_20261108/build/liblapw_response.so') as f:
    result=dict(case=case,scope='native unmixed charge map on '+str(len(f.records))+' k points; not screened stress',
        k_weights=f.model.kw.tolist(),retained=f.model.eps.shape[1],basis_sizes=[b.retained.shape[0] for b in f.bases],
        preparation_seconds=f.preparation_seconds,matrix_quadrature=f.matrix_quadrature,matrix_overlap_error=f.matrix_overlap_error,central_density_relative=rel(f.spatial.value(),f.reference.vector()),
        complement=[b.solver.info() for b in f.bases],fingerprints=f.fingerprints)
    print(case,'prepared',result['basis_sizes'],'central error',result['central_density_relative'],flush=True)
    state=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
    select,coordinates=state.hartree.native_potential_output(tuple(f.bands.spherical.slices),
        {a:op.layout.harmonics for a,(op,ks) in f.bands.ns.items()},f.bands.ir[0].representatives,
        complex_interstitial=(case=='SrVO3'))
    assert coordinates.nonspherical_shapes==f.bands.input.nonspherical_shapes
    assert coordinates.size==f.bands.input.size
    assert tuple(f.reference.rad)==state.hartree.keys
    core0=Density(old/'hxc_core/0'/(case+'.clmcor'),source/(case+'.struct'))
    raw=f.reference.copy() if hasattr(f.reference,'copy') else f.reference.from_vector(f.reference.vector())
    raw.add(core0)
    ncore=sum(f.structure.mult[a]*sum(s.occ for s in op.states) for a,op in state.core.atoms.items())
    target=float(np.trace(f.model.central_density(),axis1=1,axis2=2).sum().real)+ncore
    norm=state.normalization(raw,target)
    d=np.load('/private/tmp/edmft_hartree_periodic_20261108/'+case+'_direction.npz')['density']
    dv=select.matvec(state.potential_operator.matvec(d))
    result['potential_transpose']=f.bands.operator.check_transpose(samples=1,tolerance=3e-10)
    result['spatial_transpose']=f.spatial.operator.check_transpose(samples=1,tolerance=3e-10)
    result['branches']={}
    from material_gauge import gauge_check,gauge_cutoff_check
    result['gauge']=gauge_check(f,state,select)
    print(case,'gauge',result['gauge'],flush=True)
    if len(f.records)<=3:
        result['gauge_cutoff']=gauge_cutoff_check(f,state,select,source)
        print(case,'gauge cutoff',result['gauge_cutoff'],flush=True)
    result['normalization_alpha']=norm.alpha
    result['target_charge']=target
    # h is the input-density direction amplitude, NOT geometric strain.
    for fixed in [False,True]:
        valence=f.valence_operator(fixed_number=fixed)
        combined=(valence@select+state.core.operator)@state.potential_operator
        normalized=norm.operator@combined
        started=time.perf_counter();y=combined.matvec(d);forward=time.perf_counter()-started
        rr,cross,info=f.complement.response(*f.bands.output.unpack(f.bands.action(dv)),fixed_number=fixed)
        row=dict(metric_band_charge_derivative=info['dN'],chemical_potential_derivative_Ry=info['dmu'],forward_seconds=forward,transpose=normalized.check_transpose(samples=1,tolerance=3e-10),
            normalized_charge_derivative=float((norm.mt+norm.ir)@norm.action(y)),steps={})
        branch='fixed_n' if fixed else 'charge'
        for h in [.001,.0005]:
            if not all((data_root/case/branch/str(sign*h)/'done').exists() for sign in [1,-1]):continue
            ends=[];normalized_ends=[]
            for sign in [1,-1]:
                folder=data_root/case/branch/str(sign*h)
                end=Density(folder/(case+'.clmval'),folder/(case+'.struct'))
                end.add(Density(old/'hxc_core'/str(sign*h)/(case+'.clmcor'),folder/(case+'.struct')))
                ends.append(end.vector())
                normalized_ends.append(NativeChargeNormalization(end,f.structure,target,inversion=state.metadata['inversion']).central)
            fd=(ends[0]-ends[1])/(2*h);nfd=(normalized_ends[0]-normalized_ends[1])/(2*h)
            row['steps'][h]=dict(raw_relative=rel(y,fd),normalized_relative=rel(norm.action(y),nfd),maximum=float(np.max(abs(y-fd))))
        result['branches']['fixed_N' if fixed else 'fixed_mu']=row
        np.savez(output/(case+('_fixed_n' if fixed else '_fixed_mu')+'_compact.npz'),density=y)
        print(case,'fixed_N' if fixed else 'fixed_mu',row,flush=True)
    result['triangles'] = dict(compact_bytes=sum(op.triangle.nbytes if op.lapw else 16 for a,l,k,op in f.bands.spherical.spherical),
        former_full_bytes=sum(op.triangle.shape[0]**2*16 for a,l,k,op in f.bands.spherical.spherical))
    result['peak_process_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (output/(case+'_multik.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(case,'report saved',flush=True)
