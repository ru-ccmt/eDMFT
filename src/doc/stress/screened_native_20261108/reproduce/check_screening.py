"""Physical native hydrostatic source through the complete analytic Jacobian."""
from pathlib import Path
import sys,time,json,resource
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from wien_response_density import Density
from stress_response_adjoint import PairedOperator,ScreenedResponse
from native_reference import NativeReference
case=sys.argv[1];mesh=sys.argv[2] if len(sys.argv)>2 else 'three_k'
native=NativeReference(case,mesh);source=native.source;work=native.work
with load_native_valence_state(source,case,complex_vectors=(case=='SrVO3'),energy_policy='fixed',
                              matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so') as valence:
    print(case,mesh,'valence prepared',valence.preparation_seconds,flush=True)
    pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
    core=Density(native.xcroot/case/'hxc_core/0'/(case+'.clmcor'),source/(case+'.struct'))
    state=NativeScreenedState(valence,pc,core,target_charge=native.target,fixed_number=True)
    print(case,'screening prepared',json.dumps(state.metadata),flush=True)
    calls={'forward':0,'adjoint':0};j=state.jacobian
    def counted(x,adjoint=False):
        key='adjoint' if adjoint else 'forward';start=time.monotonic()
        y=j.rmatvec(x) if adjoint else j.matvec(x);calls[key]+=1
        print(case,key,calls[key],'seconds',time.monotonic()-start,flush=True)
        return y
    paired=PairedOperator(j.input_size,j.output_size,counted,lambda y:counted(y,True),'native charge screening')
    state.screened=ScreenedResponse(paired)
    b=np.load(work/'source_0.000125.npy')
    result=state.solve_source(b,rtol=2e-8,atol=1e-9,restart=30,max_cycles=8)
    q=result.density;np.save(work/'analytic_density.npy',q)
    np.save(work/'analytic_induced_potential.npy',result.induced_potential)
    report=dict(forward=result.report,preparation_seconds=valence.preparation_seconds,source_norm=float(np.linalg.norm(b)),response_norm=float(np.linalg.norm(q)),native_residuals={})
    (work/'screening_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(case,'FORWARD SOLVED',json.dumps(report['forward']),flush=True)
    # Independent native program actions on the solved direction, two step sizes.
    # These are NOT used by GMRES and include new projector normalization and mu.
    jq=j.matvec(q)
    for h in [1e-5,5e-6]:
        p,_=native.evaluate('verify_p_'+str(h),native.rho.from_vector(native.rho.vector()+h*q))
        m,_=native.evaluate('verify_m_'+str(h),native.rho.from_vector(native.rho.vector()-h*q))
        fd=(p-m)/(2*h);residual=b+fd-q
        report['native_residuals'][h]=dict(norm=float(np.linalg.norm(residual)),relative_to_source=float(np.linalg.norm(residual)/np.linalg.norm(b)),jacobian_relative=float(np.linalg.norm(fd-jq)/max(np.linalg.norm(jq),1e-30)))
        print(case,'NATIVE RESIDUAL',h,report['native_residuals'][h],flush=True)
    # One physical observable: integrated MT charge of the correlated species.
    # Its adjoint is NOT the Klein/stress gradient. This tests the actual large
    # nonsymmetric native transpose solver without making an energy claim.
    c=np.zeros_like(b);sl=pc.core.slices.get(0)
    if sl is None: raise ValueError('expected atom zero spherical channel')
    c[sl]=state.normalization.mt[sl]
    if mesh=='three_k':
        lam,adj=state.screened.solve(c,adjoint=True,rtol=2e-8,atol=1e-11,restart=30,max_cycles=8)
        direct=float(c@q);reverse=float(lam@b)
        report['adjoint']=dict(adj,observable='MT electron number on atom 1 including its multiplicity',forward_contraction=direct,adjoint_contraction=reverse,absolute_difference=abs(direct-reverse))
        np.save(work/'mt_charge_adjoint.npy',lam)
        print(case,'ADJOINT SOLVED',json.dumps(report['adjoint']),flush=True)
        # A second native source uses the same scalar-observable adjoint.
        b2=np.load(work/'source_0.00025.npy')
        q2,rep2=state.screened.solve(b2,initial=q,rtol=2e-8,atol=1e-9,restart=30,max_cycles=8)
        report['second_source']=dict(rep2,response_step_relative=float(np.linalg.norm(q2-q)/np.linalg.norm(q)),forward_contraction=float(c@q2),adjoint_contraction=float(lam@b2))
    report['calls']=calls;report['peak_process_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (work/'screening_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(case,'completed',flush=True)
