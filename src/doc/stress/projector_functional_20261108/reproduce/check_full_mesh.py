"""Prepare physical projector covectors on a complete configured material mesh.

Use the independently validated, already solved hydrostatic density direction.
Do not interpret this isolated induced contribution as a total stress correction.
No new screened solve, native program execution or user-file change is made.
"""
from pathlib import Path
import sys,json,time,resource
import numpy as np
root=Path(__file__).resolve().parent
if root.name=='reproduce':root=root.parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from lapw_projector_functional import PreparedProjectorFunctional
from wien_response_density import Density

case=sys.argv[1]
source=Path('/private/tmp/edmft_full_k_response_20261108/full_bz')/case/'charge/0'
previous=Path('/private/tmp/edmft_screened_native_20261108/full_bz')/case
out=root/'full_bz'/case;out.mkdir(parents=True,exist_ok=True)
with load_native_valence_state(source,case,complex_vectors=(case=='SrVO3'),energy_policy='fixed',
                              matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so') as v:
    print(case,'prepared full mesh',v.preparation_seconds,flush=True)
    pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
    core=Density(Path('/private/tmp/edmft_xc_core_response_20261108')/case/'hxc_core/0'/(case+'.clmcor'),source/(case+'.struct'))
    state=NativeScreenedState(v,pc,core,target_charge={'SrVO3':85,'FeSe':120}[case],fixed_number=True)
    start=time.monotonic();obs=PreparedProjectorFunctional(state)
    elapsed=time.monotonic()-start
    q=np.load(previous/'analytic_density.npy');potential=pc.potential_operator.matvec(q)
    value=obs.potential_contraction(potential);reverse=float(obs.charge_gradient@q)
    report=dict(case=case,kpoints=len(v.records),retained=v.model.eps.shape[1],
                scope='isolated induced-potential projector derivative only; not a stress correction',
                valence_preparation_seconds=v.preparation_seconds,gradient_preparation_seconds=elapsed,
                induced_projector_derivative_Ry=value,charge_pullback_Ry=reverse,
                pullback_error_Ry=abs(value-reverse),
                potential_gradient_norm=float(np.linalg.norm(obs.potential_gradient)),
                charge_gradient_norm=float(np.linalg.norm(obs.charge_gradient)),
                central=state.metadata,source_fingerprints=v.fingerprints,
                peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    np.savez(out/'covectors.npz',potential=obs.potential_gradient,charge=obs.charge_gradient)
    (out/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:x for k,x in report.items() if k!='source_fingerprints'}),flush=True)
