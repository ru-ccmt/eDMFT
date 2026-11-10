"""Complete analytic readout at current native central radial parameters."""
from pathlib import Path
import os,sys,time,json,resource
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from lapw_strain_charge import hydrostatic_charge_source
from lapw_klein_functional import PreparedKleinFunctional
from wien_response_density import Density
from stress_response_adjoint import PairedOperator
case=sys.argv[1];source=root/case/'current';out=root/case/'current_readout';out.mkdir(exist_ok=True)
start=time.perf_counter()
pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
core=Density('/private/tmp/edmft_xc_core_response_20261108/'+case+'/hxc_core/0/'+case+'.clmcor',source/(case+'.struct'))
report=dict(case=case,scope='current central radial basis, fixed-energy hydrostatic tangent',timings={})
with load_native_valence_state(source,case,complex_vectors=case=='SrVO3',energy_policy='fixed',matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so') as v:
 report['timings']['prepare_valence_seconds']=v.preparation_seconds
 state=NativeScreenedState(v,pc,core,target_charge={'FeSe':120,'SrVO3':85}[case],fixed_number=True)
 report['central']=state.metadata;print(case,'prepared',report,flush=True)
 before=time.perf_counter();obs=PreparedKleinFunctional(state)
 report['timings']['prepare_functional_seconds']=time.perf_counter()-before
 before=time.perf_counter();s,parts,source_report=hydrostatic_charge_source(state,source,case)
 report['timings']['source_seconds']=time.perf_counter()-before
 report['source']=source_report
 np.savez_compressed(out/'source.npz',source=s,**parts)
 actions=[];op=state.screened.jacobian
 def action(x):
  t=time.perf_counter();y=op.matvec(x);actions.append(time.perf_counter()-t);print(case,'J',len(actions),actions[-1],flush=True);return y
 state.screened.jacobian=PairedOperator(op.input_size,op.output_size,action,op.rmatvec,'timed current J')
 warm=np.load('/private/tmp/edmft_strain_charge_20261109/full_bz/'+case+'/analytic_density.npy')
 q=state.solve_source(s,initial=warm,rtol=2e-9,atol=1e-10,restart=30,max_cycles=8)
 report['screened']=q.report
 geom=parts['band_geometry'][:,None];pot=parts['potential_source'][:,None];strains=np.eye(3)[None]
 D,readout=obs.forward_contractions(geom,pot,strains,q.density[:,None])
 report['readout']=readout;report['derivative_Ry']=float(D[0]);report['mean_tensile_stress_GPa']=float(D[0]*14710.5077/(3*pc.hartree.volume));report['pressure_GPa']=-report['mean_tensile_stress_GPa']
 report['timings'].update(total_seconds=time.perf_counter()-start,J_actions=actions)
 report['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 np.savez_compressed(out/'solution.npz',density=q.density,potential=q.induced_potential,charge_gradient=obs.charge_gradient,potential_gradient=obs.potential_gradient,lattice_joint=obs.lattice.joint_gradient)
 (out/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2),flush=True)
