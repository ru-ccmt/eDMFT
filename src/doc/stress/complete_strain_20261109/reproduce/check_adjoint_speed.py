"""One complete-functional adjoint, paired-map speed and forward agreement."""
from pathlib import Path
import sys,time,json,resource,importlib.util
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from lapw_klein_functional import PreparedKleinFunctional
from wien_response_density import Density
from stress_response_adjoint import PairedOperator
from lapw_spherical_response import PreparedSphericalBand
spec=importlib.util.spec_from_file_location('original_spherical','/Users/haule/dbin/lapw_spherical_response.py');original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
case=sys.argv[1];source=root/case/'current';prior=root/case/'current_readout';out=root/case/'adjoint';out.mkdir(exist_ok=True)
t0=time.perf_counter()
pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
core=Density('/private/tmp/edmft_xc_core_response_20261108/'+case+'/hxc_core/0/'+case+'.clmcor',source/(case+'.struct'))
parts=np.load(prior/'source.npz');solution=np.load(prior/'solution.npz')
report=dict(case=case,threads=1,timings={})
with load_native_valence_state(source,case,complex_vectors=case=='SrVO3',energy_policy='fixed',matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so') as v:
 report['timings']['prepare_valence_seconds']=v.preparation_seconds
 state=NativeScreenedState(v,pc,core,target_charge={'FeSe':120,'SrVO3':85}[case],fixed_number=True)
 t=time.perf_counter();obs=PreparedKleinFunctional(state);report['timings']['prepare_readout_seconds']=time.perf_counter()-t
 print(case,'prepared',report,flush=True)
 fast=PreparedSphericalBand.pullback;outputs=[]
 try:
  for label,fn in [('original',original.PreparedSphericalBand.pullback),('optimized',fast)]:
   PreparedSphericalBand.pullback=fn;t=time.perf_counter();outputs.append(state.jacobian.rmatvec(obs.charge_gradient));dt=time.perf_counter()-t
   report['timings'][label+'_JT_seconds']=dt;print(case,label,'JT',dt,flush=True)
 finally:PreparedSphericalBand.pullback=fast
 err=float(np.linalg.norm(outputs[0]-outputs[1])/max(np.linalg.norm(outputs[0]),1e-30));report['optimized_relative_error']=err
 if err>3e-12:raise ValueError('optimized transpose differs')
 report['JT_speedup']=report['timings']['original_JT_seconds']/report['timings']['optimized_JT_seconds']
 op=state.screened.jacobian;times=[]
 def transpose(x):
  t=time.perf_counter();y=op.rmatvec(x);times.append(time.perf_counter()-t);print(case,'JT',len(times),times[-1],flush=True);return y
 state.screened.jacobian=PairedOperator(op.input_size,op.output_size,op.matvec,transpose,'timed optimized native J')
 d,lam,r=obs.adjoint_contractions(parts['band_geometry'][:,None],parts['potential_source'][:,None],np.eye(3)[None],parts['source'][:,None],rtol=2e-9,atol=1e-11,restart=30,max_cycles=8)
 forward=float(json.loads((prior/'REPORT.json').read_text())['derivative_Ry'])
 report.update(derivative_Ry=float(d[0]),forward_Ry=forward,forward_adjoint_error_Ry=float(d[0]-forward),mean_tensile_stress_GPa=float(d[0]*14710.5077/(3*pc.hartree.volume)),readout=r)
 from lapw_energy_readout import screened_energy_gradient
 t=time.perf_counter();ge=screened_energy_gradient(state,obs,lam)
 report['energy_gradients']={f'{a},{l}':g.tolist() for (a,l),g in ge.items()}
 report['timings']['energy_covector_seconds']=time.perf_counter()-t
 report['timings'].update(total_seconds=time.perf_counter()-t0,solve_action_seconds=times)
 report['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 np.savez_compressed(out/'covectors.npz',adjoint=lam,charge_gradient=obs.charge_gradient)
 (out/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2),flush=True)
