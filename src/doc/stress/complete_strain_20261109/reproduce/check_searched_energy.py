"""Complete frozen-Sigma derivative with analytic native radial-energy searches."""
from pathlib import Path
import sys,time,json,resource
import numpy as np
root=Path(__file__).resolve().parent;sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from lapw_klein_functional import PreparedKleinFunctional
from lapw_energy_response import SearchedEnergyResponse
from wien_response_density import Density
from stress_response_adjoint import PairedOperator
case=sys.argv[1];source=root/case/'current';prior=root/case/'current_readout';out=root/case/'searched_readout';out.mkdir(exist_ok=True)
start=time.perf_counter();pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
core=Density('/private/tmp/edmft_xc_core_response_20261108/'+case+'/hxc_core/0/'+case+'.clmcor',source/(case+'.struct'))
b=np.load(prior/'source.npz');fixed=np.load(prior/'solution.npz');fixed_report=json.loads((prior/'REPORT.json').read_text())
report=dict(case=case,threads=1,timings={})
with load_native_valence_state(source,case,complex_vectors=case=='SrVO3',energy_policy='fixed',matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so') as v:
 report['timings']['prepare_valence_seconds']=v.preparation_seconds
 state=NativeScreenedState(v,pc,core,target_charge={'FeSe':120,'SrVO3':85}[case],fixed_number=True)
 t=time.perf_counter();obs=PreparedKleinFunctional(state);report['timings']['prepare_functional_seconds']=time.perf_counter()-t
 print(case,'prepared',report,flush=True)
 t=time.perf_counter();response=SearchedEnergyResponse(state,obs,source,case,complex_vectors=case=='SrVO3',library='/Users/haule/dbin/liblapw_response.so',progress=lambda i,item:print(case,'energy column',i,item,flush=True))
 report['timings']['prepare_energy_policy_seconds']=time.perf_counter()-t
 report['root_checks']=response.details;report['direct_energy_gradient_Ry_per_Ry']=response.g.tolist()
 adjpath=root/case/'adjoint/REPORT.json'
 if adjpath.exists():
  ar=json.loads(adjpath.read_text());lam=np.load(root/case/'adjoint/covectors.npz')['adjoint']
  effective=response.g+lam@response.P
  reverse=np.array([sum(np.array(ar['energy_gradients'][f"{d['atom']},{d['l']}"])[d['slots']]) for d in response.details])
  report['energy_gradient_forward_reverse_error']=float(np.max(abs(effective-reverse)))
  report['effective_energy_gradient_Ry_per_Ry']=effective.tolist()
  if report['energy_gradient_forward_reverse_error']>1e-9:raise ValueError('energy forward/reverse differ')
  fd=json.loads((root/case/'energy_parameters/REPORT.json').read_text())['derivative']
  native=np.array([fd[f"{d['atom']},{d['l']}"][d['slots'][0]] for d in response.details])
  report['native_fd_energy_policy_correction_Ry']=float(effective@native)
  report['native_fd_energy_policy_mean_GPa']=float((fixed_report['derivative_Ry']+effective@native)*14710.5077/(3*pc.hartree.volume))
  estimate=response.total_energy_directions(b['potential_source'],fixed['density'])
  report['first_pass_energy_policy_mean_GPa']=float((fixed_report['derivative_Ry']+effective@estimate)*14710.5077/(3*pc.hartree.volume))
 s=response.complete_source(b['source'],b['potential_source']);J=response.jacobian;times=[]
 def action(x):
  t=time.perf_counter();y=J.matvec(x);times.append(time.perf_counter()-t);print(case,'coupled J',len(times),times[-1],flush=True);return y
 response.screened.jacobian=PairedOperator(J.input_size,J.output_size,action,J.rmatvec,'timed searched-energy Jacobian')
 q,solver=response.screened.solve(s,initial=fixed['density'],rtol=1e-10,atol=1e-11,restart=30,max_cycles=8)
 de=response.total_energy_directions(b['potential_source'],q)
 D,row=obs.forward_contractions(b['band_geometry'][:,None],b['potential_source'][:,None],np.eye(3)[None],q[:,None])
 D=float(D[0]+response.g@de);report.update(readout=row,screened=solver,derivative_Ry=D,mean_tensile_stress_GPa=D*14710.5077/(3*pc.hartree.volume),searched_energy_derivatives_Ry=de.tolist(),fixed_energy_mean_GPa=fixed_report['mean_tensile_stress_GPa'])
 report['timings'].update(J_actions=times,total_seconds=time.perf_counter()-start)
 report['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 np.savez_compressed(out/'solution.npz',density=q,A=response.A,P=response.P,g=response.g,charge_gradient=response.charge_gradient,source=s)
 (out/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
