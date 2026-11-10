"""Complete symmetry-allowed tensor, one adjoint, current central inputs."""
from pathlib import Path
import sys,time,json,resource,importlib.util
import numpy as np
root=Path(__file__).resolve().parent;sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from lapw_total_stress import PreparedTotalStress,PreparedStrainSource
from stress_response_adjoint import PairedOperator
from wien_response_density import Density
case=sys.argv[1];source=root/case/'current';out=root/case/'tensor';out.mkdir(exist_ok=True)
start=time.perf_counter();pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
core=Density('/private/tmp/edmft_xc_core_response_20261108/'+case+'/hxc_core/0/'+case+'.clmcor',source/(case+'.struct'))
timing={}
with load_native_valence_state(source,case,complex_vectors=case=='SrVO3',energy_policy='fixed',matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so') as v:
 timing['prepare_valence_seconds']=v.preparation_seconds
 state=NativeScreenedState(v,pc,core,target_charge={'FeSe':120,'SrVO3':85}[case],fixed_number=True)
 t=time.perf_counter();total=PreparedTotalStress(state,source,case,complex_vectors=case=='SrVO3',library='/Users/haule/dbin/liblapw_response.so',progress=lambda i,d:print(case,'prepared energy',i,d,flush=True));timing['prepare_functional_and_searches_seconds']=time.perf_counter()-t
 print(case,'prepared',timing,'strains',total.strains.tolist(),flush=True)
 sources=[];timing['source_seconds']=[]
 for i,e in enumerate(total.strains):
  t=time.perf_counter();saved=out/f'source{i}.npz'
  if i==0:
   b=np.load(root/case/'current_readout/source.npz');x=PreparedStrainSource(e,total.energy.complete_source(b['source'],b['potential_source']),b['band_geometry'],b['potential_source'],dict(reused_hydrostatic_source=True))
  elif saved.exists():
   b=np.load(saved);x=PreparedStrainSource(e,b['source'],b['geometry'],b['potential'],dict(reused_same_snapshot=True))
  else:x=total.prepare_source(e)
  np.savez_compressed(saved,strain=x.strain,source=x.source,geometry=x.geometry,potential=x.potential)
  timing['source_seconds'].append(time.perf_counter()-t);sources.append(x);print(case,'source',i,timing['source_seconds'][-1],flush=True)
 times=[];J=total.energy.jacobian
 def jt(x):
  t=time.perf_counter();y=J.rmatvec(x);times.append(time.perf_counter()-t);print(case,'total JT',len(times),times[-1],flush=True);return y
 total.energy.screened.jacobian=PairedOperator(J.input_size,J.output_size,J.matvec,jt,'full searched-energy transpose')
 warm=root/case/'adjoint/covectors.npz';initial=np.load(warm)['adjoint'] if warm.exists() else None
 sigma,lam,report=total.contract(sources,initial=initial)
 report['timings']=timing|dict(JT_actions=times,total_seconds=time.perf_counter()-start)
 report['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 np.savez_compressed(out/'solution.npz',adjoint=lam,charge_gradient=total.energy.charge_gradient,A=total.energy.A,P=total.energy.P,g=total.energy.g)
 old=root/case/'searched_readout/REPORT.json'
 if old.exists():report['forward_adjoint_mean_error_GPa']=report['mean_tensile_stress_GPa']-json.loads(old.read_text())['mean_tensile_stress_GPa']
 (out/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
