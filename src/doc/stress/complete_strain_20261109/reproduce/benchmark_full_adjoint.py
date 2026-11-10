"""Warm ABBA timing of complete native J^T with identical prepared inputs."""
from pathlib import Path
import sys,time,json,importlib.util
import numpy as np
root=Path(__file__).resolve().parent;sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from lapw_klein_functional import PreparedKleinFunctional
from wien_response_density import Density
import lapw_spherical_response as spherical
import lapw_nonspherical_response as nonspherical
import lapw_compact_blocks as blocks
def old(name):
 spec=importlib.util.spec_from_file_location('original_'+name,'/Users/haule/dbin/'+name+'.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
old_blocks=old('lapw_compact_blocks');old_spherical=old('lapw_spherical_response');old_spherical.right_multiply_columns=old_blocks.right_multiply_columns
case='FeSe';source=root/case/'refined';pc=load_potential_core_state(source/'potential_core.npz');core=Density(source/(case+'.clmcor'),source/(case+'.struct'))
start=time.perf_counter();rows=[];values=[];fast=spherical.PreparedSphericalBand.pullback
with load_native_valence_state(source,case,complex_vectors=False,energy_policy='fixed',matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so') as v:
 state=NativeScreenedState(v,pc,core,target_charge=120,fixed_number=True);obs=PreparedKleinFunctional(state)
 for label in ['warm','original','optimized','optimized','original']:
  if label=='original':
   spherical.PreparedSphericalBand.pullback=old_spherical.PreparedSphericalBand.pullback
   nonspherical.right_multiply_columns=old_blocks.right_multiply_columns
  else:
   spherical.PreparedSphericalBand.pullback=fast
   nonspherical.right_multiply_columns=blocks.right_multiply_columns
  t=time.perf_counter();value=state.jacobian.rmatvec(obs.charge_gradient);dt=time.perf_counter()-t
  rows.append(dict(mode=label,seconds=dt));values.append(value);print(rows[-1],flush=True)
 errors=[float(np.linalg.norm(v-values[0])/np.linalg.norm(values[0])) for v in values]
 if max(errors)>2e-12:raise ValueError(errors)
 means={name:float(np.mean([v['seconds'] for v in rows if v['mode']==name])) for name in ['original','optimized']}
 report=dict(case=case,threads=1,scope='complete J transpose at identical fast FFT grid; includes spherical product reuse and flattened RR/QR products',rows=rows,mean_seconds=means,speedup=means['original']/means['optimized'],relative_errors=errors,total_seconds=time.perf_counter()-start)
 (root/case/'adjoint_benchmark.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
