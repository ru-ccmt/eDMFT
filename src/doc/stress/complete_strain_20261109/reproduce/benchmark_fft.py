"""Same physical first-k density map on minimal and FFT-friendly grids."""
from pathlib import Path
import sys,time,json,importlib.util
import numpy as np
from scipy.fft import next_fast_len
root=Path(__file__).resolve().parent;sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from frozen_sigma_tangent import CachedFrozenSigmaResponse
from lapw_density_response import FourierStars,InterstitialDensity
from wstruct import Struct
from wien_response_density import Density
spec=importlib.util.spec_from_file_location('old_density','/Users/haule/dbin/lapw_density_response.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
for case in sys.argv[1:] or ['FeSe','SrVO3']:
 d=root/case/'current';s=Struct();s.ReadStruct(str(d/(case+'.struct')));rho=Density(d/(case+'.clmval'),d/(case+'.struct'))
 model,records=CachedFrozenSigmaResponse.from_export(d/'response');r=records[0];n=r.nbasis-r.nlo;c=r.coefficients[:n,:,0];g=r.labels[:,:n].T
 allg=[r.labels[:,:r.nbasis-r.nlo].T for r in records];extent=np.max([np.ptp(g,axis=0) for g in allg],axis=0);reps=np.array(list(rho.pw));star_extent=np.max(abs(np.einsum('sij,qj->sqi',s.timat.transpose(0,2,1),reps)),axis=(0,1));minimal=tuple(2*np.maximum(extent,star_extent)+1);fast=tuple(next_fast_len(int(n)) for n in minimal)
 variants={}
 for label,cls,shape in [('old_minimal',old.InterstitialDensity,minimal),('band_first_minimal',InterstitialDensity,minimal),('band_first_fast',InterstitialDensity,fast)]:
  stars=FourierStars(shape,reps,s.timat.transpose(0,2,1),s.tau,complex_phases=case=='SrVO3');variants[label]=cls(c,g,s.a*s.b*s.c,stars,cache_wavefunctions=False)
 y=np.random.default_rng(18).normal(size=len(reps))+1j*np.random.default_rng(19).normal(size=len(reps));rhom=np.eye(c.shape[1])*.07;output={};timing={k:[] for k in variants}
 # Alternate order after warm-up, so a first-call cache cannot fake a speedup.
 for order in [list(variants),list(variants)[::-1],list(variants),list(variants)[::-1]]:
  for name in order:
   op=variants[name];t=time.perf_counter();v=op.value(rhom);a,b=op.pullback(y,rhom);elapsed=time.perf_counter()-t;timing[name].append(elapsed);output[name]=np.r_[v.ravel(),a.ravel(),b.ravel()]
 baseline=output['old_minimal'];errors={k:float(np.linalg.norm(v-baseline)/np.linalg.norm(baseline)) for k,v in output.items()}
 if max(errors.values())>1e-12:raise ValueError(errors)
 report=dict(case=case,minimal=list(map(int,minimal)),fast=list(map(int,fast)),bands=c.shape[1],timings=timing,median_seconds={k:float(np.median(v[1:])) for k,v in timing.items()},relative_errors=errors)
 report['speedup']=report['median_seconds']['old_minimal']/report['median_seconds']['band_first_fast']
 (root/case/'fft_benchmark.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
