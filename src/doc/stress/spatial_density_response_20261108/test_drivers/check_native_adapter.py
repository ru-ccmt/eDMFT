"""Material check of the shared in-memory radial -> spatial density operator."""
from pathlib import Path
import sys,runpy,json,time
import numpy as np
ROOT=Path(__file__).resolve().parent
d=runpy.run_path(str(ROOT/'check_spatial_adapter.py'))
from lapw_radial_native import NativeRadialContext
from lapw_basis_density_response import BasisDensityResponse
case=d['case'];mat=d['mat'];s=d['s'];full=d['full'];driver=d['d']
contexts={};matching={};site=0
start=time.perf_counter()
try:
 for a in range(s.nat):
  contexts[a]={};matching[a]=[]
  for l in range(9):
   args=list(map(float,(ROOT/case/f'radial_a{a}_l{l}.input').read_text().split()))
   n,nlo,rel,ell,lapw=map(int,args[:5]);r0,dx,rmt,energy,charge,de=args[5:11]
   elo=np.array(args[11:11+nlo]);v=np.array(args[11+2*nlo:11+2*nlo+n])
   ctx=NativeRadialContext(v,r0=r0,dx=dx,rmt=rmt,energy=energy,lo_energies=elo,
          ell=ell,charge=charge,relativistic=bool(rel),lapw=bool(lapw),library=ROOT/'liblapw_response.so')
   contexts[a][l]=ctx;old=mat.context[a,l]
   np.testing.assert_allclose(ctx.primitive,old['primitive'],rtol=1e-13,atol=1e-13)
   np.testing.assert_allclose(ctx.lo_coefficients,old['lo'],rtol=1e-13,atol=1e-13)
   dp,dm=ctx.action(driver['dv'][a])
   np.testing.assert_allclose(dp,old['dprimitive'],rtol=1e-13,atol=1e-13)
   np.testing.assert_allclose(dm,old['dlo'],rtol=1e-13,atol=1e-13)
  for mu in range(s.mult[a]):
   for k in range(len(full.rho)):
    matching[a].append({l:mat.channels[a,site,k,l] for l in range(9)})
   site+=1
 op=BasisDensityResponse(full,contexts,matching,energy_policy='fixed')
 setup=time.perf_counter()-start
 direction=op.pack(driver['dv']);actual=op.action(direction)
 reference=driver['basis'].vector();y=driver['y'].vector()
 np.testing.assert_allclose(actual,reference,rtol=1e-10,atol=1e-12)
 reverse=op.pullback(y);fwd=y@actual;rev=reverse@direction
 def timed(call):
  times=[]
  for i in range(7):
   t=time.perf_counter();call();times.append(time.perf_counter()-t)
  return float(np.median(times))
 report=dict(material=case,channels=sum(len(v) for v in contexts.values()),
   reference_relative=float(np.linalg.norm(actual-reference)/np.linalg.norm(reference)),
   physical_direction_transpose_relative=float(abs(fwd-rev)/max(1.,abs(fwd),abs(rev))),
   random_transpose_relative=op.operator.check_transpose(samples=2,tolerance=1e-11),
   timing=dict(setup_seconds=setup,forward_seconds=timed(lambda:op.action(direction)),
               reverse_seconds=timed(lambda:op.pullback(y))),
   scope='in-memory radial/basis density at 3 k points; fixed E, RMT, C and geometry; not full stress')
 (ROOT/case/'native_adapter.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2),flush=True)
finally:
 for channels in contexts.values():
  for ctx in channels.values():ctx.close()
