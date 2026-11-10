"""Check the assembled file-coordinate operator against the native density probes."""
from pathlib import Path
import sys,runpy,time,json
import numpy as np
ROOT=Path(__file__).resolve().parent
case=sys.argv[1]
# Native endpoint runs are reused, not regenerated. The earlier driver also
# verifies the complete radial pullback on those same inputs.
d=runpy.run_path(str(ROOT/'check_moving_density.py'))
from lapw_density_response import MuffinTinDensitySum,SpatialValenceDensity,InterstitialDensity
mat,s=d['mat'],d['s'];rho,drho=d['rho'],d['drho']
atoms={};matching={};radial={};site=0
start=time.perf_counter()
for ia in range(s.nat):
    hs=[key[1:] for key in d['reference'].rad if key[0]==ia]
    phi={l:mat.context[ia,l]['primitive'][:,:2] for l in range(9)}
    radial[ia]={l:mat.context[ia,l]['dprimitive'][:,:2] for l in range(9)}
    entries=[];dots=[]
    for mu in range(s.mult[ia]):
        for k in range(len(rho)):
            zs=[];dzs=[]
            for l in range(9):
                ctx=mat.context[ia,l];ch=mat.channels[ia,site,k,l]
                zs.append(ch.value(ctx['primitive'][-1],ctx['lo'],lapw=ctx['lapw']).reshape(-1,rho.shape[1]))
                dzs.append(ch.action(ctx['primitive'][-1],ctx['dprimitive'][-1],ctx['dlo'],lapw=ctx['lapw']).reshape(-1,rho.shape[1]))
            entries.append((k,np.concatenate(zs),1/s.mult[ia]));dots.append(np.concatenate(dzs))
        site+=1
    atoms[ia]=MuffinTinDensitySum(entries,rho,phi,hs,d['small']);matching[ia]=dots
irs=[]
for k,r in enumerate(d['records']):
    pw=r.nbasis-r.nlo
    irs.append(InterstitialDensity(r.coefficients[:pw,:,0],d['allg'][k],s.a*s.b*s.c,d['stars']))
full=SpatialValenceDensity(d['reference'],atoms,irs)
setup=time.perf_counter()-start
value=full.value();actual=full.action(drho,matching,radial)
np.testing.assert_allclose(actual,d['joint'].vector(),atol=2e-12,rtol=1e-10)
np.testing.assert_allclose(value,d['reference'].vector(),atol=3e-12,rtol=1e-8)
y=d['y'].vector();rb,zb,pb=full.pullback(y)
reverse=np.vdot(rb,drho).real
for ia in atoms:
    reverse+=sum(np.vdot(a,b).real for a,b in zip(zb[ia],matching[ia]))
    reverse+=sum(np.vdot(pb[ia][l],radial[ia][l]).real for l in radial[ia])
forward=y@actual
op,coords=full.fixed_basis_operator()
np.testing.assert_allclose(coords.unpack(op.rmatvec(y)),rb,atol=2e-12,rtol=1e-10)
np.testing.assert_allclose(op.matvec(coords.pack(drho)),d['fixed'].vector(),atol=2e-12,rtol=1e-10)

def timing(call):
    times=[]
    for i in range(5):
        start=time.perf_counter();call();times.append(time.perf_counter()-start)
    return float(np.median(times))
times=dict(setup_seconds=setup,
    fixed_basis_action_seconds=timing(lambda:full.action(drho)),
    moving_basis_action_seconds=timing(lambda:full.action(drho,matching,radial)),
    fixed_basis_reverse_seconds=timing(lambda:op.rmatvec(y)),
    complete_reverse_seconds=timing(lambda:full.pullback(y)))
report=dict(material=case,kpoints=len(rho),bands=rho.shape[1],
    central_native_relative=float(np.linalg.norm(value-d['reference'].vector())/np.linalg.norm(d['reference'].vector())),
    joint_reassembly_maximum=float(np.max(abs(actual-d['joint'].vector()))),
    reverse_relative=float(abs(forward-reverse)/max(1.,abs(forward),abs(reverse))),
    timing=times,scope='3-k native density reconstruction, not complete screening or stress')
(ROOT/case/'spatial_adapter.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
