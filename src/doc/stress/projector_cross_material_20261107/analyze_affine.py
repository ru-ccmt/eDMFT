import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1'
from pathlib import Path
import numpy as np,json,sys
from scipy.linalg import solve
r=Path(__file__).parent/'SrVO3_affine';prior=Path('/private/tmp/projector_actual_stress_20261107');c=Path('/private/tmp/projector_stationarity_probe_20261107/SrVO3')
ns={'__file__':str(prior/'analyze.py')};exec((prior/'analyze.py').read_text().split("names=['center'")[0],ns)
vectors=ns['vectors'];read_dump=ns['read_dump'];evidence=ns['evidence']
gauge=np.load(prior/'dftks_wien_gauge.npy');adv=np.load(prior/'advection_half_potential_derivative.npy')
# Reuse independently verified potential-star expansion and step convolution.
fns={'__file__':str(prior/'interstitial.py')};exec((prior/'interstitial.py').read_text().split("v,mask=read_stars")[0],fns)
stars=fns['read_stars'];theta=fns['theta'];product=fns['product'];N=fns['N'];th=theta(0)
names=['center','p001','p0005','p00025'];iters={n:vectors(r/(n+'_basis')/'SrVO3.vector_x') for n in names};full=vectors(prior/'complete_center'/'SrVO3.vector_x')
meta=json.loads((r/'timings.json').read_text());vc,mask=stars(r/'center'/'SrVO3.vtotal');wc,_=stars(r/'center'/'SrVO3.vns')
ops={}
for name in names[1:]:
 h=meta[name]['actual_step'];vp,_=stars(r/name/'SrVO3.vtotal');wp,_=stars(r/name/'SrVO3.vns')
 it=product((vp-vc)/h,th);it[~mask]=0;ops[name]=(it,((wp-wc)/h).conj())
matrices={n:{k:[] for k in ['raw','interstitial_raw','interstitial_true']} for n in names[1:]}
for k in range(1,76):
 d=read_dump(c/f'dmft_sternheimer.k{k:08d}.bin');e0,v0,l0,rad0=next(full);inv=solve(v0,np.eye(len(v0)),check_finite=False)[:43];H={}
 for name in names:
  e,v,l,rad=next(iters[name]);assert np.array_equal(l,l0);t=inv@v;H[name]=gauge[k-1]@((t*e)@t.conj().T)@gauge[k-1].conj().T
 npw=d.nbasis-d.nlo;labels=d.labels[:,:npw];diff=(labels[:,:,None]-labels[:,None,:])%N;A=d.coefficients[:npw,:,0]
 for name in names[1:]:
  h=meta[name]['actual_step'];matrices[name]['raw'].append((H[name]-H['center'])/h)
  for key,op in zip(['interstitial_true','interstitial_raw'],ops[name]):matrices[name][key].append(A.conj().T@op[tuple(diff)]@A)
result={};fac=14710.507641078/(3*382.70904)
for name in names[1:]:
 m={k:np.array(v) for k,v in matrices[name].items()};m['MT_physical']=m['raw']-m['interstitial_raw'];m['MT_only']=m['MT_physical']+adv;m['total']=m['MT_only']+m['interstitial_true'];res={}
 for key in ['MT_physical','interstitial_true','MT_only','total']:
  val={b:float(np.einsum('kij,kji->',evidence[b],m[key]).real) for b in ['response','norm_response']};val['GPa']=(val['response']+val['norm_response'])*fac;res[key]=val
 res['corrected_analytic_GPa']=5.67688415+res['total']['GPa'];result[name]=res
 np.savez_compressed(r/(name+'_matrices.npz'),**m)
print(json.dumps(result,indent=2),flush=True);(r/'results.json').write_text(json.dumps(result,indent=2))
