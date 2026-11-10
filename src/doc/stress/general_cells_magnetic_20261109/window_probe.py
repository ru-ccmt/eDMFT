from pathlib import Path
import sys,json,time,copy
import numpy as np
sys.path.insert(0,str(Path(__file__).parent/'stage'))
from frozen_sigma_tangent import CachedFrozenSigmaResponse,real_inner
from frozen_sigma_response import herm,invsqrt
from lapw_native_snapshot import read_basis_matrix
from lapw_lattice import NativeLattice
from wstruct import Struct
ROOT=Path(__file__).parent

def aligned(rec,mat):
 npw=rec.nbasis-rec.nlo
 assert mat.nbasis==rec.nbasis and mat.nlo==rec.nlo
 index={tuple(x):i for i,x in enumerate(mat.labels[:,:npw].T)}
 order=[index[tuple(x)] for x in rec.labels[:,:npw].T]+list(range(npw,rec.nbasis))
 return mat.matrix[np.ix_(order,order)]

results={}
for material,ds,c in [('MnO',['MnO_center','MnO_m_0.001','MnO_p_0.001'],'MnO'),('MnTe',['dmft','dmft_m_0.001','dmft_p_0.001'],'0_800_MnTe')]:
 paths=[ROOT/'native'/material/d for d in ds]
 s=Struct();s.ReadStruct(str(paths[0]/(c+'.struct')));vol=NativeLattice.from_structure(s).volume
 conversion=14710.507641078/(3*vol)
 for spin in ['']+(['_dn'] if material=='MnTe' else []):
  key=material+spin;t=time.perf_counter()
  model,records=CachedFrozenSigmaResponse.from_export(paths[0]/('response'+spin),diagnostic_gap_cutoff=1e-10,cache_mb=64,precontract=False)
  h=[];o=[];gaps=[];clusters=[]
  for k,rec in enumerate(records):
   coeff=rec.coefficients[:,:,0]
   band=[]
   for stem,magic in [('native.ham',b'EDMFTHV1'),('native',b'EDMFTOV1')]:
    m,p=[aligned(rec,read_basis_matrix(q/f'{stem}.k{rec.kpoint:08d}.bin',magic)) for q in paths[1:]]
    band.append(herm(coeff.conj().T@((p-m)/.002)@coeff))
   dh,do=band;h.append(herm(dh-.5*(do*rec.energies[:,None]+do*rec.energies[None,:])));o.append(do)
   gap=model.eps[k,model.W][None,:]-model.eps[k,model.Q][:,None];gaps.append(gap)
   seen=set();within=set(model.W)
   for n in model.W:
    group=tuple(np.where(abs(rec.energies-rec.energies[n])<1e-8)[0])
    if any(i not in within for i in group) and group not in seen:
     seen.add(group);ix=np.array(group)
     clusters.append(dict(k_index=k+1,bands=(ix+1).tolist(),weight=float(model.kw[k]),raw_weight=float(np.linalg.norm(model.raw[k,ix])**2),normalized_weight=float(np.linalg.norm(model.apply_x(model.raw[k,ix][None],model.X))**2),h_block_norm=float(np.linalg.norm(h[-1][np.ix_(ix,ix)])),h_block_anisotropy=float(np.linalg.norm(h[-1][np.ix_(ix,ix)]-np.eye(len(ix))*np.trace(h[-1][np.ix_(ix,ix)])/len(ix)))))
  h=np.asarray(h);model.prepare();gradient=model._diagnostic_projector_gradient
  rows=[]
  for cutoff in [1e-12,1e-10,1e-8,1e-6,1e-4,1e-3]:
   model.diagnostic_gap_cutoff=cutoff;model._denominators=model.window_denominators()
   hb,_=model.embedding_pullback(gradient)
   val=real_inner(hb,h)
   rows.append(dict(cutoff_Ry=cutoff,omitted_pairs=sum(int(np.sum(abs(g)<=cutoff)) for g in gaps),projector_window_mean_GPa=val*conversion))
  model.diagnostic_gap_cutoff=1e-10;model._denominators=model.window_denominators()
  for k in set(x['k_index']-1 for x in clusters): model._denominators[k][:]=np.inf
  hb,_=model.embedding_pullback(gradient);allk=real_inner(hb,h)*conversion
  # Rotate a COMPLETE degenerate multiplet, then apply the SAME fixed index
  # selection and recompute global X. Compare physical U in the original gauge.
  base=model.lattice_scalar();probes=[]
  for seed in [20261109,20261110,20261111]:
   rng=np.random.default_rng(seed);q=np.zeros_like(model.raw);q[:,model.W]=model.raw[:,model.W]
   for cluster in clusters:
    k=cluster['k_index']-1;ix=np.array(cluster['bands'])-1;n=len(ix)
    rot=np.linalg.qr(rng.normal(size=(n,n)))[0]
    selected=np.array([int(i in model.W) for i in ix])
    projector=(rot*selected)@rot.T
    q[k,ix]=np.einsum('ij,jac->iac',projector,model.raw[k,ix])
   ss=model.overlap(q)+model.offset;x=np.zeros_like(model.X)
   for ic,n in enumerate(model.meta.dimensions):
    x[:n,:n,ic]=np.eye(n);act=np.flatnonzero(model.meta.active[:n,ic]);ii=np.ix_(act,act);x[:,:,ic][ii]=invsqrt(ss[:,:,ic][ii])
   u=model.meta.compact(model.apply_x(q,x))
   delta=model.lattice_scalar(u=u)-base
   probes.append(dict(seed=seed,projector_change_norm=float(np.linalg.norm(u-model.uc)),lattice_scalar_change_Ry=delta,apparent_mean_GPa_at_eta001=delta/.001*conversion))
  density=model.central_density();nr=float(np.trace(density,axis1=1,axis2=2).real.sum())
  results[key]=dict(material=material,spin=spin or 'up_or_nonmagnetic',volume_bohr3=vol,clusters=clusters,threshold_sweep=rows,omit_all_window_motion_at_affected_k_GPa=allk,extra_change_omitting_affected_k_GPa=allk-rows[1]['projector_window_mean_GPa'],rotation_probes=probes,raw_retained_electrons=nr,compressibility_per_Ry=model.compressibility,seconds=time.perf_counter()-t,scope='window+normalization chain only, contracted against finite-difference H/O from relaxed endpoint potentials; not total stress')
  (ROOT/'window_probe.json').write_text(json.dumps(results,indent=2)+'\n')
  np.savez(ROOT/(key+'_window_direction.npz'),h=h,o=np.asarray(o))
  print(key,json.dumps({k:v for k,v in results[key].items() if k!='clusters'}),flush=True)
