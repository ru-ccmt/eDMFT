"""Analytic central-gauge H/O/raw geometry vs independent native endpoints."""
from pathlib import Path
import json,sys,time
import numpy as np
root=Path(__file__).resolve().parent
if root.name=='reproduce':root=root.parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_geometry_response import native_band_geometry_source
from lapw_valence_state import load_native_valence_state
from lapw_native_snapshot import read_basis_matrix
from lapw_strain_potential import unequal_centered_derivative
from frozen_sigma_response import FrozenSigmaResponse
from lapw_projector_functional import projector_functional_gradient

def run(case):
 source=Path('/private/tmp/edmft_full_k_response_20261108')/case/'charge/0'
 work=root/case
 with load_native_valence_state(source,case,complex_vectors=(case=='SrVO3'),energy_policy='fixed',matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so') as v:
  cell=np.array([v.structure.a,v.structure.b,v.structure.c]);reports={}
  g,coords=projector_functional_gradient(v.complement)
  for step in [.00025,.000125]:
   folders=[work/f'{label}_{step}' for label in ['m','p']]
   endpoints=[json.loads((p/'endpoint.json').read_text()) for p in folders]
   hm,hp=[r['eta'] for r in endpoints]
   e=np.diag((np.array(endpoints[1]['cell'])-endpoints[0]['cell'])/(cell*(hp-hm)))
   start=time.monotonic();analytic,report=native_band_geometry_source(v,source,case,e)
   report['geometry_seconds']=time.monotonic()-start
   ah,ao,ap,ahq,aoq=v.bands.output.unpack(analytic)
   h=[];o=[];hq=[];oq=[]
   for k,base in enumerate(v.bases):
    c=np.column_stack((base.retained,base.complement));nr=base.retained.shape[1]
    outputs=[]
    for prefix,magic in [('native.ham',b'EDMFTHV1'),('native',b'EDMFTOV1')]:
     center=read_basis_matrix(source/f'{prefix}.k{k+1:08d}.bin',magic)
     mats=[read_basis_matrix(p/f'{prefix}.k{k+1:08d}.bin',magic) for p in folders]
     for m in mats:np.testing.assert_array_equal(m.labels,center.labels)
     # Native covariant H/O are independent of the endpoint eigenvectors.
     derivative=unequal_centered_derivative(mats[0].matrix.real,center.matrix.real,mats[1].matrix.real,hm,hp)+1j*unequal_centered_derivative(mats[0].matrix.imag,center.matrix.imag,mats[1].matrix.imag,hm,hp)
     outputs.append(c.conj().T@(derivative@c[:,:nr]))
    h.append(outputs[0][:nr]);hq.append(outputs[0][nr:]);o.append(outputs[1][:nr]);oq.append(outputs[1][nr:])
   models=[FrozenSigmaResponse.from_export(p/'response')[0] for p in folders]
   raw=(unequal_centered_derivative(models[0].raw.real,v.model.raw.real,models[1].raw.real,hm,hp)+
        1j*unequal_centered_derivative(models[0].raw.imag,v.model.raw.imag,models[1].raw.imag,hm,hp))
   numeric=v.bands.output.pack(np.array(h),np.array(o),raw,hq,oq)
   report['relative']={}
   def norms(label,left,right):
    av=np.concatenate([x.ravel() for x in left]);bv=np.concatenate([x.ravel() for x in right])
    report['relative'][label]=dict(error=float(np.linalg.norm(av-bv)/np.linalg.norm(bv)),norm=float(np.linalg.norm(bv)),max_error=float(np.max(abs(av-bv))))
   norms('H_RR',ah,h);norms('H_QR',ahq,hq);norms('O_RR',ao,o);norms('O_QR',aoq,oq);norms('raw',ap,raw)
   report['projector_geometry_only_Ry']=dict(analytic=float(g@analytic),native=float(g@numeric),difference=float(g@(analytic-numeric)))
   report['step']=step;report['actual_steps']=[hm,hp]
   np.savez(work/f'comparison_{step}.npz',analytic=analytic,native=numeric,gradient=g)
   reports[str(step)]=report
   print(case,step,json.dumps(report),flush=True)
  a,b=[reports[h] for h in ['0.00025','0.000125']]
  ha=(a['actual_steps'][1]-a['actual_steps'][0])/2;hb=(b['actual_steps'][1]-b['actual_steps'][0])/2
  d1=np.load(work/'comparison_0.00025.npz');d2=np.load(work/'comparison_0.000125.npz');ratio=(ha/hb)**2
  extrap=(ratio*d2['native']-d1['native'])/(ratio-1)
  result=dict(case=case,kpoints=len(v.records),preparation_seconds=v.preparation_seconds,steps=reports,source_fingerprints=v.fingerprints,
    richardson_joint_relative=float(np.linalg.norm(extrap-d2['analytic'])/np.linalg.norm(extrap)),
    richardson_projector_scalar_error_Ry=float(d2['gradient']@(d2['analytic']-extrap)))
  (work/'REPORT.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':run(sys.argv[1])
