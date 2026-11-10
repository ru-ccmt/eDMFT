"""Analyse the third native step without repeating analytic preparation."""
from pathlib import Path
import json,sys
import numpy as np
root=Path(__file__).resolve().parent
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from lapw_strain_potential import unequal_centered_derivative
from wien_response_density import Density
from wstruct import Struct
for case in ['SrVO3','FeSe']:
 work=root/'three_k'/case
 ref=Path('/private/tmp/edmft_screened_native_20261108/three_k')/case
 source=Path('/private/tmp/edmft_full_k_response_20261108')/case/'charge/0'
 template=Density(source/(case+'.clmsum'),source/(case+'.struct'))
 folders=[root/'native_extra'/case/f'source_{label}_6.25e-05' for label in ['m','p']]
 offsets=[json.loads((p/'done.json').read_text())['actual_scale']-1 for p in folders]
 derivative={}
 for suffix in ['.clmval','.clmcor','normalized.npy']:
  def read(folder):
   if suffix.endswith('.npy'):return np.load(folder/suffix)
   return template.zeros().add(Density(folder/(case+suffix),folder/(case+'.struct'))).vector()
  derivative[suffix]=unequal_centered_derivative(read(folders[0]),read(ref/'central'),read(folders[1]),*offsets)
 np.savez_compressed(work/'native_6.25e-05.npz',**derivative)
 report=json.loads((work/'REPORT.json').read_text())
 analytic=np.load(work/'analytic.npz');source_derivative=analytic['source']
 old=np.load(work/'native_0.000125.npz')['normalized.npy']
 old_offsets=report['steps']['0.000125']['actual_offsets']
 ratio=((old_offsets[1]-old_offsets[0])/(offsets[1]-offsets[0]))**2
 extrap=(ratio*derivative['normalized.npy']-old)/(ratio-1)
 nrad=sum(map(len,template.rad.values()))
 out=dict(actual_offsets=offsets,source_relative=float(np.linalg.norm(source_derivative-derivative['normalized.npy'])/np.linalg.norm(derivative['normalized.npy'])),
          source_PW_relative=float(np.linalg.norm((source_derivative-derivative['normalized.npy'])[nrad:])/np.linalg.norm(derivative['normalized.npy'][nrad:])),
          richardson_relative=float(np.linalg.norm(source_derivative-extrap)/np.linalg.norm(extrap)),
          richardson_PW_relative=float(np.linalg.norm((source_derivative-extrap)[nrad:])/np.linalg.norm(extrap[nrad:])))
 s0=Struct();s0.ReadStruct(str(source/(case+'.struct')));cells=[]
 for p in folders:
  s=Struct();s.ReadStruct(str(p/(case+'.struct')));cells.append(np.array([s.a,s.b,s.c]))
 out['strain']=np.diag((cells[1]-cells[0])/(np.array([s0.a,s0.b,s0.c])*(offsets[1]-offsets[0]))).tolist()
 assert np.max(abs(np.array(out['strain'])-report['strain']))<1e-9
 report['third_step']=out
 (work/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
 print(case,json.dumps(out),flush=True)
