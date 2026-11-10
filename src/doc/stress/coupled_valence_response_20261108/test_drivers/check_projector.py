from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'stage'),'/Users/haule/projects/eDMFT_new/src/python','/private/tmp/edmft_augmentation_response_20261108']
from material_augmentation import Material
sys.path.insert(0,str(ROOT/'stage'))
from lapw_radial_native import NativeRadialContext
from lapw_augmentation_response import PreparedAugmentationChannel
from frozen_sigma_response import FrozenSigmaResponse
from indmffile import Indmfl
case=sys.argv[1];central=Path('/private/tmp/edmft_spatial_response_20261108')/case/'central'
mat=Material(case,central);s=mat.s
model,records=FrozenSigmaResponse.from_export(central/'response')
inl=Indmfl(str(central/case));inl.read()
projectors={}
with (central/'projectorw.dat').open() as f:
 for _ in range(int(f.readline().split()[1])):
  nr,nr0,a,l=map(int,f.readline().split()[1:]);rows=np.array([list(map(float,f.readline().split())) for _ in range(nr)])
  assert nr==nr0
  projectors[a-1,l]=rows[:,1:3]
raw=np.zeros_like(model.raw);operations={};contexts={};rows=[]
site_atom=np.repeat(np.arange(s.nat),s.mult);first=np.r_[0,np.cumsum(s.mult)[:-1]]
for ic,orbs in inl.cix.items():
 assert len(orbs)==1
 site,l,qsplit=orbs[0];site-=1;a=int(site_atom[site]);n=s.jrj[a]
 if (a,l) not in contexts:
  e=mat.header[a][l];lapw=bool(e<150);e=e if lapw else e-200;elo=mat.lo[a][l,:mat.ilo[a][l]]
  ctx=NativeRadialContext(mat.v[a][0,0],r0=s.r0[a],dx=np.log(s.rmt[a]/s.r0[a])/(n-1),rmt=s.rmt[a],energy=e,
    lo_energies=elo,ell=l,charge=s.Znuc[a],relativistic=bool('RELA' in s.mode),lapw=lapw,library=ROOT/'liblapw_response.so')
  contexts[a,l]=ctx
 ctx=contexts[a,l]
 rix=np.einsum('rci,rc,rc->i',ctx.primitive[:,:2],ctx.projector_metric[:,:2],projectors[a,l])
 locrot,axes,shift,rmt2=inl.atoms[site+1]
 assert shift is None
 axes=np.array(axes) if len(axes) else s.rotloc[a]
 assert axes.shape==(3,3)
 for rotation,tau in zip(s.timat,s.tau):
  diff=rotation@s.pos[first[a]]+tau-s.pos[site];diff-=np.rint(diff)
  if max(abs(diff))<2e-7:break
 effective=rotation@s.pos[first[a]]+tau
 for k,(ev,c,labels) in enumerate(mat.cv):
  with (central/f'native.k{k+1:08d}.bin').open('rb') as f:f.seek(44);kf=np.fromfile(f,'f8',3)
  frac=labels.T+kf;q=2*np.pi*frac/np.array([s.a,s.b,s.c]);qn=np.linalg.norm(q,axis=1);npw=len(c)-mat.nlo_total
  phase=np.exp(2j*np.pi*(frac@effective));qloc=(axes@rotation.T@q.T).T
  groups=[npw+mat.lo_indices[a,l,j] for j in range(ctx.nlo)]
  ch=PreparedAugmentationChannel(l,s.rmt[a],s.a*s.b*s.c,qloc,qn,phase,c,npw,groups)
  z=ch.value(ctx.primitive[-1],ctx.lo_coefficients,lapw=ctx.lapw)
  p=np.einsum('min,i,ma->na',z.conj(),rix,inl.cftrans[ic].T)
  raw[k,:,:(2*l+1),ic-1]=p;operations[ic-1,k]=(ch,ctx,projectors[a,l],inl.cftrans[ic].T)
  ref=model.raw[k,:,:(2*l+1),ic-1]
  row=dict(cix=int(ic),k=k,relative=float(np.linalg.norm(p-ref)/np.linalg.norm(ref)),maximum=float(np.max(abs(p-ref))))
  rows.append(row);print(row,flush=True)
(ROOT/(case+'_projector.json')).write_text(json.dumps(rows,indent=2)+'\n')
