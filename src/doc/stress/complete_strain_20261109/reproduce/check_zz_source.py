"""Independent native axial strain source, including searched radial energies."""
from pathlib import Path
import sys,os,json,shutil,subprocess
import numpy as np
root=Path(__file__).resolve().parent;sys.path.insert(0,'/Users/haule/dbin')
from lapw_potential_core_state import load_potential_core_state
from lapw_charge_response import NativeChargeNormalization,native_mixer_inversion
from wien_response_density import Density
from wien_frozen_sigma_response import vector_signature
from wstruct import Struct
from scipy.io import FortranFile,FortranEOFError
import struct
from lapw_valence_state import _snapshot_inputs
def physical_signature(path,nat,nlo):
 result=[]
 with FortranFile(path) as f:
  for _ in range(2*nat):f.read_record('u1')
  while True:
   try:head=f.read_record('u1')
   except FortranEOFError:break
   kx,ky,kz,name,nb,n,wt=struct.unpack('<ddd10siid',head.tobytes());g=f.read_record('i4').reshape(3,nb,order='F')
   result.append((kx,ky,kz,wt,n,nb,nlo,tuple(sorted(map(tuple,g[:,:nb-nlo].T)))))
   for _ in range(n):f.read_record('u1');f.read_record('u1')
 return result
case=sys.argv[1];pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz');s=pc.structure;target=float(np.dot(s.mult,s.Znuc));src=root/case/'native_live/center';out=root/case/'zz_source';out.mkdir(exist_ok=True)
env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))};env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
old=Path('/private/tmp/edmft_xc_core_response_20261108/native');jobs=[('lapw0',old/'SRC_lapw0/lapw0'),('lapw1',Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if case=='SrVO3' else 'lapw1')),('dmft2',Path('/Users/haule/dbin/dmft2_response')),('lcore',old/'SRC_lcore/lcore')]
analytic=np.load(root/case/'tensor/source1.npz')['source'];rows=[]
for h in [.00025,.000125]:
 values=[];steps=[]
 for sign in [-1,1]:
  d=out/(str(h)+('_m' if sign<0 else '_p'));d.mkdir(exist_ok=True)
  names=[case+'.'+z for z in ['struct','in0','in1','in1c','in2','in2c','inm','inc','indmfl','indmfi','indmf2','klist','kgen','rotlm']]+['EF.dat','Eorb_imp.dat','projectorw.dat','lapw0.def','lapw1.def','dmft2.def','lcore.def']+[p.name for p in src.glob('sig.inp*')]
  for name in names:
   if (src/name).exists():shutil.copy2(src/name,d/name)
  lines=(d/(case+'.struct')).read_text().splitlines();c=float(f'{s.c*(1+sign*h):.6f}');eta=c/s.c-1;steps.append(eta);lines[3]=lines[3][:20]+f'{c:10.6f}'+lines[3][30:];(d/(case+'.struct')).write_text('\n'.join(lines)+'\n');scale=(1+eta)**(1/3)
  rho=pc.input_density.transport(scale);template=Density(src/(case+'.clmsum'),src/(case+'.struct'));template.rad={k:rho.rad[k].copy() for k in template.rad};template.pw={k:rho.pw[k] for k in template.pw};template.write(d/(case+'.clmsum'))
  shutil.copy2('/private/tmp/edmft_screened_native_20261108/full_bz/'+case+'/central/lapw1.fixgvec',d/'lapw1.fixgvec');(d/'lapw1.fixgvec.read').touch();(d/'core_native.bin').unlink(missing_ok=True)
  for name,exe in jobs:
   with (d/(name+'.log')).open('w') as f:subprocess.run([str(exe),name+'.def'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
   if name in ['lapw1','dmft2'] and name.upper()+' END' not in (d/(name+'.log')).read_text():raise ValueError(str(d))
  ss=Struct();ss.ReadStruct(str(d/(case+'.struct')));raw=Density(d/(case+'.clmval'),d/(case+'.struct'));raw.add(Density(d/(case+'.clmcor'),d/(case+'.struct')))
  norm=NativeChargeNormalization(raw,ss,target,inversion=native_mixer_inversion(s)).central
  pulled=template.zeros().add(raw.from_vector(norm)).transport(scale,True)
  norm=NativeChargeNormalization(pulled,s,target,inversion=native_mixer_inversion(s)).central;values.append(norm)
  window=tuple(map(int,(d/(case+'.indmfl')).read_text().splitlines()[0].split()[:2]));same=vector_signature(d/(case+'.vector'),s.nat,window)==vector_signature(root/case/'current'/(case+'.vector'),s.nat,window)
  if not same:
   # LO columns carry auxiliary plane-wave labels chosen by LAPW1. Their
   # choice can change under axial strain although the full (2l+1)-dimensional
   # local-orbital subspace is identical. Compare the actual PW basis and LO
   # count, retaining the original window-gap check above.
   nlo=_snapshot_inputs(case,root/case/'current',case=='SrVO3').nlo_total
   if physical_signature(d/(case+'.vector'),s.nat,nlo)!=physical_signature(root/case/'current'/(case+'.vector'),s.nat,nlo):raise ValueError('physical PW/LO basis changed')
  np.save(d/'normalized.npy',norm)
 fd=(values[1]-values[0])/(steps[1]-steps[0]);row=dict(step=h,actual_interval=steps[1]-steps[0],relative_error=float(np.linalg.norm(fd-analytic)/np.linalg.norm(analytic)),absolute_error=float(np.linalg.norm(fd-analytic)));rows.append(row);print(case,row,flush=True)
 (out/'REPORT.json').write_text(json.dumps(dict(case=case,scope='explicit zz unmixed charge source; not independently relaxed axial stress',rows=rows),indent=2)+'\n')
