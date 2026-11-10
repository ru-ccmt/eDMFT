from pathlib import Path
import re,json,hashlib,struct,subprocess,shutil,os
import numpy as np
root=Path('/private/tmp/dft_after_remaining_dftks_20261107')
def vectors(folder,case):
 nat=int((folder/(case+'.struct')).read_text().splitlines()[1][27:30]);s=(folder/(case+'.output1')).read_text();nlos=set(map(int,re.findall(r'number of local orbitals, nlo \(hamilt\)\s+(\d+)',s)))
 if len(nlos)!=1: raise ValueError((folder,nlos))
 nlo=nlos.pop();out=[]
 with (folder/(case+'.vector')).open('rb') as f:
  def rec():
   raw=f.read(4)
   if not raw:return None
   n=struct.unpack('<i',raw)[0];data=f.read(n);assert len(data)==n and f.read(4)==raw;return data
  for i in range(2*nat):rec()
  while (h:=rec()) is not None:
   x,y,z,name,n,nb,w=struct.unpack('<3d10s2id',h);g=rec();assert len(g)==12*n
   pw=np.frombuffer(g,dtype='<i4').reshape(n,3)[:n-nlo];rows=sorted(map(tuple,pw.tolist()));digest=hashlib.sha256(np.array(rows,dtype='<i4').tobytes()).hexdigest()
   for j in range(nb):rec();rec()
   out.append(dict(k=[x,y,z],weight=w,basis=n,bands=nb,pw_hash=digest))
 return out
res=json.loads((root/'results.json').read_text())
env=os.environ|{'WIENROOT':'/Users/haule/wien2k_24','WIEN_DMFT_ROOT':'/Users/haule/dbin','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'}
for label,v in res.items():
 case=v['case'];d=root/label;ref=[Path(p) for p in v['reference_dirs']]
 data=[vectors(p,case) for p in ref]+[vectors(d,case)];keys=['k','weight','basis','bands','pw_hash'];v['basis_consistency']={key:all([a[key] for a in data[j]]==[a[key] for a in data[3]] for j in [0,2]) for key in keys};v['original_center_vector_nk']=len(data[1]);v['nk']=len(data[3]);v['nbands']=sorted(set(a['bands'] for a in data[3]));v['basis_range']=[min(a['basis'] for a in data[3]),max(a['basis'] for a in data[3])]
 assert all(v['basis_consistency'].values()),(label,v['basis_consistency'])
 for p in (d/'after').iterdir():
  if p.is_file():shutil.copy2(p,d/p.name)
 for name,cmd in [('merge',['/Users/haule/dbin/merge_scf.py']),('mixer',['/Users/haule/wien2k_24/x','mixer','-f',case])]:
  with (d/(name+'.log')).open('w') as f:subprocess.run(cmd,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 s=(d/(case+'.scfm')).read_text();pressure=float(re.findall(r'^:PRESSURE:\s*([-+0-9.]+)',s,re.M)[-1]);v['mixer_mean_GPa']=pressure;v['fresh_charge_distance']=float(re.findall(r'^:DIS.*\s([0-9.]+)\s*$',s,re.M)[-1]);assert abs(pressure-v['after']['mean_GPa'])<1e-5
 print(label, 'mean', pressure,'FD',v['numerical_GPa'],'difference',v['after_minus_numeric_GPa'],'DIS',v['fresh_charge_distance'],'basis',v['basis_consistency'],flush=True)
(root/'validated_results.json').write_text(json.dumps(res,indent=2)+'\n')
