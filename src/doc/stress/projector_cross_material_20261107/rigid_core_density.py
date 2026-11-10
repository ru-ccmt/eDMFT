from pathlib import Path
import re,concurrent.futures,json
import numpy as np
from scipy.interpolate import CubicSpline
import affine_density as m
m.r=Path(__file__).parent/'SrVO3_rigid_core';m.r.mkdir(exist_ok=True)
base_affine=m.affine
lines=(m.c/'SrVO3.clmcor').read_text().splitlines();cores=[];i=0
while i<len(lines):
 if 'CLM(R) FOR L' not in lines[i]:i+=1;continue
 i+=1
 while not lines[i].strip():i+=1
 atom=len(cores);v=[]
 while len(v)<len(m.rad[atom]):v.extend(map(float,m.pat.findall(lines[i])));i+=1
 cores.append(np.array(v))

def rigid(text,lam):
 a=base_affine(text,lam).splitlines();i=0;iat=-1
 while i<len(a):
  if 'ATOMNUMBER' in a[i]:iat=int(a[i].split('=')[1])-1
  if not re.search(r'CLM\(R\) FOR L\s+0\s+M=\s*0',a[i]):i+=1;continue
  i+=1
  while not a[i].strip():i+=1
  start=i;v=[]
  while len(v)<len(m.rad[iat]):v.extend(map(float,m.pat.findall(a[i])));i+=1
  core=cores[iat];correction=core-CubicSpline(m.rad[iat],core)(m.rad[iat]-np.log(lam))/lam;v=np.array(v)+correction
  a[start:i]=['   '+''.join(f'{x:19.12E}' for x in v[j:j+4]) for j in range(0,len(v),4)]
 return '\n'.join(a)+'\n'
m.affine=rigid
# Central density/potential is identical to the affine experiment; reuse it.
import shutil
for name in ['center','center_basis']:shutil.copytree(Path(__file__).parent/'SrVO3_affine'/name,m.r/name,dirs_exist_ok=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:out=dict(pool.map(m.run,[('p0005',.0005),('p00025',.00025)]))
out['center']=json.loads((m.r/'center'/'metadata.json').read_text());(m.r/'timings.json').write_text(json.dumps(out,indent=2))
print('core radial endpoints', [float(c[-1]) for c in cores],flush=True)
