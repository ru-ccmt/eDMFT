from pathlib import Path
import os,sys,re,shutil,subprocess,time,json,concurrent.futures
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1'
import numpy as np
from scipy.interpolate import CubicSpline
r=Path(__file__).parent/'SrVO3_affine';r.mkdir(exist_ok=True)
c=Path('/private/tmp/projector_stationarity_probe_20261107/SrVO3')
prior=Path('/private/tmp/projector_actual_stress_20261107')
struct=(c/'SrVO3.struct').read_text().splitlines();cell=np.array([float(struct[3][i:i+10]) for i in range(0,30,10)]);rad=[]
for line in struct:
 m=re.search(r'NPT=\s*(\d+)\s+R0=\s*([\d.]+)\s+RMT=\s*([\d.]+)',line)
 if m:rad.append(np.linspace(np.log(float(m[2])),np.log(float(m[3])),int(m[1])))
pat=re.compile(r'[-+]?\d+\.\d+E[-+]\d+')
def affine(text,lam):
 a=text.splitlines();iat=-1;i=0
 while i<len(a):
  if 'ATOMNUMBER' in a[i]:iat=int(a[i].split('=')[1])-1
  if 'CLM(R) FOR L' not in a[i]:i+=1;continue
  i+=1
  while not a[i].strip():i+=1
  start=i;v=[]
  while len(v)<len(rad[iat]):v.extend(map(float,pat.findall(a[i])));i+=1
  assert len(v)==len(rad[iat]);v=np.array(v)
  # Stored radial density is r^2 rho_lm(r), giving lambda^-1 here.
  out=CubicSpline(rad[iat],v)(rad[iat]-np.log(lam))/lam
  a[start:i]=['   '+''.join(f'{x:19.12E}' for x in out[j:j+4]) for j in range(0,len(out),4)]
 j=next(i for i,x in enumerate(a) if 'NUMBER OF PW' in x);n=int(a[j].split()[0])
 for i in range(j+1,j+1+n):
  line=a[i];z=np.array([float(line[18:37]),float(line[37:56])])/lam**3
  a[i]=line[:18]+''.join(f'{x:19.12E}' for x in z)
 return '\n'.join(a)+'\n'
env=os.environ|{'WIENROOT':'/Users/haule/wien2k_24','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'}
def run(spec):
 name,h=spec;d=r/name;d.mkdir(exist_ok=True)
 for suf in ['struct','clmsum','clmcor','in0','in1c','inm','inc','klist','kgen','rotlm']:
  p=c/('SrVO3.'+suf)
  if p.exists():shutil.copy2(p,d/p.name)
 ss=struct.copy();new=cell*(1+h);ss[3]=''.join(f'{x:10.6f}' for x in new)+ss[3][30:]
 actual=np.array([float(ss[3][i:i+10]) for i in range(0,30,10)]);lam=float(np.mean(actual/cell));(d/'SrVO3.struct').write_text('\n'.join(ss)+'\n')
 p=d/'SrVO3.clmsum';p.write_text(affine(p.read_text(),lam))
 p=d/'SrVO3.in0';a=p.read_text().splitlines();a[0]='TOT '+a[0][4:];a[1]=a[1].replace('NR2V','R2V ');p.write_text('\n'.join(a)+'\n')
 t=time.monotonic()
 with (d/'lapw0.log').open('w') as f:p=subprocess.run(['/Users/haule/wien2k_24/x','lapw0','-f','SrVO3'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 assert p.returncode==0;time0=time.monotonic()-t
 # A fixed-geometry LAPW1 probe measures the potential matrix element only.
 b=r/(name+'_basis');b.mkdir(exist_ok=True)
 for fname in ['SrVO3.struct','SrVO3.in1c','SrVO3.klist','SrVO3.kgen','SrVO3.rotlm','lapw1.def']:shutil.copy2(c/fname,b/fname)
 for suf in ['vsp','vns']:shutil.copy2(d/('SrVO3.'+suf),b/('SrVO3.'+suf))
 t=time.monotonic()
 with (b/'lapw1.log').open('w') as f:p=subprocess.run(['/Users/haule/wien2k_24/lapw1c','lapw1.def'],cwd=b,env=env|{'WIEN_LAPW1_BASIS_VSP':str(c/'SrVO3.vsp')},stdout=f,stderr=subprocess.STDOUT)
 assert p.returncode==0
 out={'nominal_step':h,'actual_step':lam-1,'lapw0_seconds':time0,'lapw1_probe_seconds':time.monotonic()-t}
 (d/'metadata.json').write_text(json.dumps(out,indent=2));print(name,out,flush=True);return name,out
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as p:out=dict(p.map(run,[('center',0.),('p001',.001),('p0005',.0005),('p00025',.00025)]))
 (r/'timings.json').write_text(json.dumps(out,indent=2))
