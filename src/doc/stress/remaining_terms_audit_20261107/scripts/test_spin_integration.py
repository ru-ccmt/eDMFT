from pathlib import Path
import subprocess,shutil,re,json,os,numpy as np
r=Path(__file__).resolve().parent;src=r/'FeSe_DFT';d=r/'spin_integration';d.mkdir(exist_ok=True)
for p in src.iterdir():
 if p.is_file() and (p.name.startswith('FeSe.') or p.name=='lapw0.def'):shutil.copy2(p,d/p.name)
rho=(src/'FeSe.clmsum').read_text().splitlines();npt=[int(t[15:20]) for t in (src/'FeSe.struct').read_text().splitlines() if 'NPT=' in t];pw=next(i for i,s in enumerate(rho) if 'NUMBER OF PW' in s);npw=int(rho[pw].split()[0])
def scale(f):
 rr=rho.copy();atom=0
 for i,t in enumerate(rho[:pw]):
  if 'ATOMNUMBER' in t:atom=int(t.split('=')[-1])-1
  if 'CLM(R)' in t:
   k=i+2;left=npt[atom]
   while left:
    count=min(4,left);v=[float(rr[k][3+19*j:3+19*(j+1)]) for j in range(count)];rr[k]='   '+''.join(f'{x*f:19.12E}' for x in v);k+=1;left-=count
 for i in range(pw+1,pw+npw+1):
  t=rho[i];rr[i]=t[:18]+''.join(f'{float(t[18+19*j:18+19*(j+1)])*f:19.12E}' for j in range(2))
 return '\n'.join(rr)+'\n'
def tensor(s):return np.array([list(map(float,re.findall('^:STR_DENS'+f'{i:03d}:'+r'\s*(.*)',s,re.M)[-1].split())) for i in [1,2,3]])
out={};ref=tensor((src/'stage/FeSe.scf0').read_text());tensors={}
for name,factor in [('unpolarized',.5),('polarized',.65),('swapped',.35)]:
 (d/'FeSe.clmup').write_text(scale(factor));(d/'FeSe.clmdn').write_text(scale(1-factor))
 with (d/(name+'.log')).open('w') as f:subprocess.run(['/Users/haule/wien2k_24/lapw0','lapw0.def'],cwd=d,env=os.environ|{'OMP_NUM_THREADS':'1'},stdout=f,stderr=subprocess.STDOUT,check=True)
 s=(d/'FeSe.scf0').read_text();(d/(name+'.scf0')).write_text(s);tensors[name]=tensor(s)
 out[name]={'stress':tensors[name].tolist()}
err=float(np.max(abs(tensors['unpolarized']-ref))*14710.5164057);spinerr=float(np.max(abs(tensors['polarized']-tensors['swapped']))*14710.5164057)
out['unpolarized_representation_error_GPa']=err;out['spin_exchange_error_GPa']=spinerr
assert err<1e-6 and spinerr<1e-6,out
(r/'spin_integration_results.json').write_text(json.dumps(out,indent=2)+'\n');print(err,spinerr)
