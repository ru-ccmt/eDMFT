from pathlib import Path
import os,shutil,subprocess,json,re,numpy as np
r=Path(__file__).parent
src=r/'Fe_atom1'
env=os.environ|{'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','WIENROOT':'/Users/haule/wien2k_24','WIEN_DMFT_ROOT':'/Users/haule/dbin'}
for key in list(env):
 if key.startswith(('EDMFT_STRESS_','WIEN_STRESS_','EDMFT_PROJECTOR_')):del env[key]
out={}
for label,exe in [('before',r/'before/dmft2'),('after',Path('/Users/haule/projects/eDMFT_new/src/dmft2/dmft2'))]:
 d=r/('DFT_control_'+label);d.mkdir(exist_ok=True)
 for p in src.iterdir():
  if p.is_file() and not p.is_symlink():shutil.copy2(p,d/p.name)
 p=d/'FeSe.indmf2';ls=p.read_text().splitlines();v=ls[1].split();v[2]='-2';ls[1]=' '.join(v);p.write_text('\n'.join(ls)+'\n')
 with (d/'control.log').open('w') as log:
  proc=subprocess.run(['/opt/homebrew/bin/mpirun','-n','1',str(exe),'dmft2.def'],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
 assert proc.returncode==0,(label,proc.returncode)
 s=(d/'FeSe.scf2').read_text();rows=[]
 for i in range(1,4):rows.append(list(map(float,re.findall(f'^:STR_VAL{i:03d}:\\s+([^\\n]+)',s,re.M)[-1].split()[:3])))
 out[label]=rows
out['max_difference_RyBohr3']=float(np.max(abs(np.array(out['before'])-np.array(out['after']))))
assert out['max_difference_RyBohr3']<1e-11,out
(r/'dft_control.json').write_text(json.dumps(out,indent=2)+'\n')
print('DFT-only one-k FeSe control',out['max_difference_RyBohr3'])
