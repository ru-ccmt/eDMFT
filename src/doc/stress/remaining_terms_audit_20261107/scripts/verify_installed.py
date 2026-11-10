from pathlib import Path
import os,subprocess,re,numpy as np,json
r=Path(__file__).resolve().parent;out={}
for label in ['SrVO3_DFT','SrVO3_DMFT','FeSe_DFT','FeSe_DMFT']:
 d=r/label;c=label.split('_')[0];out[label]={}
 for exe,suffix,tag in [('lapw0','scf0','STR_DENS'),('lcore','scfc','STR_CORE')]:
  with (d/('installed_'+exe+'.log')).open('w') as f:subprocess.run(['/Users/haule/wien2k_24/'+exe,exe+'.def'],cwd=d,env=os.environ|{'OMP_NUM_THREADS':'1','WIEN_STRESS_DEBUG':'1'},stdout=f,stderr=subprocess.STDOUT,check=True)
  def get(p):
   s=p.read_text();return np.array([list(map(float,re.findall('^:'+tag+f'{i:03d}:'+r'\s*(.*)',s,re.M)[-1].split())) for i in [1,2,3]])
  err=float(np.max(abs(get(d/(c+'.'+suffix))-get(d/'stage'/(c+'.'+suffix)))))
  assert err<1e-12,(label,exe,err);out[label][exe+'_max_difference_RyBohr3']=err
 print(label,out[label],flush=True)
(r/'installed_verification.json').write_text(json.dumps(out,indent=2)+'\n')
