from pathlib import Path
import os,subprocess,json,re,numpy as np,shutil
r=Path(__file__).resolve().parent;out={}
def tensor(p,tag):
 s=p.read_text();return np.array([list(map(float,re.findall('^:'+tag+f'{i:03d}:'+r'\s*(.*)',s,re.M)[-1].split())) for i in [1,2,3]])
for label in ['SrVO3_DFT','SrVO3_DMFT','FeSe_DFT','FeSe_DMFT']:
 d=r/label;case=label.split('_')[0];st=r/'SRC_lapw0/lapw0';co=r/'SRC_lcore/lcore';out[label]={}
 for mode,exe,env in [('lcore',co,{}),('legacy',st,{'WIEN_STRESS_GGA_MODE':'LEGACY'}),('omp2',st,{'OMP_NUM_THREADS':'2'})]:
  with (d/(mode+'.log')).open('w') as f:subprocess.run([str(exe),'lcore.def' if mode=='lcore' else 'lapw0.def'],cwd=d,env=os.environ|{'OMP_NUM_THREADS':'1','WIEN_STRESS_DEBUG':'1'}|env,stdout=f,stderr=subprocess.STDOUT,check=True)
  if mode=='lcore':
   out[label]['core_max_change']=float(np.max(abs(tensor(d/(case+'.scfc'),'STR_CORE')-tensor(d/'baseline'/(case+'.scfc'),'STR_CORE'))))
   out[label]['core_density_identical']=(d/(case+'.clmcor')).read_bytes()==(d/'baseline'/(case+'.clmcor')).read_bytes()
   shutil.copy2(d/(case+'.scfc'),d/'stage'/(case+'.scfc'))
  else:
   ref=d/('stage' if mode=='omp2' else 'baseline')/(case+'.scf0');delta=tensor(d/(case+'.scf0'),'STR_DENS')-tensor(ref,'STR_DENS');out[label][mode+'_max_delta_GPa']=float(np.max(abs(delta))*14710.5164057)
   out[label][mode+'_mean_delta_GPa']=float(np.trace(delta)*14710.5164057/3)
   out[label][mode+'_potential_bytes_identical']={ext:(d/(case+'.'+ext)).read_bytes()==(d/'baseline'/(case+'.'+ext)).read_bytes() for ext in ['vsp','vns']}
 print(label,out[label],flush=True)
(r/'runtime_mode_results.json').write_text(json.dumps(out,indent=2)+'\n')
