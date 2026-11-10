from pathlib import Path
import subprocess,shutil,os,json,re,hashlib
root=Path('/private/tmp/stress_remaining_audit_20261107');results={}
for label in ['SrVO3_DFT','SrVO3_DMFT','FeSe_DFT','FeSe_DMFT']:
 d=root/label;case='SrVO3' if label.startswith('SrVO3') else 'FeSe';baseline=d/'baseline';baseline.mkdir(exist_ok=True)
 for suf in ['scf0','output0','vsp','vns','scfc','outputc','clmcor']:
  p=d/(case+'.'+suf)
  if p.exists() and not (baseline/p.name).exists():shutil.copy2(p,baseline/p.name)
 env=os.environ.copy();env.update(WIEN_STRESS_DEBUG='1',OMP_NUM_THREADS='1')
 with (d/'stage_lapw0.log').open('w') as f:subprocess.run([str(root/'SRC_lapw0/lapw0'),'lapw0.def'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 def tensor(p):
  s=p.read_text();return [[float(x) for x in re.findall(r'^:STR_DENS'+f'{i:03d}'+r':\s*(.*)',s,re.M)[-1].split()[:3]] for i in [1,2,3]]
 old=tensor(baseline/(case+'.scf0'));new=tensor(d/(case+'.scf0'))
 delta=[(new[i][i]-old[i][i])*14710.5164057 for i in range(3)]
 results[label]={'delta_diagonal_GPa':delta,'delta_mean_GPa':sum(delta)/3,'potential_bytes_identical':{s:(baseline/(case+'.'+s)).read_bytes()==(d/(case+'.'+s)).read_bytes() for s in ['vsp','vns']}}
 print(label,results[label])
 (d/'stage').mkdir(exist_ok=True)
 for suf in ['scf0','output0','vsp','vns']:shutil.copy2(d/(case+'.'+suf),d/'stage'/(case+'.'+suf))
(root/'density_stage_results.json').write_text(json.dumps(results,indent=2)+'\n')
