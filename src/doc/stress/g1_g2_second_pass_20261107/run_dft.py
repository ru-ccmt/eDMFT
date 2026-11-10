from pathlib import Path
import subprocess,shutil,os,json,time
from concurrent.futures import ThreadPoolExecutor
root=Path(__file__).resolve().parent
wien=Path('/Users/haule/wien2k_24')
env=dict(os.environ,OMP_NUM_THREADS='2',VECLIB_MAXIMUM_THREADS='2',OPENBLAS_NUM_THREADS='1',WIENROOT=str(wien))
for k in list(env):
 if k.startswith(('WIEN_STRESS_','EDMFT_STRESS_')):del env[k]
def one(case,source,exe):
 folder=root/'runs'/('DFT_'+case);shutil.copytree(source,folder,ignore=shutil.ignore_patterns('audit'))
 (folder/'audit').mkdir()
 for mode,program,extra in [('current',wien/exe,{})]:
  t=time.monotonic();print('START',case,mode,flush=True)
  with (folder/(mode+'.log')).open('w') as f:
   p=subprocess.run([str(program),'lapw2.def'],cwd=folder,env=env|extra,stdout=f,stderr=subprocess.STDOUT)
  assert p.returncode==0,(case,mode,p.returncode)
  assert 'LAPW2 END' in (folder/(mode+'.log')).read_text(),(case,mode)
  out=folder/'audit'/mode;out.mkdir()
  for suff in ['scf0','scfc','scf2','clmval','output2']:
   f=folder/(case+'.'+suff)
   if f.exists():shutil.copy2(f,out/f.name)
  print('DONE',case,mode,round(time.monotonic()-t,2),flush=True)
 return str(folder)
items=[('FeSe',Path('/private/tmp/lapw2_stress_audit_20261107/runs/FeSe_DFT_RELA'),'lapw2'),('SrVO3',Path('/private/tmp/lapw2_stress_audit_20261107/runs/SrVO3_DFT_center'),'lapw2c')]
with ThreadPoolExecutor(max_workers=2) as p:out=list(p.map(lambda x:one(*x),items))
(root/'dft_runs.json').write_text(json.dumps(out,indent=2)+'\n')
