from pathlib import Path
import os,shutil,subprocess,time,json,hashlib,concurrent.futures
root=Path(__file__).parent
prior=Path('/private/tmp/projector_reaudit_20261107')
exe=Path('/Users/haule/projects/eDMFT_new/src/dmft2/dmft2')
env=os.environ|{'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1',
 'WIENROOT':'/Users/haule/wien2k_24','WIEN_DMFT_ROOT':'/Users/haule/dbin','PYTHONDONTWRITEBYTECODE':'1'}
for key in list(env):
 if key.startswith(('EDMFT_STRESS_','WIEN_STRESS_','EDMFT_PROJECTOR_')):del env[key]
env['EDMFT_STRESS_PROJECTOR_DIAG']='1'
def run(name):
 src=prior/name;dst=root/name;dst.mkdir(exist_ok=True)
 for p in src.iterdir():
  if p.is_file() and not p.is_symlink():shutil.copy2(p,dst/p.name)
 ranks=1 if name.startswith('Fe_atom') else 4
 print('START',name,flush=True);start=time.monotonic()
 with (dst/'new.log').open('w') as log:
  p=subprocess.run(['/opt/homebrew/bin/mpirun','-n',str(ranks),str(exe),'dmft2.def'],cwd=dst,env=env,stdout=log,stderr=subprocess.STDOUT)
 print('END',name,p.returncode,round(time.monotonic()-start,2),flush=True)
 if p.returncode:raise RuntimeError(name)
 print(name,(dst/'projector_strain_fd.dat').read_text(),flush=True)
 return name,{'seconds':time.monotonic()-start,'exit':p.returncode}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 result=dict(pool.map(run,['SrVO3','FeSe','Fe_atom1','Fe_atom2']))
result['executable_sha256']=hashlib.sha256(exe.read_bytes()).hexdigest()
(root/'runs.json').write_text(json.dumps(result,indent=2)+'\n')
