from pathlib import Path
import os,shutil,subprocess,time,json,hashlib
root=Path(__file__).resolve().parent
items=[('MnO','/private/tmp/edmft_mno_comparison_20261107/rerun_controls/fixed_dense','dmft2.def'),('MnTe_up','/private/tmp/edmft_mnte_comparison_20261107/up','dmft2.def'),('MnTe_dn','/private/tmp/edmft_mnte_comparison_20261107/dn','dmft2dn.def')]
env=dict(os.environ,OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1',WIENROOT='/Users/haule/wien2k_24',WIEN_DMFT_ROOT='/Users/haule/dbin')
for k in list(env):
 if k.startswith(('WIEN_STRESS_','EDMFT_STRESS_')):del env[k]
records={}
for name,src,definition in items:
 folder=root/('DMFT_'+name);shutil.copytree(src,folder)
 # Definitions must remain local; do not permit an output path into original runs.
 for line in (folder/definition).read_text().splitlines():
  parts=line.split("'")
  if len(parts)>1:assert not Path(parts[1]).is_absolute(),line
 start=time.monotonic();print('START',name,flush=True)
 with (folder/'current_run.log').open('w') as log:
  proc=subprocess.run(['/opt/homebrew/bin/mpirun','-n','4','/Users/haule/dbin/dmft2',definition],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT)
 assert proc.returncode==0,name
 assert 'DMFT2 END' in (folder/'current_run.log').read_text(),name
 records[name]={'source':src,'definition':definition,'seconds':time.monotonic()-start}
 print('DONE',name,records[name]['seconds'],flush=True)
(root/'dmft_runs.json').write_text(json.dumps(records,indent=2)+'\n')
