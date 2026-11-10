from pathlib import Path
import os,shutil,subprocess
root=Path(__file__).resolve().parent
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
for case in ['FeSe','SrVO3']:
 for step in [0,.01,-.01,.005,-.005]:
  folder=root/case/'hxc_core'/str(step);folder.mkdir(parents=True,exist_ok=True)
  for name in [case+'.struct',case+'.vsp',case+'.vns']:shutil.copy2(root/case/'xc'/str(step)/name,folder/name)
  for name in [case+'.inc','lcore.def']:shutil.copy2(root/case/'core'/'0'/name,folder/name)
  (folder/'core_native.bin').unlink(missing_ok=True)
  with (folder/'run.log').open('w') as log:subprocess.run([str(root/'native/SRC_lcore/lcore'),'lcore.def'],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
  print(case,step,'completed',flush=True)
