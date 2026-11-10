"""Native one-step endpoints with common input density and fixed radial energies."""
from pathlib import Path
import os,sys,shutil,subprocess
root=Path(__file__).resolve().parent
case=sys.argv[1];source=Path('/private/tmp/edmft_coupled_valence_response_20261108')/dict(FeSe='FeSe_b509_one',SrVO3='SrVO3_b703_one')[case]/'central'
env=os.environ.copy()
for key in list(env):
 if key.startswith(('EDMFT_','DFTKS_','WIEN_')):del env[key]
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
exe=Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if case=='SrVO3' else 'lapw1')
for step in [0,.01,-.01,.005,-.005,.001,-.001,.0005,-.0005]:
 folder=root/case/'charge'/str(step);folder.mkdir(parents=True,exist_ok=True)
 if (folder/'done').exists():continue
 for p in source.iterdir():
  if p.is_file() and (p.name.startswith(case+'.') or p.name in ['EF.dat','Eorb_imp.dat','projectorw.dat','dmft2.def','lapw1.def'] or p.name.startswith('sig.inp')):
   if p.suffix in ['.vector','.clmval','.output1','.output2'] or p.name.startswith(case+'.scf'):continue
   shutil.copy2(p,folder/p.name)
 for name in [case+'.vsp',case+'.vns',case+'.clmsum']:shutil.copy2(root/case/'xc'/str(step)/name,folder/name)
 for name,cmd,extra in [('lapw1',[str(exe),'lapw1.def'],{'WIEN_LAPW1_MATRIX_DUMP':str(folder/'native')}),('dmft2',['/Users/haule/dbin/dmft2_response','dmft2.def'],{'EDMFT_RESPONSE_EXPORT':str(folder/'response')})]:
  with (folder/(name+'.log')).open('w') as log:subprocess.run(cmd,cwd=folder,env=env|extra,stdout=log,stderr=subprocess.STDOUT,check=True)
  if name=='dmft2' and 'DMFT2 END' not in (folder/(name+'.log')).read_text():raise RuntimeError(str(folder))
 (folder/'done').touch();print(case,step,'done',flush=True)
