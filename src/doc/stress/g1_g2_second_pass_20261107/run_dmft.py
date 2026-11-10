from pathlib import Path
import subprocess,shutil,os,time,json
root=Path(__file__).resolve().parent
old=Path('/Users/haule/dbin/dmft2.before_g1_audit_20261107');new=Path('/Users/haule/projects/eDMFT_new/src/dmft2/dmft2')
env=dict(os.environ,OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1',WIENROOT='/Users/haule/wien2k_24',WIEN_DMFT_ROOT='/Users/haule/dbin')
for k in list(env):
 if k.startswith(('WIEN_STRESS_','EDMFT_STRESS_')):del env[k]
for case,src in [('SrVO3','SrVO3_clean'),('FeSe','FeSe_final')]:
 source=Path('/private/tmp/edmft_stress_fix_20261107/tests')/src
 folder=root/'runs'/('DMFT_'+case);shutil.copytree(source,folder,dirs_exist_ok=True)
 audit=folder/'audit';audit.mkdir(exist_ok=True)
 for mode,program,extra in [('current',Path('/Users/haule/dbin/dmft2'),{})]:
  # Restore the same input generation (including EF and normalization) for every mode.
  for initial in source.iterdir():
   if initial.is_file():shutil.copy2(initial,folder/initial.name)
  start=time.monotonic();print('START',case,mode,flush=True)
  with (folder/(mode+'.log')).open('w') as f:
   p=subprocess.run(['/opt/homebrew/bin/mpirun','-n','4',str(program),'dmft2.def'],cwd=folder,env=env|extra,stdout=f,stderr=subprocess.STDOUT)
  assert p.returncode==0,(case,mode,p.returncode)
  log=(folder/(mode+'.log')).read_text()
  assert 'DMFT2 END' in log,(case,mode,log[-1000:])
  out=audit/mode;out.mkdir(exist_ok=True)
  for f in folder.iterdir():
   if f.is_file() and (f.suffix in ['.scf0','.scfc','.scf2','.outputdmf2','.clmval'] or f.name.startswith(('dmft2_info','SOlapm','Eorb_imp','EF.dat'))):shutil.copy2(f,out/f.name)
  print('DONE',case,mode,round(time.monotonic()-start,2),flush=True)
print('ALL DMFT TESTS DONE',flush=True)
