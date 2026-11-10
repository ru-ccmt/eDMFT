from pathlib import Path
import shutil,subprocess,os,time,re,json
root=Path(__file__).resolve().parent
src=Path('/Users/haule/projects/materials/SrVO3/V_plus_force')
dst=root/'SrVO3_bands';dst.mkdir(exist_ok=True)
base=root/'SrVO3_DMFT'
for p in base.iterdir():
 if p.is_file() and p.name.startswith('SrVO3.') and p.stat().st_size<20000000: shutil.copy2(p,dst/p.name)
for name in ['SrVO3.klist','SrVO3.indmf2','SrVO3.indmfl','SrVO3.indmfi','SrVO3.in1c','SrVO3.in2c','SrVO3.kgen','SOlapm.dat','projectorw.dat','sig.inp','sig.inp1','EF.dat','Eorb_imp.dat','mpi_prefix.dat','lapw1.def','dmft2.def']:
 if (src/name).exists():shutil.copy2(src/name,dst/name)
env=dict(os.environ,OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1',WIENROOT='/Users/haule/wien2k_24',WIEN_DMFT_ROOT='/Users/haule/dbin')
for k in list(env):
 if k.startswith(('WIEN_STRESS_','EDMFT_STRESS_')):del env[k]
orig=(src/'SrVO3.in1c').read_text()
for nb in [400,680]:
 (dst/'SrVO3.in1c').write_text(re.sub(r'43\s+FIXED_NBANDS',f'{nb}   FIXED_NBANDS',orig))
 for name in ['SOlapm.dat','SrVO3.indmf2','EF.dat','sig.inp1']:shutil.copy2(src/name,dst/name)
 save=dst/f'results_{nb}';save.mkdir(exist_ok=True)
 for step in ['lapw1','dmft2']:
  cmd=['/opt/homebrew/bin/mpirun','-n','4','/Users/haule/dbin/'+('lapw1c' if step=='lapw1' else step),step+'.def']
  t=time.monotonic();print('START',nb,step,flush=True)
  with (save/(step+'.log')).open('w') as f:p=subprocess.run(cmd,cwd=dst,env=env,stdout=f,stderr=subprocess.STDOUT)
  print('END',nb,step,p.returncode,round(time.monotonic()-t,1),flush=True)
  if p.returncode:raise RuntimeError(cmd)
 for p in dst.iterdir():
  if p.is_file() and (p.name in ['SrVO3.scf2','SrVO3.scf1','SOlapm.dat','EF.dat','SrVO3.in1c'] or p.name.startswith(('dmft2_info','SrVO3.energy_'))):shutil.copy2(p,save/p.name)
print('ALL DONE',flush=True)
