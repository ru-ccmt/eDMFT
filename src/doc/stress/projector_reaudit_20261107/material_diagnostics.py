from pathlib import Path
import os,shutil,subprocess,time,json,hashlib,concurrent.futures
root=Path(__file__).parent
env=os.environ|{'WIENROOT':'/Users/haule/wien2k_24','WIEN_DMFT_ROOT':'/Users/haule/dbin',
 'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'}
for key in list(env):
 if key.startswith(('EDMFT_STRESS_','WIEN_STRESS_','EDMFT_PROJECTOR_')):del env[key]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(spec):
 case,src=spec;src=Path(src);d=root/case;d.mkdir(exist_ok=True)
 sufs=['struct','clmsum','clmup','clmdn','in0','in1','in1c','in2','in2c','inc','inm','klist','kgen','rotlm','vsp','vns','clmcor','r2v','indmfl','indmfi','indmf2']
 names=[case+'.'+s for s in sufs]+['EF.dat','Eorb_imp.dat','projectorw.dat','SOlapm.dat','sig.inp']
 names+=sorted(p.name for p in src.glob('sig.inp[0-9]*'))
 hashes={}
 for name in names:
  p=src/name
  if p.exists():shutil.copy2(p,d/name);hashes[name]=sha(p)
 (d/'mpi_prefix.dat').write_text('/opt/homebrew/bin/mpirun -n 4\n')
 def cmd(label,command,e=env):
  print('START',case,label,flush=True);start=time.monotonic()
  with (d/(label+'.log')).open('w') as log:r=subprocess.run(command,cwd=d,env=e,stdout=log,stderr=subprocess.STDOUT)
  print('END',case,label,r.returncode,round(time.monotonic()-start,2),flush=True)
  if r.returncode:raise RuntimeError((case,label,r.returncode))
 cmd('lapw1',['/Users/haule/dbin/x_dmft.py','-p','lapw1'])
 cmd('prepare',['/Users/haule/dbin/x_dmft.py','-d','-p','dmft2'])
 shutil.copy2(src/(case+'.indmf2'),d/(case+'.indmf2'))
 # Preserve precisely the saved sparse self-energy after wrapper preparation.
 for name in names:
  if name.startswith('sig.inp'):shutil.copy2(src/name,d/name)
 (d/'projector_strain_fd.in').touch()
 (d/'normalization_strain_fd.in').touch()
 diag=env|{'EDMFT_STRESS_PROJECTOR_DIAG':'1','EDMFT_STRESS_DETAILED_DIAGNOSTICS':'1'}
 for mode in ['DIAG','APPLY']:
  cmd(mode,['/opt/homebrew/bin/mpirun','-n','4','/Users/haule/dbin/dmft2','dmft2.def'],diag|{'EDMFT_STRESS_NORM_SYMMETRY':mode})
  out=d/mode;out.mkdir(exist_ok=True)
  for name in [case+'.scf2',case+'.outputdmf2','projector_strain_fd.dat','normalization_strain_derivative.dat','normalization_moment.dat','dmft2_info.out']:
   if (d/name).exists():shutil.copy2(d/name,out/name)
 result={'source':str(src),'source_unchanged':all(sha(src/n)==h for n,h in hashes.items()),'input_hashes':hashes}
 (d/'provenance.json').write_text(json.dumps(result,indent=2)+'\n')
 return case,result
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 results=dict(pool.map(run,[('SrVO3','/Users/haule/projects/materials/SrVO3/V_plus_force'),('FeSe','/Users/haule/projects/materials/FeSe_again/FeSe_stress/FeSe_RELA')]))
(root/'materials_provenance.json').write_text(json.dumps(results,indent=2)+'\n')
