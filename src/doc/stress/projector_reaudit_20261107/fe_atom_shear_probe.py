"""Isolate the first-k raw-projector geometry check by correlated Fe atom.

Only scratch input files change. These one-k calculations are NOT material
stress estimates. The test observable precedes self-energy integration.
"""
from pathlib import Path
import os,shutil,subprocess,concurrent.futures
r=Path(__file__).parent;s=r/'FeSe'
env=os.environ|{'WIENROOT':'/Users/haule/wien2k_24','WIEN_DMFT_ROOT':'/Users/haule/dbin',
'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'}
for k in list(env):
 if k.startswith(('EDMFT_STRESS_','WIEN_STRESS_','EDMFT_PROJECTOR_')):del env[k]
def run(atom):
 d=r/f'Fe_atom{atom}';d.mkdir(exist_ok=True)
 names=['FeSe.'+x for x in ['struct','in0','in1','in1c','in2','in2c','inc','inm','klist','kgen','rotlm','vsp','vns','clmcor','r2v','indmfi','indmf2']]+['EF.dat','projectorw.dat','sig.inp','sig.inp1']
 for name in names:
  if (s/name).exists():shutil.copy2(s/name,d/name)
 lines=(s/'FeSe.indmfl').read_text().splitlines()
 firstblock=lines[15:]
 end=next(i for i,x in enumerate(firstblock[1:],1) if x.strip().startswith('2     5   5'))
 block=lines[3:8];block[0]=f'{atom} 1 -1'
 new=lines[:2]+['1 # one correlated atom for raw-projector diagnostic']+block+[lines[13],'1 5 5']+firstblock[:end]
 (d/'FeSe.indmfl').write_text('\n'.join(new)+'\n')
 kl=(s/'FeSe.klist').read_text().splitlines()
 (d/'FeSe.klist').write_text(kl[0]+'\nEND\n')
 (d/'mpi_prefix.dat').write_text('/opt/homebrew/bin/mpirun -n 1\n')
 def cmd(label,command):
  with (d/(label+'.log')).open('w') as log:
   p=subprocess.run(command,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
  if p.returncode:raise RuntimeError((atom,label,p.returncode))
 cmd('lapw1',['/Users/haule/dbin/x_dmft.py','-p','lapw1'])
 # dftKS writes serial vector names when MPI size is one.
 cmd('prepare',['/Users/haule/dbin/x_dmft.py','-d','dmft2'])
 text=(s/'FeSe.indmf2').read_text().splitlines()
 q=text[1].split();q[2]='0';text[1]=' '.join(q)
 (d/'FeSe.indmf2').write_text('\n'.join(text)+'\n')
 shutil.copy2(s/'sig.inp1',d/'sig.inp1')
 (d/'projector_strain_fd.in').touch()
 cmd('dmft2',['/opt/homebrew/bin/mpirun','-n','1','/Users/haule/dbin/dmft2','dmft2.def'])
 print('ATOM',atom,(d/'projector_strain_fd.dat').read_text(),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,[1,2]))
