from pathlib import Path
import sys,os,json,shutil,subprocess,re
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from wien_response_density import Density
case=sys.argv[1]
report=json.loads((ROOT/(case+'_charge_check.json')).read_text())
target=report['central_band_number']
base=Path('/private/tmp/edmft_general_potential_response_20261108')/case/'joint'
env=os.environ.copy()
for key in list(env):
 if key.startswith(('EDMFT_','DFTKS_','WIEN_')):del env[key]
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
rows=[]
for step in [.001,.0005]:
 densities=[];mus=[]
 for sign in ['minus','plus']:
  src=base/(sign+'_'+str(step));dst=ROOT/case/'fixed_N'/(sign+'_'+str(step));dst.mkdir(parents=True,exist_ok=True)
  if not (dst/'done').exists():
   for p in src.iterdir():
    if p.is_file() and (p.name.startswith(case+'.') or p.name in ['EF.dat','projectorw.dat','Eorb_imp.dat','dmft2.def'] or p.name.startswith('sig.inp')):
     if p.name.startswith((case+'.output',case+'.scf')) or p.suffix=='.clmval':continue
     shutil.copy2(p,dst/p.name)
   file=dst/(case+'.indmf2');lines=file.read_text().splitlines();i=next(j for j,l in enumerate(lines) if l.strip() and not l.startswith('#'));fields=lines[i].split();fields[2]='1';lines[i]=' '.join(fields);file.write_text('\n'.join(lines)+'\n')
   in2_name=next(l for l in (dst/'dmft2.def').read_text().splitlines() if int(l.split(',')[0])==5).split("'")[1]
   file=dst/in2_name;lines=file.read_text().splitlines();fields=lines[1].split();fields[1]=f'{target:.16g}';lines[1]=' '.join(fields);file.write_text('\n'.join(lines)+'\n')
   with (dst/'dmft2.log').open('w') as log:subprocess.run(['/Users/haule/dbin/dmft2_response','dmft2.def'],cwd=dst,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
   if 'DMFT2 END' not in (dst/'dmft2.log').read_text():raise RuntimeError(str(dst))
   (dst/'done').touch()
  densities.append(Density(dst/(case+'.clmval'),dst/(case+'.struct')).vector())
  mus.append(float((dst/'EF.dat').read_text().split()[0]))
  print(case,step,sign,'done',mus[-1],flush=True)
 fd=(densities[1]-densities[0])/(2*step)
 analytic=np.load(ROOT/(case+'_joint_fixedN.npz'))['density']
 row=dict(step=step,density_relative=float(np.linalg.norm(fd-analytic)/np.linalg.norm(fd)),density_maximum=float(np.max(abs(fd-analytic))),dmu_eV=(mus[1]-mus[0])/(2*step),target_number=target)
 rows.append(row);print(row,flush=True)
(ROOT/(case+'_fixed_N_endpoints.json')).write_text(json.dumps(rows,indent=2)+'\n')
