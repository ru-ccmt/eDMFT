"""Independent nonlinear N(mu) root, then native dmft2 at that fixed mu.

The ordinary native EF solver permits an absolute/relative Ry tolerance of
1e-9, which can dominate division by a small finite-difference amplitude.
This diagnostic never uses the analytic response to select endpoint mu.
"""
from pathlib import Path
import sys,json,shutil,os,subprocess
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from frozen_sigma_response import FrozenSigmaResponse
from wien_response_density import Density
ROOT=Path(__file__).resolve().parent;case=sys.argv[1]
base=Path('/private/tmp/edmft_general_potential_response_20261108')/case/'joint'
target=json.loads((ROOT/(case+'_charge_check.json')).read_text())['central_band_number']
retained={'FeSe':66,'SrVO3':43}[case]
env=os.environ.copy()
for key in list(env):
 if key.startswith(('EDMFT_','DFTKS_','WIEN_')):del env[key]
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
rows=[]
for step in [.001,.0005]:
 ds=[];mus=[];residuals=[]
 for sign in ['minus','plus']:
  src=base/(sign+'_'+str(step));old=ROOT/case/'fixed_N'/(sign+'_'+str(step));dst=ROOT/case/'fixed_N_polished'/(sign+'_'+str(step));dst.mkdir(parents=True,exist_ok=True)
  m,records=FrozenSigmaResponse.from_export(src/'response')
  conversion=float((src/'EF.dat').read_text().split()[0])/m.meta.mu
  # The discarded high states are thermally empty; this probe keeps all
  # occupied states and the entire window, and verifies N(mu) directly.
  m.eps=m.eps[:,:retained];m.uc=m.uc[:,:retained]
  def residual(mu):
   m.meta.mu=mu
   return float(np.trace(m.density(),axis1=1,axis2=2).sum().real-target)
  initial=float((old/'EF.dat').read_text().split()[0])/conversion
  mu=brentq(residual,initial-1e-6,initial+1e-6,xtol=5e-15,rtol=1e-14)
  resid=residual(mu);residuals.append(resid);mus.append(mu*conversion)
  if not (dst/'done').exists():
   for p in old.iterdir():
    if p.is_file() and (p.name.startswith(case+'.') or p.name in ['EF.dat','projectorw.dat','Eorb_imp.dat','dmft2.def'] or p.name.startswith('sig.inp')):
     if p.name.startswith((case+'.output',case+'.scf')) or p.suffix=='.clmval':continue
     shutil.copy2(p,dst/p.name)
   file=dst/(case+'.indmf2');lines=file.read_text().splitlines();i=next(j for j,l in enumerate(lines) if l.strip() and not l.startswith('#'));fields=lines[i].split();fields[1]=f'{mu*conversion:.16g}';fields[2]='0';lines[i]=' '.join(fields);file.write_text('\n'.join(lines)+'\n')
   (dst/'EF.dat').write_text(f'{mu*conversion:.16g}\n')
   with (dst/'dmft2.log').open('w') as log:subprocess.run(['/Users/haule/dbin/dmft2_response','dmft2.def'],cwd=dst,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
   if 'DMFT2 END' not in (dst/'dmft2.log').read_text():raise RuntimeError(str(dst))
   (dst/'done').touch()
  ds.append(Density(dst/(case+'.clmval'),dst/(case+'.struct')).vector())
  print(case,step,sign,'root residual',resid,'mu correction eV',(mu-initial)*conversion,flush=True)
 fd=(ds[1]-ds[0])/(2*step);analytic=np.load(ROOT/(case+'_joint_fixedN.npz'))['density']
 row=dict(step=step,density_relative=float(np.linalg.norm(fd-analytic)/np.linalg.norm(fd)),density_maximum=float(np.max(abs(fd-analytic))),dmu_eV=(mus[1]-mus[0])/(2*step),root_residuals=residuals)
 rows.append(row);print(row,flush=True)
(ROOT/(case+'_fixed_N_polished.json')).write_text(json.dumps(rows,indent=2)+'\n')
