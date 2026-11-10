"""Independent N(mu) root followed by native dmft2; no analytic dmu input."""
from pathlib import Path
import sys,os,shutil,subprocess,json
import numpy as np
from scipy.optimize import brentq
root=Path(__file__).resolve().parent/'full_bz'
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from frozen_sigma_response import FrozenSigmaResponse
from number_reference import NumberReference
case=sys.argv[1];nr={'FeSe':66,'SrVO3':43}[case]
central=root/case/'charge/0';base,_=FrozenSigmaResponse.from_export(central/'response')
base.eps=base.eps[:,:nr];base.uc=base.uc[:,:nr]
target=NumberReference(base)(base.meta.mu)
env=os.environ.copy()
for key in list(env):
 if key.startswith(('EDMFT_','DFTKS_','WIEN_')):del env[key]
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
report=dict(target_band_number=target,steps={})
for step in [.001,-.001,.0005,-.0005]:
 src=root/case/'charge'/str(step);dst=root/case/'fixed_n'/str(step);dst.mkdir(parents=True,exist_ok=True)
 m,_=FrozenSigmaResponse.from_export(src/'response');conversion=float((src/'EF.dat').read_text().split()[0])/m.meta.mu
 m.eps=m.eps[:,:nr];m.uc=m.uc[:,:nr]
 oracle=NumberReference(m)
 def residual(mu):
  return oracle(mu)-target
 if (dst/'done').exists():
  oldmu=float((dst/'EF.dat').read_text().split()[0])/conversion
  if abs(residual(oldmu))>5e-12:raise ValueError('completed endpoint has wrong target')
  report['steps'][step]=dict(mu_eV=oldmu*conversion,number_residual=residual(oldmu),reused_dense_reference=True)
  print(case,step,report['steps'][step],flush=True)
  continue
 mu=brentq(residual,base.meta.mu-.02,base.meta.mu+.02,xtol=5e-15,rtol=1e-14)
 for p in src.iterdir():
  if p.is_file() and (p.name.startswith(case+'.') or p.name in ['EF.dat','projectorw.dat','Eorb_imp.dat','dmft2.def'] or p.name.startswith('sig.inp')):
   if p.name.startswith((case+'.output',case+'.scf')) or p.suffix=='.clmval':continue
   shutil.copy2(p,dst/p.name)
 p=dst/(case+'.indmf2');lines=p.read_text().splitlines();i=next(j for j,l in enumerate(lines) if l.strip() and not l.startswith('#'));fields=lines[i].split();fields[1]=f'{mu*conversion:.16g}';fields[2]='0';lines[i]=' '.join(fields);p.write_text('\n'.join(lines)+'\n');(dst/'EF.dat').write_text(f'{mu*conversion:.16g}\n')
 with (dst/'dmft2.log').open('w') as log:subprocess.run(['/Users/haule/dbin/dmft2_response','dmft2.def'],cwd=dst,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 if 'DMFT2 END' not in (dst/'dmft2.log').read_text():raise RuntimeError(str(dst))
 (dst/'done').touch()
 row=dict(mu_eV=mu*conversion,number_residual=residual(mu));report['steps'][step]=row;print(case,step,row,flush=True)
(root/(case+'_fixed_n_roots.json')).write_text(json.dumps(report,indent=2)+'\n')
