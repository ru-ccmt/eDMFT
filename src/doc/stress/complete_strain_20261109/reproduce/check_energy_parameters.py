"""Read actual native searched radial parameters on the converged endpoints."""
from pathlib import Path
import sys,os,time,subprocess,shutil,json
import numpy as np
from scipy.io import FortranFile
root=Path(__file__).resolve().parent
case=sys.argv[1]
base,dirs={'FeSe':(Path('/Users/haule/projects/materials/FeSe_again/FeSe_stress'),['FeSe_RELA_strain_m001','FeSe_RELA','FeSe_RELA_strain_p001']),'SrVO3':(Path('/Users/haule/projects/materials/SrVO3'),['V_strain_m001','V_plus_force','V_strain_p001'])}[case]
exe=Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if case=='SrVO3' else 'lapw1')
env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))};env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
out=root/case/'energy_parameters';out.mkdir(exist_ok=True);rows=[]
for sign,name in zip(['minus','center','plus'],dirs):
 src=base/name;d=out/sign;d.mkdir(exist_ok=True)
 if sign=='center':d=root/case/'current'
 else:
  for suf in ['struct','in1','in1c','vsp','vns','klist','kgen']:
   p=src/(case+'.'+suf)
   if p.exists():shutil.copy2(p,d/p.name)
  shutil.copy2(root/case/'current/lapw1.def',d/'lapw1.def')
  with (d/'lapw1.log').open('w') as f:subprocess.run([str(exe),'lapw1.def'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
  if 'LAPW1 END' not in (d/'lapw1.log').read_text():raise ValueError(str(d))
 nat=2 if case=='FeSe' else 4;params={}
 with FortranFile(d/(case+'.vector')) as f:
  for a in range(nat):
   e=f.read_record('f8');elo=f.read_record('f8').reshape(-1,4).T
   for l in range(11):
    main=e[l]-200 if e[l]>150 else e[l]
    lo=elo[l,elo[l]<995] if l<4 else np.zeros(0)
    params[f'{a},{l}']=np.r_[main,lo].tolist()
 rows.append(dict(source=str(src),directory=str(d),parameters=params))
cols=base/dirs[1]/(case+'.struct');c0=np.array([float(cols.read_text().splitlines()[3][i:i+10]) for i in [0,10,20]])
scales=[]
for name in [dirs[0],dirs[2]]:
 line=(base/name/(case+'.struct')).read_text().splitlines()[3];scales.append(float(np.mean(np.array([float(line[i:i+10]) for i in [0,10,20]])/c0)))
deriv={k:((np.array(rows[2]['parameters'][k])-rows[0]['parameters'][k])/(scales[1]-scales[0])).tolist() for k in rows[0]['parameters']}
report=dict(case=case,rows=rows,strain_interval=scales[1]-scales[0],derivative=deriv)
(out/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n');print(case,json.dumps(deriv),flush=True)
