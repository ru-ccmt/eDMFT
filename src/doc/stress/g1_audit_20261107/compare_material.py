from pathlib import Path
import re,json,hashlib
import numpy as np
root=Path(__file__).resolve().parent
FAC=14710.5164057
def tensor(p,key):
 s=p.read_text();return np.array([[float(x) for x in re.findall(r'^:'+key+f'{a:03d}:'+r'\s*(.*)',s,re.M)[-1].split()[:3]] for a in [1,2,3]])
def ordinary(p):
 return [x for x in p.read_text().splitlines() if x.startswith(':') and not x.startswith((':STR',':DBG'))]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
report={}
dftfd=json.loads(Path('/Users/haule/projects/eDMFT_new/src/doc/stress/lapw2_audit_20261107/comparison.json').read_text())
for level in ['DFT','DMFT']:
 for case in ['SrVO3','FeSe']:
  folder=root/'runs'/f'{level}_{case}'
  if not (folder/'audit/corrected'/(case+'.scf2')).exists():continue
  ans={}
  for mode in ['before','radial_only','corrected']:
   f=folder/'audit'/mode
   if not (f/(case+'.scf2')).exists():continue
   val=tensor(f/(case+'.scf2'),'STR_VAL')
   den=tensor(f/(case+'.scf0'),'STR_DENS');core=tensor(f/(case+'.scfc'),'STR_CORE')
   t=(den-core-val)*FAC
   ans[mode]={'tensor_gpa':t.tolist(),'mean_gpa':float(np.trace(t)/3),'clmval_hash':digest(f/(case+'.clmval'))}
  if level=='DFT':fd=dftfd[case]['numerical_mean_gpa']
  else:
   fd=json.loads((root/(case+'_existing_fd.json')).read_text())['stress_derivative_gpa']/3
  ans['numerical_mean_gpa']=fd
  ans['difference_gpa']=ans['corrected']['mean_gpa']-fd
  ans['charge_bytes_identical']=ans['before']['clmval_hash']==ans['corrected']['clmval_hash']
  ans['ordinary_scf2_records_identical']=ordinary(folder/'audit/before'/(case+'.scf2'))==ordinary(folder/'audit/corrected'/(case+'.scf2'))
  ans['angular_delta_gpa']=(np.array(ans['corrected']['tensor_gpa'])-np.array(ans['radial_only']['tensor_gpa'])).tolist()
  ans['radial_delta_gpa']=(np.array(ans['radial_only']['tensor_gpa'])-np.array(ans['before']['tensor_gpa'])).tolist()
  report[f'{level}_{case}']=ans
(root/'material_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
for name,x in report.items():
 print(name)
 for mode in ['before','radial_only','corrected']:print(' ',mode,np.diag(x[mode]['tensor_gpa']),x[mode]['mean_gpa'])
 print(' FD',x['numerical_mean_gpa'],'difference',x['difference_gpa'],'charge equal',x['charge_bytes_identical'],'ordinary records',x['ordinary_scf2_records_identical'])
