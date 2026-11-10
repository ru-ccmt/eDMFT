from pathlib import Path
import re,json,hashlib
import numpy as np
root=Path(__file__).resolve().parent;prior=Path('/private/tmp/edmft_g1_audit_20261107');fac=14710.5164057
prev=json.loads((prior/'material_comparison.json').read_text())
def tensor(path,key):
 text=path.read_text()
 return np.array([[float(x) for x in re.findall('^:'+key+f'{a:03d}:'+r'\s*(.*)',text,re.M)[-1].split()[:3]] for a in [1,2,3]])
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def records(path):return [s for s in path.read_text().splitlines() if s.startswith(':') and not s.startswith((':STR',':DBG'))]
out={}
for level in ['DFT','DMFT']:
 for case in ['SrVO3','FeSe']:
  label=f'{level}_{case}';p=root/'runs'/label/'audit/current';old=prior/'runs'/label/'audit/corrected'
  stress=(tensor(p/(case+'.scf0'),'STR_DENS')-tensor(p/(case+'.scfc'),'STR_CORE')-tensor(p/(case+'.scf2'),'STR_VAL'))*fac
  oldstress=np.array(prev[label]['corrected']['tensor_gpa']);fd=prev[label]['numerical_mean_gpa']
  assert np.isfinite(stress).all()
  out[label]={'tensor_gpa':stress.tolist(),'mean_gpa':float(np.trace(stress)/3),'existing_fd_mean_gpa':fd,'analytic_minus_fd_gpa':float(np.trace(stress)/3-fd),'max_change_from_first_pass_gpa':float(abs(stress-oldstress).max()),'charge_identical_to_first_pass':digest(p/(case+'.clmval'))==digest(old/(case+'.clmval')),'ordinary_scf2_records_identical':records(p/(case+'.scf2'))==records(old/(case+'.scf2'))}
  assert out[label]['max_change_from_first_pass_gpa']<1e-8
  assert out[label]['charge_identical_to_first_pass'] and out[label]['ordinary_scf2_records_identical']
(root/'material_comparison.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
