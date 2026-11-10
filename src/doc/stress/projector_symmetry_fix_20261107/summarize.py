from pathlib import Path
import json,re,hashlib
import numpy as np
r=Path(__file__).parent;old=Path('/private/tmp/projector_reaudit_20261107')
fac=14710.507641078

def tensor(p,tag):
 s=p.read_text();a=[]
 for i in range(1,4):
  row=re.findall('^'+re.escape(tag)+f'{i:03d}'+r':?\s+([^\n]+)',s,re.M)[-1]
  a.append([float(t) for t in row.split()[:3]])
 return np.array(a)

result={}
for case in ['SrVO3','FeSe']:
 d=r/case;src=Path(json.loads((old/case/'provenance.json').read_text())['source'])
 new=tensor(d/(case+'.scf2'),':STR_VAL')
 baseline=tensor(old/case/'DIAG'/(case+'.scf2'),':STR_VAL')
 previous_apply=tensor(old/case/'APPLY'/(case+'.scf2'),':STR_VAL')
 dc=tensor(src/(case+'.scf0'),':STR_DENS')-tensor(src/(case+'.scfc'),':STR_CORE')
 flags={}
 for tag in [':SUM',':VSUM',':YSUM',':DRHO']:
  rows=lambda p:[s for s in p.read_text().splitlines() if s.startswith(tag)]
  flags[tag]=rows(d/(case+'.scf2'))==rows(old/case/'DIAG'/(case+'.scf2'))
 raw={}
 for line in (d/'projector_strain_fd.dat').read_text().splitlines():
  if line.startswith('#'):continue
  a=line.split();raw[a[0]]=float(a[4])
 result[case]={
  'mean_stress_GPa':float(np.trace(dc-new)*fac/3),
  'shift_from_previous_default_GPa':float(np.trace(baseline-new)*fac/3),
  'difference_from_old_APPLY_switch_GPa':float(np.trace(previous_apply-new)*fac/3),
  'scalar_outputs_unchanged':flags,'raw_projector_relative_errors':raw}
for atom in [1,2]:
 raw={}
 for line in (r/f'Fe_atom{atom}/projector_strain_fd.dat').read_text().splitlines():
  if line.startswith('#'):continue
  a=line.split();raw[a[0]]=float(a[4])
 result[f'Fe_atom{atom}']=raw
result['unit']=json.loads((r/'unit/results.json').read_text())
result['timings']=json.loads((r/'runs.json').read_text())
result['source_and_binary_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
 for p in [Path('/Users/haule/projects/eDMFT_new/src/dmft2')/name for name in
 ['symmetrize.f90','l2main.f90','l2mainSubs.f90','lomain.f90','Makefile','dmft2']]}
(r/'results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({c:result[c] for c in ['SrVO3','FeSe','Fe_atom1','Fe_atom2']},indent=2))
