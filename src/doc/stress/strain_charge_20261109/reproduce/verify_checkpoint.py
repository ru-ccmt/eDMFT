"""Validate recorded source convergence, screened residuals and patch preimages."""
from pathlib import Path
import ast,hashlib,json
root=Path(__file__).resolve().parent
if root.name=='reproduce':root=root.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for label in ['repository','dbin']:
 log=(root/f'tests_{label}.log').read_text()
 assert 'Ran 215 tests' in log and '\nOK\n' in log
 assert 'Verified staged source/field helpers and tests, with installed dependencies.' in log
for mesh,nks in [('three_k',[3,3]),('full_bz',[75,140])]:
 for case,nk in zip(['SrVO3','FeSe'],nks):
  r=json.loads((root/mesh/case/'REPORT.json').read_text())
  assert r['metadata']['kpoints']==nk
  for sector in ['valence','core','raw','source']:
   aa,bb=[r['steps'][h]['comparisons'][sector]['all']['relative'] for h in ['0.00025','0.000125']]
   assert 3.8<aa/bb<4.2,(case,mesh,sector,aa/bb)
  assert r['richardson_source_relative']<1e-7
  assert r['richardson_PW_relative']<2e-5
  assert r['hydrostatic_field_relative']<1e-12
  assert abs(r['source_charge'])<1e-10
  sol=r['screened'];assert sol['converged'] and sol['residual_norm']<=sol['threshold']
  assert sol['relative_to_numerical_source_solution']<2e-5
  assert r['includes_window_normalization'] and not r['stress_applied']
  if mesh=='three_k':
   assert r['third_step']['richardson_relative']<1e-8
for path in (root/'stage').rglob('*.py'):ast.parse(path.read_text(),filename=str(path))
for item in json.loads((root/'installation_manifest.json').read_text()):
 source,dest=Path(item['source']),Path(item['destination'])
 assert sha(source)==item['source_sha256']
 assert (sha(dest) if dest.exists() else None)==item['before_sha256']
print('215-test suites, native source convergence, screened residuals and unchanged targets verified.')
