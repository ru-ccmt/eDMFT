"""Check recorded derivative convergence, tests and unchanged install targets."""
from pathlib import Path
import ast
import hashlib
import json

root=Path(__file__).resolve().parent
if root.name=='reproduce':root=root.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for label in ['repository','dbin']:
 log=(root/f'tests_{label}.log').read_text()
 assert 'Ran 209 tests' in log and '\nOK\n' in log
 assert 'Verified staged geometry helper/test and installed dependencies.' in log
for case in ['SrVO3','FeSe']:
 report=json.loads((root/case/'REPORT.json').read_text())
 assert report['kpoints']==3
 a,b=[report['steps'][step] for step in ['0.00025','0.000125']]
 for sector in ['H_RR','H_QR','O_RR','O_QR','raw']:
  aa,bb=[item['relative'][sector]['error'] for item in [a,b]]
  assert 3.8<aa/bb<4.2,(case,sector,aa/bb)
  assert bb<8e-7,(case,sector,bb)
 assert report['richardson_joint_relative']<1e-9
 assert abs(report['richardson_projector_scalar_error_Ry'])<1e-10
 assert b['central_H_relative']<1e-8 and b['central_O_relative']<1e-11
 assert not b['stress_applied']
for path in (root/'stage').rglob('*.py'):ast.parse(path.read_text(),filename=str(path))
for item in json.loads((root/'installation_manifest.json').read_text()):
 source,dest=Path(item['source']),Path(item['destination'])
 assert sha(source)==item['source_sha256']
 assert (sha(dest) if dest.exists() else None)==item['before_sha256']
print('209-test suites, native matrix/projector convergence and unchanged targets verified.')
