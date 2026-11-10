from pathlib import Path
import subprocess,shutil,json,hashlib
r=Path(__file__).resolve().parent;d=r/'installer_fixture';d.mkdir(exist_ok=True);m=json.loads((r/'source_manifest_before.json').read_text())
for name,v in m.items():
 rel=Path(name).relative_to('/Users/haule/wien2k_24');p=d/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(v['backup'],p)
p=d/'SRC_lapw0/Makefile';s=p.read_text().replace('FC = gfortran','FC = SITE_COMPILER');p.write_text(s)
for args in [[],['--apply'],[]]:subprocess.run(['python3',str(r/'apply_remaining_stress.py'),str(d)]+args,check=True)
assert 'FC = SITE_COMPILER' in p.read_text()
print('PASS: check, apply, idempotent check, and site compiler preservation')
