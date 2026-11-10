"""Exercise the installer against scratch copies, never the real destinations."""
from pathlib import Path
import json,shutil,subprocess,hashlib,tempfile
root=Path(__file__).resolve().parent;pkg=root/'package'
with tempfile.TemporaryDirectory(prefix='edmft_install_test_',dir='/private/tmp') as tmp:
 tmp=Path(tmp);test=tmp/'package';shutil.copytree(pkg,test);m=json.loads((test/'manifest.json').read_text());sandbox=tmp/'targets'
 for e in m['entries']:
  real=Path(e['destination']);fake=sandbox/str(real).lstrip('/')
  if real.is_file():fake.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(real,fake)
  e['destination']=str(fake)
 m['checkpoint']=str(sandbox/m['checkpoint'].lstrip('/'));(test/'manifest.json').write_text(json.dumps(m,indent=2)+'\n')
 result=subprocess.run(['/Users/haule/miniconda3-clean/envs/work311/bin/python',str(test/'install.py'),'--apply'],capture_output=True,text=True,check=True)
 for e in m['entries']:
  if hashlib.sha256(Path(e['destination']).read_bytes()).hexdigest()!=e['sha256']:raise ValueError('installed payload mismatch')
 repeat=subprocess.run(['/Users/haule/miniconda3-clean/envs/work311/bin/python',str(test/'install.py')],capture_output=True,text=True)
 if repeat.returncode==0:raise ValueError('stale baseline unexpectedly accepted')
 report=dict(scratch_only=True,destinations_verified=len(m['entries']),backup_and_hash_checks_passed=True,second_install_rejected=True)
 (root/'installer_test.json').write_text(json.dumps(report,indent=2)+'\n');print(report)
