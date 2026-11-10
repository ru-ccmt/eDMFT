"""Apply only after approval: verify every preimage, back up, copy atomically."""
from pathlib import Path
import hashlib,json,shutil,datetime,os
root=Path(__file__).resolve().parent
rows=json.loads((root/'installation_manifest.json').read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
for path,value in json.loads((root/'expected_before_installation.json').read_text()).items():
 if sha(Path(path))!=value:raise RuntimeError('Dependency changed since review: '+path)
for row in rows:
 if sha(Path(row['source']))!=row['new_sha256']:raise RuntimeError('Staged source changed: '+row['source'])
 if sha(Path(row['target']))!=row['old_sha256']:raise RuntimeError('Destination changed: '+row['target'])
tag=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
checkpoint=Path('/Users/haule/projects/eDMFT_new/src/doc/stress/xc_core_response_20261108')
for row in rows:
 source,target=Path(row['source']),Path(row['target']);target.parent.mkdir(parents=True,exist_ok=True)
 if target.exists():
  backup=checkpoint/'before_installation'/target.relative_to('/');backup.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(target,backup);row['backup']=str(backup)
 temporary=target.with_name(target.name+'.xc_core_'+tag+'.tmp');shutil.copy2(source,temporary);os.replace(temporary,target)
 if sha(target)!=row['new_sha256']:raise RuntimeError('Verification failed: '+str(target))
record=dict(time=tag,files=rows,production_native_executables_changed=False)
(checkpoint/'INSTALLATION.json').write_text(json.dumps(record,indent=2)+'\n')
(root/'INSTALLATION.json').write_text(json.dumps(record,indent=2)+'\n')
print('Installed and verified',len(rows),'files. No native production executable changed.')
