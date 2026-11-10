from pathlib import Path
import json,hashlib,shutil
root=Path(__file__).resolve().parent
manifest=json.loads((root/'source_changes.json').read_text())
for name,item in manifest.items():
 p=Path(name)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==item['before_sha256'],f'Source changed during audit: {p}'
 backup=p.with_name(p.name+'.before_g1_audit_20261107')
 assert not backup.exists(),backup
 shutil.copy2(p,backup);shutil.copy2(item['prepared'],p)
 item['installed_backup']=str(backup)
for p in [Path('/Users/haule/wien2k_24')/n for n in ['lapw2','lapw2c','SRC_lapw2/lapw2','SRC_lapw2/lapw2c']]+[Path('/Users/haule/dbin/dmft2'),Path('/Users/haule/projects/eDMFT_new/src/dmft2/dmft2')]:
 backup=p.with_name(p.name+'.before_g1_audit_20261107')
 assert not backup.exists(),backup
 shutil.copy2(p,backup)
manifest['binary_backups']=[str(p) for p in [Path('/Users/haule/wien2k_24/lapw2.before_g1_audit_20261107'),Path('/Users/haule/wien2k_24/lapw2c.before_g1_audit_20261107'),Path('/Users/haule/dbin/dmft2.before_g1_audit_20261107')]]
(root/'installation.json').write_text(json.dumps(manifest,indent=2)+'\n')
