"""Guarded installation of the approved phase-2 increment; preserve preimages."""
from pathlib import Path
import hashlib,json,os,shutil,tempfile,datetime
root=Path('/private/tmp/edmft_charge_closure_20261108')
repo=Path('/Users/haule/projects/eDMFT_new');dbin=Path('/Users/haule/dbin')
doc=repo/'src/doc/stress/charge_closure_20261108'
records=json.loads((root/'installation_manifest.json').read_text())
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
for item in records:
    src,dst=Path(item['source']),Path(item['destination'])
    if not src.is_relative_to(root) or not (dst.is_relative_to(repo/'src/python') or dst.is_relative_to(repo/'src/doc/stress') or dst.parent==dbin):raise RuntimeError('unexpected path')
    if dst.is_symlink():raise RuntimeError('refusing unexpected destination symlink: '+str(dst))
    if digest(src)!=item['sha256'] or digest(dst)!=item['before_sha256']:raise RuntimeError('preimage/source changed: '+str(dst))
if doc.exists():raise RuntimeError('checkpoint already exists')
doc.mkdir(parents=True)
for item in records:
    dst=Path(item['destination'])
    if dst.exists():
        backup=doc/'before_installation'/dst.relative_to('/Users/haule')
        backup.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(dst,backup)
        item['backup']=str(backup)
def publish(src,dst,mode):
    dst.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.'+dst.name+'.',dir=dst.parent)
    try:
        with os.fdopen(fd,'wb') as f:f.write(src.read_bytes())
        os.chmod(tmp,mode);os.replace(tmp,dst)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)
written=[]
try:
    for item in records:
        src,dst=Path(item['source']),Path(item['destination'])
        mode=(dst.stat().st_mode&0o777) if dst.exists() else (src.stat().st_mode&0o777)
        publish(src,dst,mode);written.append(item)
        if digest(dst)!=item['sha256']:raise RuntimeError('postcopy mismatch: '+str(dst))
except Exception:
    for item in reversed(written):
        dst=Path(item['destination'])
        if item.get('backup'):
            backup=Path(item['backup']);publish(backup,dst,backup.stat().st_mode&0o777)
        else:dst.unlink()
    raise
(doc/'installation_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
record=dict(status='installed',time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),payload_files=len(records),authorization='User: Are we still in phase 2? OK to proceed; preceding approval to apply patches/install helpers.',source_root=str(repo),helpers=str(dbin),all_hashes_verified=True,native_executable_replaced=False)
(doc/'INSTALLATION.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
