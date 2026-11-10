from pathlib import Path
import hashlib,json,os,shutil,datetime
r=Path(__file__).resolve().parent
repo=Path('/Users/haule/projects/eDMFT_new');bin_dir=Path('/Users/haule/dbin')
evidence=repo/'src/doc/stress/general_cells_magnetic_20261109'
items=json.loads((r/'manifest.json').read_text());actions=[]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
for x in items:
 dest=repo/x['path'];rel=Path(x['path']).relative_to('src/python');src=r/'stage'/rel
 assert sha(dest)==x['base_sha256'],f'Source drift: {dest}'
 assert sha(src)==x['staged_sha256'],f'Staged drift: {src}'
 actions.append((src,dest,Path('repo')/x['path']))
 if rel.parent==Path('.'):
  dest=bin_dir/rel.name
  assert sha(dest) in (None,x['base_sha256']),f'Installed helper drift: {dest}'
  actions.append((src,dest,Path('dbin')/rel.name))
evidence.mkdir(parents=True,exist_ok=False)
for src,dest,relative in actions:
 if dest.exists():
  back=evidence/'backup'/relative;back.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(dest,back)
records=[]
for src,dest,relative in actions:
 temp=dest.with_name(dest.name+'.install_general_cells.tmp')
 assert not temp.exists()
 mode=dest.stat().st_mode if dest.exists() else src.stat().st_mode
 shutil.copy2(src,temp);os.chmod(temp,mode);os.replace(temp,dest)
 assert sha(src)==sha(dest)
 records.append(dict(path=str(dest),sha256=sha(dest)))
for p in r.iterdir():
 if p.is_file() and (p.suffix in ('.json','.py','.patch') or p.name in ('REPORT.txt','regression_final.log')):
  shutil.copy2(p,evidence/p.name)
for rel in ['SrVO3_regression/report.json']:
 dest=evidence/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(r/rel,dest)
record=dict(installed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),authorization='User: Please install.',files=records)
(evidence/'INSTALL.json').write_text(json.dumps(record,indent=2)+'\n')
(evidence/'INSTALL.txt').write_text('Installed the approved 35-source-file patch and 25 Python helpers.\nThe original REPORT.txt documents the pre-install review; INSTALL.json records the completed installation.\nNative snapshots remain in '+str(r/'native')+'.\n')
print('Installed',len(items),'source/test files and',len(actions)-len(items),'helpers. Evidence:',evidence)
