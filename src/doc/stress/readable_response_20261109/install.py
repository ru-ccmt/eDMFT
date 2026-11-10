from pathlib import Path
import hashlib,json,os,shutil,tempfile
root=Path(__file__).resolve().parent
source=Path('/Users/haule/projects/eDMFT_new/src/python')
bin=Path('/Users/haule/dbin')
doc=source.parent/'doc/stress/readable_response_20261109'
files=(root/'changed_files.txt').read_text().splitlines()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for name in files:
 base=root/'base'/name
 for target in [source/name]+([] if name.startswith('tests/') else [bin/name]):
  if base.exists():assert target.is_file() and sha(target)==sha(base),'Destination changed: '+str(target)
  else:assert not target.exists(),'Destination appeared: '+str(target)
assert not doc.exists(),'Checkpoint already exists'
doc.mkdir(parents=True)
manifest=[]
# Install the new module before hooks that import it.
files.sort(key=lambda x:x!='stress_response_summary.py')
for name in files:
 staged=root/'stage'/name
 for label,target in [('source',source/name)]+([] if name.startswith('tests/') else [('dbin',bin/name)]):
  if target.exists():
   before=doc/'before'/label/name;before.parent.mkdir(parents=True,exist_ok=True)
   shutil.copy2(target,before)
  with tempfile.NamedTemporaryFile(dir=target.parent,delete=False) as stream:tmp=Path(stream.name)
  shutil.copy2(staged,tmp);os.replace(tmp,target)
  assert sha(target)==sha(staged)
  manifest.append(dict(path=str(target),sha256=sha(target)))
for name in ['REPORT.txt','component_tests.log','readable_response.patch','changed_files.txt','install.py']:
 shutil.copy2(root/name,doc/name)
(doc/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Installed',len(files),'source/test files and',len([x for x in files if not x.startswith('tests/')]),'helpers')
print(doc)
