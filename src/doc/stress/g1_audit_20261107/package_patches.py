from pathlib import Path
import difflib,json,shutil,hashlib,subprocess
root=Path(__file__).resolve().parent;manifest=json.loads((root/'source_changes.json').read_text())
for kind in ['lapw2','dmft2']:
 patch=[];fixture=root/('patch_check_'+kind);fixture.mkdir(exist_ok=True)
 for name,item in manifest.items():
  if ('wien2k_24' in name)!=(kind=='lapw2'):continue
  p=Path(name);rel=Path('SRC_lapw2' if kind=='lapw2' else 'src/dmft2')/p.name
  before=Path(item['backup']).read_text();after=p.read_text()
  patch.extend(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='a/'+str(rel),tofile='b/'+str(rel)))
  f=fixture/rel;f.parent.mkdir(parents=True,exist_ok=True);f.write_text(before)
 p=root/(kind+'-g1-audit-20261107.patch');p.write_text(''.join(patch))
 with p.open() as f:subprocess.run(['/usr/bin/patch','--batch','-p1'],cwd=fixture,stdin=f,check=True,stdout=subprocess.PIPE)
 for name,item in manifest.items():
  if ('wien2k_24' in name)!=(kind=='lapw2'):continue
  rel=Path('SRC_lapw2' if kind=='lapw2' else 'src/dmft2')/Path(name).name
  assert hashlib.sha256((fixture/rel).read_bytes()).hexdigest()==item['after_sha256']
print('Both portable patches reproduce exact installed sources.')
