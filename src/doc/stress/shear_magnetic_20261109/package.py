from pathlib import Path
import hashlib,json,difflib
r=Path(__file__).resolve().parent
repo=Path('/Users/haule/projects/eDMFT_new/src/python');bin_dir=Path('/Users/haule/dbin')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
items=[];diff=[]
for p in sorted((r/'stage').rglob('*.py')):
    relative=p.relative_to(r/'stage');base=r/'base'/relative
    if sha(base)==sha(p):continue
    target=repo/relative
    assert sha(target)==sha(base),f'Source changed since staging: {target}'
    if relative.parent==Path('.'):
        assert sha(bin_dir/relative.name) in (None,sha(base)),f'Installed helper changed: {relative}'
    path='src/python/'+str(relative)
    items.append(dict(path=path,base_sha256=sha(base),staged_sha256=sha(p)))
    old=base.read_text().splitlines(keepends=True) if base.is_file() else []
    diff.extend(difflib.unified_diff(old,p.read_text().splitlines(keepends=True),fromfile='a/'+path if base.exists() else '/dev/null',tofile='b/'+path))
(r/'manifest.json').write_text(json.dumps(items,indent=2)+'\n')
(r/'shear_magnetic_integration.patch').write_text(''.join(diff))
print(len(items),'changed/new source and test files')
for item in items:print(item['path'])
