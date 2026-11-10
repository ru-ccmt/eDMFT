from pathlib import Path
import hashlib,json,os,shutil,datetime
r=Path(__file__).resolve().parent
repo=Path('/Users/haule/projects/eDMFT_new');bin_dir=Path('/Users/haule/dbin')
evidence=repo/'src/doc/stress/shear_magnetic_20261109'
items=json.loads((r/'manifest.json').read_text());actions=[]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
for x in items:
    dest=repo/x['path'];rel=Path(x['path']).relative_to('src/python');src=r/'stage'/rel
    assert sha(dest)==x['base_sha256'],f'Source drift: {dest}'
    assert sha(src)==x['staged_sha256'],f'Staged drift: {src}'
    actions.append((src,dest,Path('repo')/x['path']))
    if rel.parent==Path('.'):
        dest=bin_dir/rel.name
        assert sha(dest) in (None,x['base_sha256']),f'Installed helper drift: {dest}'
        actions.append((src,dest,Path('dbin')/rel.name))
assert 'Ran 294 tests' in (r/'tests_final.log').read_text()
assert (r/'tests_final.log').read_text().rstrip().endswith('OK')
report=json.loads((r/'MnTe_total/report.json').read_text())
assert report['converged'] and report['full_tensor_within_symmetry'] and report['diagnostic_only']
assert len(report['strains'])==2
assert all(report['physical_policy'][k] for k in ['magnetic','fixed_number','normalization_enabled'])
evidence.mkdir(parents=True,exist_ok=False)
for src,dest,relative in actions:
    if dest.exists():
        back=evidence/'backup'/relative;back.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(dest,back)
records=[]
for src,dest,relative in actions:
    tmp=dest.with_name(dest.name+'.install_shear_magnetic.tmp');assert not tmp.exists()
    mode=dest.stat().st_mode if dest.exists() else src.stat().st_mode
    shutil.copy2(src,tmp);os.chmod(tmp,mode);os.replace(tmp,dest)
    assert sha(src)==sha(dest)
    records.append(dict(path=str(dest),sha256=sha(dest)))
for p in r.iterdir():
    if p.is_file() and p.suffix in ('.json','.py','.patch','.txt','.log'):
        shutil.copy2(p,evidence/p.name)
for relative in ['MnTe_total/report.json','native_MnTe/preparation.json','native_Ferro/preparation.json','native_Ferro/fermif.log']:
    dest=evidence/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(r/relative,dest)
record=dict(installed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),authorization='User explicitly requested installation and implementation of general shear and automatic magnetic response.',files=records)
(evidence/'INSTALL.json').write_text(json.dumps(record,indent=2)+'\n')
print('Installed',len(items),'source/test files and',len(actions)-len(items),'helpers. Evidence:',evidence)
