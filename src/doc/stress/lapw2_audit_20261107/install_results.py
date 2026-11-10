"""Install only the authorized lapw2 correction, preserving prior binaries."""
from pathlib import Path
import shutil,hashlib,json,subprocess

root=Path(__file__).resolve().parent
wien=Path('/Users/haule/wien2k_24')
dest=Path('/Users/haule/projects/eDMFT_new/src/doc/stress/lapw2_audit_20261107')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert (root/'relaxed_runs.json').exists()
assert (root/'REPORT.txt').exists()
report=json.loads((root/'comparison.json').read_text())
assert set(report)=={'FeSe','SrVO3'}
assert all(all(v['basis_consistency'].values()) and v['clmval_identical'] for v in report.values())
source=wien/'SRC_lapw2/val_corr_stress.F'
assert sha(source)==sha(root/'val_corr_stress.before.F'),'Source changed during tests; do not overwrite'
paths=[source,wien/'lapw2',wien/'lapw2c',wien/'SRC_lapw2/lapw2',wien/'SRC_lapw2/lapw2c']
manifest={}
for p in paths:
    backup=p.with_name(p.name+'.before_g2_audit_20261107')
    assert not backup.exists(),f'Backup already exists: {backup}'
    shutil.copy2(p,backup)
    manifest[str(p)]={'before_sha256':sha(p),'backup':str(backup)}
shutil.copy2(root/'SRC_lapw2/val_corr_stress.F',source)
for name in ['lapw2','lapw2c']:
    for p in [wien/name,wien/'SRC_lapw2'/name]:
        shutil.copy2(root/'SRC_lapw2'/name,p)
for p in paths:manifest[str(p)]['after_sha256']=sha(p)
dest.mkdir(parents=True,exist_ok=True)
for name in ['REPORT.txt','comparison.json','initial_comparison.json','geometry_results.json',
             'geometry_probe.f90','geometry.log','build.log','runs.log','relax.log','relax_continued.log',
             'lapw2-g2-audit.patch','check_geometry.py','compare_cases.py','run_cases.py',
             'relax_cases.py','patch_source.py','install_results.py','make_report.py','finalize_checks.py',
             'val_corr_stress.before.F','FeSe_vector_manifest.json','SrVO3_vector_manifest.json',
             'binaries_used.json','source_manifest.json','compiler.txt']:
    p=root/name
    if p.exists():shutil.copy2(p,dest/name)
# Retain the actual input and output records, not just derived numbers. Large
# vector and density files remain in the isolated test directories in /tmp.
for folder in sorted((root/'runs').iterdir()):
    out=dest/'results'/folder.name
    out.mkdir(parents=True,exist_ok=True)
    for p in folder.iterdir():
        if p.is_file() and (p.suffix in ['.def','.struct','.in0','.in1','.in1c','.in2','.in2c','.inc','.inm','.scf0','.scf1','.scf2','.scfc','.scfm','.log'] or p.name.endswith('.output2')):
            shutil.copy2(p,out/p.name)
    for p in (folder/'audit').rglob('*'):
        if p.is_file() and ('scf' in p.suffix or p.suffix in ['.txt','.json'] or p.name.startswith('unused_')):
            target=out/'audit'/p.relative_to(folder/'audit')
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(p,target)
(dest/'installation.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
print('Saved audit:',dest)
