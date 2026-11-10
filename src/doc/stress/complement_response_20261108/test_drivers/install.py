"""Install the validated optional complement module with guards and snapshots."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,datetime

scratch=Path('/private/tmp/edmft_complement_response_20261108')
repo=Path('/Users/haule/projects/eDMFT_new');dbin=Path('/Users/haule/dbin')
doc=repo/'src/doc/stress/complement_response_20261108'
python='/Users/haule/miniconda3-clean/envs/work311/bin/python'

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
expected=json.loads((scratch/'expected.json').read_text())
for name,old in expected.items():
    p=Path(name)
    if (digest(p) if p.is_file() else None)!=old:raise RuntimeError('changed since staging: '+name)
if doc.exists():raise RuntimeError('checkpoint already exists')
items=[(scratch/'stage/lapw_complement_response.py',repo/'src/python/lapw_complement_response.py'),
       (scratch/'stage/lapw_complement_response.py',dbin/'lapw_complement_response.py'),
       (scratch/'stage/test_lapw_complement_response.py',repo/'src/python/tests/test_lapw_complement_response.py'),
       (scratch/'stage/stress_response_roadmap_20261108.txt',repo/'src/doc/stress/stress_response_roadmap_20261108.txt')]
for source,target in items:
    if not source.is_file() or str(target) not in expected:raise RuntimeError('unguarded target '+str(target))
manifest=dict(created=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=[],validation={},
              production_binaries_before={str(dbin/n):digest(dbin/n) for n in ['dmft2','dmft2_response']},
              scope='optional static empty-complement response; no production stress switch')
doc.mkdir(parents=True)
def save():(doc/'INSTALLATION.json').write_text(json.dumps(manifest,indent=2)+'\n')
save()
for source,target in items:
    rel=Path('repo')/target.relative_to(repo) if target.is_relative_to(repo) else Path('dbin')/target.name
    old=digest(target) if target.exists() else None
    if old:
        backup=doc/'backup'/rel;backup.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(target,backup)
    snapshot=doc/'source'/rel;snapshot.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,snapshot)
    temp=target.with_name(target.name+'.complement_install_'+str(os.getpid()))
    if temp.exists():raise RuntimeError('unexpected temporary installation file')
    shutil.copy2(source,temp);os.replace(temp,target)
    manifest['files'].append(dict(target=str(target),before=old,after=digest(target),source=str(snapshot.relative_to(doc))))
    save()
for name in ['REPORT.txt','staged_tests.log']:shutil.copy2(scratch/name,doc/name)
env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
for label,path,pattern in [('source',repo/'src/python','test_*.py'),('dbin',dbin,'test_lapw*')]:
    command=[python,'-m','unittest','discover','-s',str(repo/'src/python/tests'),'-p',pattern,'-q']
    result=subprocess.run(command,cwd=scratch,env=env|dict(PYTHONPATH=str(path),EDMFT_RADIAL_LIBRARY=str(dbin/'liblapw_response.so')),
                          text=True,capture_output=True)
    out=result.stdout+result.stderr;(doc/(label+'_tests.log')).write_text(out)
    manifest['validation'][label]=out.strip();save()
    if result.returncode:raise RuntimeError('installed tests failed: '+label)
results=doc/'results';results.mkdir()
for name in ['FeSe_complement.json','SrVO3_complement.json','FeSe_reduced_material.json','SrVO3_reduced_material.json']:
    shutil.copy2(scratch/name,results/name)
drivers=doc/'test_drivers';drivers.mkdir()
for name in ['check_material.py','build_reduced_material.py','check_reduced_material.py','install.py']:
    shutil.copy2(scratch/name,drivers/name)
fingerprints={}
old=Path('/private/tmp/edmft_coupled_valence_response_20261108')
for name in ['FeSe_b509_one','SrVO3_b703_one']:
    for p in (old/name/'central').glob('native*.bin'):
        fingerprints[str(p)]=dict(size=p.stat().st_size,sha256=digest(p))
    p=old/(name+'_coupled.npz')
    fingerprints[str(p)]=dict(size=p.stat().st_size,sha256=digest(p))
for name in ['build_material.py']:
    p=old/name;fingerprints[str(p)]=dict(size=p.stat().st_size,sha256=digest(p))
(doc/'MATERIAL_FINGERPRINTS.json').write_text(json.dumps(fingerprints,indent=2)+'\n')
for name,sha in manifest['production_binaries_before'].items():
    if digest(Path(name))!=sha:raise RuntimeError('production binary changed: '+name)
manifest['production_binaries_unchanged']=True;manifest['complete']=True;save()
print(json.dumps(dict(checkpoint=str(doc),validation=manifest['validation'],production_binaries_unchanged=True),indent=2))
