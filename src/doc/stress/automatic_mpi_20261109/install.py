from pathlib import Path
import difflib, hashlib, json, os, shutil, stat, tempfile
root=Path(__file__).resolve().parent
source=Path('/Users/haule/projects/eDMFT_new/src/python')
bin=Path('/Users/haule/dbin')
doc=source.parent/'doc/stress/automatic_mpi_20261109'
files=['response_mpi.py','dmft2_step.py','x_dmft.py','tests/test_analytic_stress_integration.py','tests/test_response_mpi.py']
helpers=files[:3]
sha=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
# Verify ALL destinations before modifying ANY destination.
for name in files:
    base=root/'base'/name
    targets=[source/name]+([bin/name] if name in helpers else [])
    for target in targets:
        if base.exists():
            assert target.is_file() and sha(target)==sha(base), 'Changed since staging: '+str(target)
        else:
            assert not target.exists(), 'New destination appeared: '+str(target)
configs=[]
for case,folder,n in [('SrVO3','/Users/haule/projects/materials/SrVO3/V_plus_force',4),('FeSe','/Users/haule/projects/materials/FeSe_again/FeSe_stress/FeSe_RELA',6)]:
    path=Path(folder)/(case+'.stress_response.json')
    expected={'solver':{'backend':'analytic','mpi':['/Users/haule/miniconda3-clean/envs/work311/bin/mpiexec','-genv','FI_PROVIDER','tcp','-n',str(n)],'native_mpi':['/opt/homebrew/bin/mpirun','-n',str(n)]}}
    assert json.loads(path.read_text())==expected,'Config changed: '+str(path)
    configs.append(path)
assert not doc.exists(), 'Checkpoint already exists'
doc.mkdir(parents=True)
manifest=[]
patches=[]
for name in files:
    base=root/'base'/name
    staged=root/'stage'/name
    patches.extend(difflib.unified_diff(base.read_text().splitlines(True) if base.exists() else [],staged.read_text().splitlines(True),fromfile='a/src/python/'+name if base.exists() else '/dev/null',tofile='b/src/python/'+name))
    for label,target in [('source',source/name)]+([('dbin',bin/name)] if name in helpers else []):
        if target.exists():
            backup=doc/'before'/label/name
            backup.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(target,backup)
        target.parent.mkdir(parents=True,exist_ok=True)
        # Atomic replacement prevents another launcher seeing a partial module.
        with tempfile.NamedTemporaryFile(dir=target.parent,delete=False) as stream:
            tmp=Path(stream.name)
        shutil.copy2(staged,tmp)
        os.replace(tmp,target)
        assert sha(target)==sha(staged)
        manifest.append({'path':str(target),'sha256':sha(target)})
for path in configs:
    backup=doc/'before/cases'/path.name
    backup.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(path,backup)
    cfg=json.loads(path.read_text())
    del cfg['solver']['mpi']
    del cfg['solver']['native_mpi']
    with tempfile.NamedTemporaryFile(mode='w',dir=path.parent,delete=False) as stream:
        tmp=Path(stream.name)
        json.dump(cfg,stream,indent=2)
        stream.write('\n')
    tmp.chmod(stat.S_IMODE(path.stat().st_mode))
    os.replace(tmp,path)
    manifest.append({'path':str(path),'sha256':sha(path)})
(doc/'automatic_mpi.patch').write_text(''.join(patches))
(doc/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
for name in ['REPORT.txt','component_tests.log','real_mpi_tests.json','validate_launches.py','install.py']:
    shutil.copy2(root/name,doc/name)
print('Installed 5 source/test files and 3 helpers; migrated 2 case MPI configs.')
print(doc)
