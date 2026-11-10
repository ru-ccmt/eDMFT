"""Guarded, backed-up installation of the coupled spherical-input milestone."""

from pathlib import Path
import datetime
import hashlib
import json
import os
import shutil
import subprocess

scratch = Path('/private/tmp/edmft_coupled_valence_response_20261108')
repo = Path('/Users/haule/projects/eDMFT_new')
dbin = Path('/Users/haule/dbin')
doc = repo / 'src/doc/stress/coupled_valence_response_20261108'
python = '/Users/haule/miniconda3-clean/envs/work311/bin/python'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


expected = json.loads((scratch / 'install_expected.json').read_text())
for name, old in expected.items():
    path = Path(name)
    actual = digest(path) if path.is_file() else None
    if actual != old:
        raise RuntimeError('changed since staging: ' + name)
if doc.exists():
    raise RuntimeError('checkpoint already exists')

items = []
for name in ['lapw_spherical_response.py', 'lapw_nonspherical_response.py',
             'lapw_projector_response.py', 'lapw_spherical_potential_response.py',
             'lapw_radial_native.py']:
    items.extend([(scratch/'stage'/name, repo/'src/python'/name),
                  (scratch/'stage'/name, dbin/name)])
for name in ['test_lapw_potential_assembly.py', 'test_lapw_radial_native.py']:
    items.append((scratch/'stage'/name, repo/'src/python/tests'/name))
items.extend([
    (scratch/'stage/lapw_radial_capi.f90', repo/'src/dmft2/atpar/lapw_radial_capi.f90'),
    (scratch/'stage/stress_response_roadmap_20261108.txt',
     repo/'src/doc/stress/stress_response_roadmap_20261108.txt'),
    (scratch/'clean_build/liblapw_response.so', repo/'src/dmft2/atpar/liblapw_response.so'),
    (scratch/'clean_build/liblapw_response.so', dbin/'liblapw_response.so'),
])
for source, target in items:
    if not source.is_file() or str(target) not in expected:
        raise RuntimeError('unguarded installation item: ' + str(target))
for name in ['outwin.f90', 'outwin_response.f90', 'lapw_radial_response.f90']:
    if digest(scratch/'clean_build'/name) != digest(repo/'src/dmft2/atpar'/name):
        raise RuntimeError('native dependency changed: ' + name)

manifest = dict(
    created=datetime.datetime.now(datetime.timezone.utc).isoformat(), files=[],
    validation={}, scope='optional coupled spherical-potential valence response',
    production_binaries_before={str(dbin/name): digest(dbin/name)
                                for name in ['dmft2', 'dmft2_response']},
)
doc.mkdir(parents=True)


def save():
    (doc/'INSTALLATION.json').write_text(json.dumps(manifest, indent=2)+'\n')


def atomic_copy(source, target):
    # Replacing the pathname leaves any running process's old shared-library
    # mapping intact. Never truncate a potentially mapped library in place.
    temp = target.with_name(target.name + '.coupled_install_' + str(os.getpid()))
    if temp.exists():
        raise RuntimeError('unexpected temporary installation file')
    shutil.copy2(source, temp)
    os.replace(temp, target)


save()
for source, target in items:
    rel = (Path('repo')/target.relative_to(repo) if target.is_relative_to(repo)
           else Path('dbin')/target.name)
    old = digest(target) if target.exists() else None
    if old:
        backup = doc/'backup'/rel
        backup.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(target, backup)
    snapshot = doc/'source'/rel
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, snapshot)
    atomic_copy(source, target)
    manifest['files'].append(dict(target=str(target), before=old,
                                 after=digest(target), source=str(snapshot.relative_to(doc))))
    save()

for name in ['REPORT.txt', 'clean_build.log', 'staged_tests.log']:
    shutil.copy2(scratch/name, doc/name)
env = os.environ.copy()
env.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1',
           OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
for label, location, pattern, library in [
    ('source', repo/'src/python', 'test_*.py', repo/'src/dmft2/atpar/liblapw_response.so'),
    ('dbin', dbin, 'test_lapw*', dbin/'liblapw_response.so'),
]:
    command = [python, '-m', 'unittest', 'discover', '-s',
               str(repo/'src/python/tests'), '-p', pattern, '-q']
    result = subprocess.run(command, cwd=scratch, text=True, capture_output=True,
                            env=env | dict(PYTHONPATH=str(location),
                                           EDMFT_RADIAL_LIBRARY=str(library)))
    output = result.stdout + result.stderr
    (doc/(label+'_tests.log')).write_text(output)
    manifest['validation'][label] = output.strip()
    save()
    if result.returncode:
        raise RuntimeError('installed test failure: '+label)

results = doc/'results'
results.mkdir()
for pattern in ['*_coupled.json', '*_projector.json', '*_spherical.json']:
    for path in scratch.glob(pattern):
        shutil.copy2(path, results/path.name)
drivers = doc/'test_drivers'
drivers.mkdir()
for name in ['build_material.py', 'setup_endpoints.py', 'check_coupled.py',
             'check_projector.py', 'check_spherical.py', 'install.py']:
    shutil.copy2(scratch/name, drivers/name)
for path in [Path('/private/tmp/edmft_augmentation_response_20261108/material_augmentation.py'),
             Path('/private/tmp/edmft_radial_all_l_20261108/test_wien_matrices.py')]:
    shutil.copy2(path, drivers/path.name)

fingerprints = {}
for case in ['FeSe', 'SrVO3']:
    for folder in scratch.glob(case+'*/central'):
        if 'emax200_control' in str(folder):
            continue
        for path in folder.iterdir():
            if path.is_file() and (path.name in [case+'.struct', case+'.in1', case+'.in1c',
                    case+'.in2', case+'.in2c', case+'.vector', case+'.vsp', case+'.vns',
                    case+'.indmfl', 'projectorw.dat', 'SOlapm.dat', 'EF.dat', 'dmft2.def']
                    or path.name.startswith(('response', 'sig.inp'))):
                fingerprints[str(path)] = dict(size=path.stat().st_size, sha256=digest(path))
for path in [dbin/'dmft2_response',
             Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1/lapw1'),
             Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1/lapw1c')]:
    fingerprints[str(path)] = dict(size=path.stat().st_size, sha256=digest(path))
(doc/'MATERIAL_FINGERPRINTS.json').write_text(json.dumps(fingerprints, indent=2)+'\n')
for name, old in manifest['production_binaries_before'].items():
    if digest(Path(name)) != old:
        raise RuntimeError('production binary changed: '+name)
manifest['production_binaries_unchanged'] = True
manifest['complete'] = True
save()
print(json.dumps(dict(checkpoint=str(doc), validation=manifest['validation'],
                     production_binaries_unchanged=True), indent=2))
