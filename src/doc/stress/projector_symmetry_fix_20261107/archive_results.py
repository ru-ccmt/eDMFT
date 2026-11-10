from pathlib import Path
import difflib
import hashlib
import json
import shutil

run = Path('/private/tmp/projector_symmetry_fix_20261107')
repo = Path('/Users/haule/projects/eDMFT_new')
out = run / 'archive'
out.mkdir(exist_ok=True)

def copy(src, rel):
    dst = out / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)

for name in ['REPORT.txt', 'results.json', 'runs.json', 'installation.json',
             'dft_control.json', 'run_material_checks.py', 'dft_control.py',
             'summarize.py', 'build.log', 'existing_checks.log', 'unit.log',
             'archive_results.py']:
    copy(run / name, name)
for sub, names in {
    'unit': ['results.json', 'symmetry_probe.f90'],
    'harness': ['results.json', 'source_hashes.json', 'rank_deficient_C.json',
                'extracted_routines.f90'],
}.items():
    for name in names:
        copy(run / sub / name, Path(sub) / name)

changes = {
    'src/dmft2/symmetrize.f90': 'symmetrize.f90',
    'src/dmft2/l2main.f90': 'l2main.f90',
    'src/dmft2/l2mainSubs.f90': 'l2mainSubs.f90',
    'src/dmft2/lomain.f90': 'lomain.f90',
    'src/dmft2/Makefile': 'Makefile',
    'src/dmft2/tests/README.txt': 'tests_README.txt',
    'src/dmft2/tests/check_projector_symmetry.py': None,
}
patch = []
for rel, old_name in changes.items():
    current = repo / rel
    old = (run / 'before' / old_name).read_text().splitlines(True) if old_name else []
    new = current.read_text().splitlines(True)
    patch.extend(difflib.unified_diff(old, new,
        fromfile=('a/' + rel if old_name else '/dev/null'), tofile='b/' + rel))
    copy(current, Path('source') / Path(rel).relative_to('src/dmft2'))
(out / 'projector_symmetry_incremental.patch').write_text(''.join(patch))

for name in ['SrVO3', 'FeSe', 'Fe_atom1', 'Fe_atom2',
             'DFT_control_before', 'DFT_control_after']:
    case = 'SrVO3' if name == 'SrVO3' else 'FeSe'
    filenames = [case + '.scf2', case + '.outputdmf2', 'new.log',
                 'projector_strain_fd.dat', 'dmft2.def', case + '.struct',
                 case + '.indmfl', case + '.in1c', case + '.in2c',
                 'control.log', 'EF.dat']
    for filename in filenames:
        src = run / name / filename
        if src.is_file():
            copy(src, Path(name) / filename)

manifest = {}
for p in sorted(out.rglob('*')):
    if p.is_file() and p.name != 'MANIFEST.json':
        manifest[str(p.relative_to(out))] = {
            'bytes': p.stat().st_size,
            'sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
        }
(out / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(f'Archived {len(manifest)} files, {sum(x["bytes"] for x in manifest.values())} bytes')
