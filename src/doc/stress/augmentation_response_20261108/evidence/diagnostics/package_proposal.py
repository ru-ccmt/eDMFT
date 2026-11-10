"""Archive a reviewable patch and evidence; NEVER modify the target tree."""
from pathlib import Path
import difflib
import hashlib
import json
import shutil
import subprocess

WORK = Path(__file__).resolve().parent
ROOT = Path('/Users/haule/projects/eDMFT_new')
STAGE = WORK / 'stage'
OUT = WORK / 'proposal'

def digest(path):
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def copy(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)

def main():
    OUT.mkdir(exist_ok=False)
    paths = {
        'Makefile': 'src/dmft2/atpar/Makefile',
        'lapw_radial_potential.f90': 'src/dmft2/atpar/lapw_radial_potential.f90',
        'check_radial_potential.py': 'src/dmft2/tests/check_radial_potential.py',
        'test_radial_potential.f90': 'src/dmft2/tests/test_radial_potential.f90',
        'lapw_augmentation_response.py': 'src/python/lapw_augmentation_response.py',
        'test_lapw_augmentation_response.py': 'src/python/tests/test_lapw_augmentation_response.py',
    }
    patch, manifest = [], []
    for name, relative in paths.items():
        old, new = ROOT / relative, STAGE / name
        if name != 'Makefile' and old.exists():
            raise RuntimeError(f'new target now exists: {old}')
        copy(new, OUT / 'source' / relative)
        if old.exists():
            copy(old, OUT / 'baseline' / relative)
        patch.extend(difflib.unified_diff(
            old.read_text().splitlines(True) if old.exists() else [],
            new.read_text().splitlines(True),
            fromfile='a/' + relative if old.exists() else '/dev/null',
            tofile='b/' + relative))
        manifest.append(dict(target=str(old), baseline_sha256=digest(old),
                             proposed_sha256=digest(new), source='source/' + relative))
    dbin = Path('/Users/haule/dbin/lapw_augmentation_response.py')
    if dbin.exists():
        raise RuntimeError('unexpected existing dbin helper')
    manifest.append(dict(target=str(dbin), baseline_sha256=None,
                         proposed_sha256=digest(STAGE / dbin.name),
                         source='source/src/python/' + dbin.name))
    (OUT / 'augmentation_response.patch').write_text(''.join(patch))
    copy(WORK / 'REPORT.txt', OUT / 'REPORT.txt')
    (OUT / 'INSTALLATION_PROPOSAL.json').write_text(json.dumps(dict(
        status='prepared_not_applied', root=str(ROOT),
        authorization='Awaiting user approval under supplied AGENTS instruction.',
        targets=manifest), indent=2) + '\n')

    evidence = OUT / 'evidence'
    for name in ['python_regression.log', 'staged_build.log', 'staged_validation.json',
                 'energy_cancellation.json', 'profile.json']:
        copy(WORK / name, evidence / name)
    copy(WORK / 'native_pair_tests/report.json', evidence / 'native_checked_report.json')
    for case in ['FeSe', 'SrVO3']:
        for name in ['vns_probe/report.json', 'moving_probe/report.json',
                     'moving_probe/complete_matrices.json',
                     'moving_probe/complete_matrices.npz',
                     'moving_probe/central_response.npz']:
            copy(WORK / case / name, evidence / case / name)
    for pattern in ['*.py', '*.f90']:
        for source in WORK.glob(pattern):
            copy(source, evidence / 'diagnostics' / source.name)

    native = WORK / 'native/SRC_lapw1'
    original = Path('/Users/haule/wien2k_24/SRC_lapw1')
    native_patch = []
    for source in sorted(native.iterdir()):
        if source.suffix not in ('.F', '.f', '.f90', '.F90', '.inc'):
            continue
        base = original / source.name
        if not base.exists() or source.read_bytes() == base.read_bytes():
            continue
        native_patch.extend(difflib.unified_diff(
            base.read_text().splitlines(True), source.read_text().splitlines(True),
            fromfile='a/SRC_lapw1/' + source.name, tofile='b/SRC_lapw1/' + source.name))
    (evidence / 'diagnostics/native_probe_only.patch').write_text(''.join(native_patch))

    input_hashes = {}
    refs = Path('/private/tmp/edmft_native_operator_20261108/native_probe')
    for case in ['FeSe', 'SrVO3']:
        for source in (refs / case / 'center').iterdir():
            if source.is_file() and (source.suffix in (
                    '.struct', '.vsp', '.vns', '.in1', '.in1c', '.klist', '.def', '.vector')
                    or source.name.startswith('native.')):
                input_hashes[str(source)] = digest(source)
    for source in [
        Path('/Users/haule/projects/materials/SrVO3/V_plus_force/projectorw.dat'),
        Path('/Users/haule/projects/materials/FeSe_again/FeSe_stress/FeSe_RELA/projectorw.dat'),
        native / 'lapw1', native / 'lapw1c', WORK / 'export_radial', WORK / 'probe_spherical',
        refs / 'SRC_lapw1/lapw1', refs / 'SRC_lapw1/lapw1c']:
        input_hashes[str(source)] = digest(source)
    for source in sorted((ROOT / 'src/dmft2/atpar').glob('*.f90')):
        input_hashes[str(source)] = digest(source)
    for source in sorted(native.glob('*.F')):
        input_hashes[str(source)] = digest(source)
    (evidence / 'provenance.json').write_text(json.dumps(input_hashes, indent=2) + '\n')

    # Apply only in a private scratch tree to check the exact reviewable patch.
    check = WORK / 'patch_application_check'
    check.mkdir(exist_ok=False)
    for name, relative in paths.items():
        target = check / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if (ROOT / relative).exists():
            copy(ROOT / relative, target)
    result = subprocess.run(['patch', '-p1', '-i', str(OUT / 'augmentation_response.patch')],
                            cwd=check, check=True, capture_output=True, text=True)
    for name, relative in paths.items():
        assert digest(check / relative) == digest(STAGE / name)
    for row in manifest:
        assert digest(Path(row['target'])) == row['baseline_sha256']
    (evidence / 'patch_application_check.log').write_text(result.stdout + result.stderr)
    entries = [f'{digest(p)}  {p.relative_to(OUT)}' for p in sorted(OUT.rglob('*')) if p.is_file()]
    (OUT / 'MANIFEST.sha256').write_text('\n'.join(entries) + '\n')
    print(json.dumps(dict(package=str(OUT), files=len(entries),
                          patch_lines=len(patch), targets_unchanged=True), indent=2))

if __name__ == '__main__':
    main()
