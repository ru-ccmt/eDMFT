"""Install only the approved, hashed patch; preserve unrelated working changes."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import shutil
import subprocess
import sys


WORK = Path(__file__).resolve().parent
PROPOSAL = WORK / 'proposal'
ROOT = Path('/Users/haule/projects/eDMFT_new')
ARCHIVE = ROOT / 'src/doc/stress/augmentation_response_20261108'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def main():
    proposal = json.loads((PROPOSAL / 'INSTALLATION_PROPOSAL.json').read_text())
    for line in (PROPOSAL / 'MANIFEST.sha256').read_text().splitlines():
        expected, relative = line.split('  ', 1)
        if sha(PROPOSAL / relative) != expected:
            raise RuntimeError('Proposal archive changed: ' + relative)
    for row in proposal['targets']:
        if sha(Path(row['target'])) != row['baseline_sha256']:
            raise RuntimeError('Target baseline changed: ' + row['target'])
        if sha(PROPOSAL / row['source']) != row['proposed_sha256']:
            raise RuntimeError('Proposed source changed: ' + row['source'])
    if ARCHIVE.exists():
        raise RuntimeError('Installation archive already exists')

    binary_paths = [ROOT / 'src/dmft2/dmft2', Path('/Users/haule/dbin/dmft2')]
    binary_hashes = {str(path): sha(path) for path in binary_paths}
    shutil.copytree(PROPOSAL, ARCHIVE)
    shutil.copy2(ARCHIVE / 'REPORT.txt', ARCHIVE / 'PROPOSAL_REPORT.txt')
    for name in ['apply_approved.py', 'verify_installed.py']:
        shutil.copy2(WORK / name, ARCHIVE / name)
    native = ROOT / 'src/dmft2/atpar'
    backup = ARCHIVE / 'build_backup'
    backup.mkdir()
    for path in native.iterdir():
        if path.is_file() and (path.suffix in ('.o', '.mod') or path.name == 'libatpar.a'):
            shutil.copy2(path, backup / path.name)
    record = dict(status='installing', authorized_by='User: ok, apply',
                  started_utc=datetime.now(timezone.utc).isoformat(),
                  targets=proposal['targets'], production_binaries_before=binary_hashes)
    status = ARCHIVE / 'INSTALLATION.json'
    status.write_text(json.dumps(record, indent=2) + '\n')
    try:
        for row in proposal['targets']:
            target = Path(row['target'])
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(PROPOSAL / row['source'], target)
            if sha(target) != row['proposed_sha256']:
                raise RuntimeError('Installed hash mismatch: ' + str(target))
        with (ARCHIVE / 'installed_build.log').open('w') as log:
            subprocess.run(['make', '-j2', 'libatpar.a',
                            'PF90=/opt/homebrew/bin/mpif90'], cwd=native,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        env = os.environ.copy()
        env.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1',
                   OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
        subprocess.run([sys.executable, str(ARCHIVE / 'verify_installed.py'),
                        '--root', str(ROOT), '--work', str(WORK / 'installed_native_check'),
                        '--evidence', str(ARCHIVE)], env=env, check=True)
        for row in proposal['targets']:
            if sha(Path(row['target'])) != row['proposed_sha256']:
                raise RuntimeError('Source changed during verification: ' + row['target'])
        if {str(path): sha(path) for path in binary_paths} != binary_hashes:
            raise RuntimeError('Production dmft2 executable changed during installation')
        record.update(status='installed_and_verified',
                      completed_utc=datetime.now(timezone.utc).isoformat(),
                      validation=json.loads((ARCHIVE / 'installed_validation.json').read_text()))
        text = (ARCHIVE / 'REPORT.txt').read_text()
        old = '''The next layer of the native response engine has been implemented and tested
in scratch. The proposed patch is NOT installed. The user's newly supplied
instruction says "Do not change my source code before asking"; installation
of this concrete patch therefore awaits approval. Repository sources, dbin
helpers, production executables, and material SCF inputs were not changed by
this checkpoint. All new builds and material probes used scratch directories.'''
        new = '''This patch was approved by the user ("ok, apply"), installed in
/Users/haule/projects/eDMFT_new, and verified. The matching Python helper is
installed in /Users/haule/dbin. The radial libatpar.a library was rebuilt.
Production dmft2 executables and material SCF inputs were not changed.
All material probes used scratch directories; the production library and
installed helpers were additionally checked after installation.'''
        if old not in text:
            raise RuntimeError('Report status paragraph changed')
        text = text.replace(old, new, 1)
        text += '''
POST-INSTALLATION VERIFICATION

All seven source/helper targets match the approved SHA256 hashes. The radial
library was rebuilt with the existing Homebrew MPI compiler and normal
O3/OpenMP flags. 79 installed-source Python tests and nine tests importing
the dbin helper pass. All 32 mixed-radial cases pass when linked to the actual
installed libatpar.a. Results and hashes are in installed_validation.json.

Existing Makefile contents are preserved in baseline/; previous library,
object and module files are in build_backup/. INSTALLATION.json records the
approved targets and verification. PROPOSAL_REPORT.txt preserves the original
pre-installation note. No changes were committed by this installation.
'''
        (ARCHIVE / 'REPORT.txt').write_text(text)
    except BaseException as error:
        record.update(status='failed', error=repr(error))
        raise
    finally:
        status.write_text(json.dumps(record, indent=2) + '\n')
    entries = [f'{sha(p)}  {p.relative_to(ARCHIVE)}'
               for p in sorted(ARCHIVE.rglob('*'))
               if p.is_file() and p.name != 'MANIFEST.sha256']
    (ARCHIVE / 'MANIFEST.sha256').write_text('\n'.join(entries) + '\n')
    print('Approved patch installed and verified. Report:', ARCHIVE / 'REPORT.txt')


if __name__ == '__main__':
    main()
