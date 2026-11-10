"""Validate the approved installation and preserve installed-only test evidence."""
from pathlib import Path
from datetime import datetime
import hashlib
import json
import shutil

root = Path(__file__).resolve().parent
doc = Path('/Users/haule/projects/eDMFT_new/src/doc/stress/core_functional_20261109')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
manifest = json.loads((root / 'installation_manifest.json').read_text())
for item in manifest:
    dest = Path(item['destination'])
    if dest.is_symlink() or sha(dest) != item['source_sha256']:
        raise RuntimeError(f'Installed hash mismatch: {dest}')
hashes = json.loads((root / 'evidence/SHA256.json').read_text())
for name, digest in hashes.items():
    if sha(doc / name) != digest:
        raise RuntimeError(f'Installed documentation hash mismatch: {name}')
logs = ['post_install_repository.log', 'post_install_dbin.log']
for name in logs:
    text = (root / name).read_text()
    if not all(marker in text for marker in [
        'Ran 239 tests', '\nOK\n',
        'Verified installed-only imports for helpers and tests; 239 tests, zero skips.',
    ]):
        raise RuntimeError(f'Incomplete installed test evidence: {name}')
    if any(marker in text for marker in ['FAILED (', 'ERROR:', 'skipped=']):
        raise RuntimeError(f'Failed installed tests: {name}')
result = json.loads((root / 'installation_result.json').read_text())
result.update(
    status='installed and verified',
    verified_at=datetime.now().astimezone().isoformat(),
    verified_targets=len(manifest),
    verified_documentation_files=len(hashes),
    post_install_tests=dict(repository=239, dbin=239, skipped=0, installed_imports_only=True),
    stress_applied=False,
    phase=4,
    remaining='Electrostatic and explicit geometric functional terms, complete baseline subtraction, and total-stress validation.',
)
for name in logs:
    shutil.copy2(root / name, doc / name)
shutil.copy2(root / 'post_install_tests.py', doc / 'reproduce/post_install_tests.py')
shutil.copy2(root / 'record_installation.py', doc / 'reproduce/record_installation.py')
(root / 'installation_result.json').write_text(json.dumps(result, indent=2) + '\n')
(doc / 'INSTALLATION.json').write_text(json.dumps(result, indent=2) + '\n')
post_files = logs + ['reproduce/post_install_tests.py', 'reproduce/record_installation.py', 'INSTALLATION.json']
(doc / 'POST_INSTALL_SHA256.json').write_text(json.dumps(
    {name: sha(doc / name) for name in post_files}, indent=2
) + '\n')
print(f'Installed and verified: {len(manifest)} targets, {len(hashes)} evidence files; 239 tests passed in each installed location.')
