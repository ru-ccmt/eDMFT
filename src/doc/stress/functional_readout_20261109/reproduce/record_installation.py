"""Record the approved installation and successful installed-only validation."""
from pathlib import Path
from datetime import datetime
import hashlib
import json
import shutil

root = Path(__file__).resolve().parent
doc = Path('/Users/haule/projects/eDMFT_new/src/doc/stress/functional_readout_20261109')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
manifest = json.loads((root / 'installation_manifest.json').read_text())
for item in manifest:
    assert sha(Path(item['destination'])) == item['source_sha256'], item['destination']
checks = {}
for label in ['repository', 'dbin']:
    name = 'post_install_' + label + '.log'
    text = (root / name).read_text()
    assert 'Ran 225 tests' in text and '\nOK\n' in text
    assert 'Verified all checked helpers and tests load from installed locations.' in text
    shutil.copy2(root / name, doc / name)
    checks[label] = dict(passed=225, failures=0, errors=0, skipped=0,
                         verified_installed_imports=True, log=name)
(doc / 'reproduce').mkdir(exist_ok=True)
shutil.copy2(root / 'post_install_tests.py', doc / 'reproduce/post_install_tests.py')
report = (doc / 'REPORT.txt').read_text()
old = 'This NEW\nfunctional-readout checkpoint is staged for review, not installed.'
assert old in report
report = report.replace(old,
    'The functional-readout checkpoint is installed in eDMFT_new and dbin after explicit user approval.\n'
    'Post-install validation: 225/225 tests passed in each location, zero skips;\n'
    'all checked helpers and tests imported from installed paths.', 1)
(doc / 'REPORT.txt').write_text(report)
result = json.loads((root / 'installation_result.json').read_text())
result.update(status='installed and verified',
              verified_at=datetime.now().astimezone().isoformat(),
              post_install_checks=checks,
              installed_hashes_verified=len(manifest),
              current_phase=4, lattice_scalar_readout_validated=True,
              next_step='Reconcile physical-potential cancellation and Fortran baseline; complete Klein correction before six strains and production integration',
              limitations=['Complete Klein/core/baseline reconciliation remains pending',
                           'Broader geometry/spin support remains pending',
                           'MPI/streaming remains needed for large meshes',
                           'No new production stress or pressure agreement is claimed'])
(root / 'installation_result.json').write_text(json.dumps(result, indent=2) + '\n')
(doc / 'INSTALLATION.json').write_text(json.dumps(result, indent=2) + '\n')
checksums = {str(p.relative_to(doc)): sha(p) for p in sorted(doc.rglob('*'))
             if p.is_file() and p.name != 'SHA256.json'}
(doc / 'SHA256.json').write_text(json.dumps(checksums, indent=2) + '\n')
assert all(sha(doc / p) == digest for p, digest in checksums.items())
print(f'Validated all {len(manifest)} installed targets and {len(checksums)} documentation files.')
print(doc / 'REPORT.txt')
