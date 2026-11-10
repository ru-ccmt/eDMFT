"""Guarded installation of the reviewed explicit LAPW band-geometry patch. Preview is the default.

Run with --apply only after user approval. All source/evidence hashes and all
current destination preimages are checked before any destination is written.
Preimages are copied to scratch backups. No git operation or deletion is used.
"""
from pathlib import Path
from tempfile import NamedTemporaryFile
from datetime import datetime
import argparse
import hashlib
import json
import os
import shutil

root = Path(__file__).resolve().parent
repo = Path('/Users/haule/projects/eDMFT_new')
dbin = Path('/Users/haule/dbin')
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--apply', action='store_true')
args = parser.parse_args()
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
items = json.loads((root / 'installation_manifest.json').read_text())
evidence = root / 'evidence'
evidence_hashes = json.loads((evidence / 'SHA256.json').read_text())
doc = repo / 'src/doc/stress/band_geometry_20261108'
for rel, digest in evidence_hashes.items():
    p = evidence / rel
    if sha(p) != digest:
        raise RuntimeError(f'Changed evidence: {p}')
for item in items:
    source, dest = Path(item['source']), Path(item['destination'])
    if not (source.is_relative_to(root) and
            (dest.is_relative_to(repo) or dest.is_relative_to(dbin))):
        raise RuntimeError(f'Unexpected path: {item}')
    if sha(source) != item['source_sha256']:
        raise RuntimeError(f'Changed staged file: {source}')
    before = sha(dest) if dest.exists() else None
    if dest.is_symlink() or before != item['before_sha256']:
        raise RuntimeError(f'Changed destination preimage: {dest}')
if doc.exists():
    raise RuntimeError(f'Documentation target already exists; inspect before replacing: {doc}')
print(f'Preflight passed: {len(items)} source/helper/roadmap files and {len(evidence_hashes)} evidence files.')
if not args.apply:
    print('Preview only: nothing installed. --apply requires user approval.')
    raise SystemExit(0)
backup = root / ('before_install_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
backup.mkdir()
for i, item in enumerate(items):
    p = Path(item['destination'])
    if p.exists():
        shutil.copy2(p, backup / f'{i:02d}_{p.name}')
(backup / 'manifest.json').write_text(json.dumps(items, indent=2) + '\n')
installed = []
for item in items:
    source, dest = Path(item['source']), Path(item['destination'])
    dest.parent.mkdir(parents=True, exist_ok=True)
    mode = (dest.stat().st_mode & 0o777) if dest.exists() else 0o644
    with NamedTemporaryFile(dir=dest.parent, prefix='.' + dest.name + '.', delete=False) as f:
        tmp = Path(f.name)
        f.write(source.read_bytes())
        f.flush()
        os.fsync(f.fileno())
    try:
        tmp.chmod(mode)
        os.replace(tmp, dest)
    finally:
        tmp.unlink(missing_ok=True)
    if sha(dest) != item['source_sha256']:
        raise RuntimeError(f'Installed hash mismatch: {dest}; backup is {backup}')
    installed.append(str(dest))
shutil.copytree(evidence, doc)
log = dict(installed_at=datetime.now().astimezone().isoformat(),
           backup=str(backup), destinations=installed,
           evidence=str(doc), status='installed; post-install checks still required')
(root / 'installation_result.json').write_text(json.dumps(log, indent=2) + '\n')
(doc / 'INSTALLATION.json').write_text(json.dumps(log, indent=2) + '\n')
print(f'Installed {len(installed)} files. Preserved preimages in {backup}.')
print('Run post-install repository and dbin suites before reporting completion.')
