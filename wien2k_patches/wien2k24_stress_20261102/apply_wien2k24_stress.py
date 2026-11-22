#!/usr/bin/env python3
"""Install exact, reviewed WIEN2k-24 source transitions, preserving site flags.

Every Fortran file is recognized by its complete SHA-256 digest. We assemble
and verify the entire prospective result in a temporary directory before any
WIENROOT file is changed. This also supports mixtures of the documented old
releases/supplements: all paths must independently have a known provenance.
Makefiles are deliberately not hash-pinned to a compiler installation. Only
two object-list tokens are inserted, with all other bytes preserved.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import uuid

BUNDLE = Path(__file__).resolve().parent


class BundleError(RuntimeError):
    pass


def sha(data):
    return hashlib.sha256(data).hexdigest() if data is not None else None


def read_regular(path):
    if path.is_symlink():
        raise BundleError(f"refusing to replace a source symlink: {path}")
    if not path.exists():
        return None
    if not path.is_file():
        raise BundleError(f"expected a regular source file: {path}")
    return path.read_bytes()


def safe_relative(value):
    path = Path(value)
    if path.is_absolute() or '..' in path.parts or not path.parts:
        raise BundleError(f"invalid bundle path: {value}")
    return path


def verify_integrity():
    """Detect stale/damaged package contents before interpreting the registry."""
    manifest = BUNDLE / 'MANIFEST.sha256'
    if not manifest.is_file():
        raise BundleError('missing bundle integrity manifest')
    checked = set()
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        digest, relative = line.split(None, 1)
        path = BUNDLE / safe_relative(relative.strip())
        if not path.is_file() or sha(path.read_bytes()) != digest:
            raise BundleError(f"bundle checksum mismatch: {relative}")
        checked.add(str(path.relative_to(BUNDLE)))
    required = {'sources.json', 'apply_wien2k24_stress.py'}
    if not required <= checked:
        raise BundleError('integrity manifest omits required installer files')
    return checked


def add_object(data, rule):
    """Modify exactly one OBJS assignment, never compiler/library settings.

    A matching token in a comment or a dependency rule is insufficient. The
    object must be in the actual continued OBJS assignment. Refuse ambiguous
    Makefiles instead of appending an assignment that might not affect linking.
    """
    text = data.decode('utf-8')
    lines = text.splitlines(keepends=True)
    pattern = re.compile(r'^\s*' + re.escape(rule['variable']) + r'\s*(?::|\+|\?)?=')
    starts = [i for i, line in enumerate(lines) if pattern.match(line)]
    if len(starts) != 1:
        raise BundleError(f"{rule['path']}: expected one {rule['variable']} assignment")
    begin = end = starts[0]
    while lines[end].rstrip('\r\n').rstrip().endswith('\\'):
        end += 1
        if end >= len(lines):
            raise BundleError(f"{rule['path']}: unterminated object list")
    block = ''.join(lines[begin:end + 1])
    token = re.compile(r'(?<![\w.])' + re.escape(rule['object']) + r'(?![\w.])')
    active = '\n'.join(line.split('#', 1)[0] for line in lines[begin:end + 1])
    count = len(token.findall(active))
    if count == 1:
        return data
    if count:
        raise BundleError(f"{rule['path']}: duplicate {rule['object']}")
    anchor = re.compile(r'(?<![\w.])' + re.escape(rule['anchor']) + r'(?![\w.])')
    if len(anchor.findall(active)) != 1:
        raise BundleError(f"{rule['path']}: missing/ambiguous anchor {rule['anchor']}")
    changed = False
    for i in range(begin, end + 1):
        body, sep, comment = lines[i].partition('#')
        if anchor.search(body):
            lines[i] = anchor.sub(rule['anchor'] + ' ' + rule['object'], body, count=1) + sep + comment
            changed = True
            break
    assert changed
    return ''.join(lines).encode('utf-8')


def prepare(root, registry, with_fftw=False):
    """Return fully verified candidate bytes without writing to WIENROOT."""
    plan = []
    selected = {p: item for p, item in registry['sources'].items()
                if with_fftw or not item['optional']}
    with tempfile.TemporaryDirectory(prefix='wien-stress-preflight-') as directory:
        stage = Path(directory)
        for relative, item in sorted(selected.items()):
            path = root / safe_relative(relative)
            before = read_regular(path)
            old_hash = sha(before)
            if old_hash == item['target_sha256']:
                continue
            matches = [x for x in item['variants'] if x['sha256'] == old_hash]
            if len(matches) != 1:
                raise BundleError(f"unknown source version: {relative} (sha256={old_hash}); "
                                  'preserve and review this file before applying')
            variant = matches[0]
            target = stage / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            if before is not None:
                target.write_bytes(before)
            patch = BUNDLE / safe_relative(variant['patch'])
            result = subprocess.run(
                ['patch', '--batch', '--forward', '-F', '0', '-p1', '-d', str(stage), '-i', str(patch)],
                capture_output=True, text=True)
            if result.returncode or not target.is_file() or sha(target.read_bytes()) != item['target_sha256']:
                raise BundleError(f"preflight failed for {relative}:\n{result.stdout}{result.stderr}")
            plan.append({'path': relative, 'before': before, 'after': target.read_bytes(),
                         'from': variant['labels']})
    for rule in registry['makefiles']:
        path = root / safe_relative(rule['path'])
        before = read_regular(path)
        if before is None:
            raise BundleError(f"missing build file: {path}; use a full WIEN2k source tree")
        after = add_object(before, rule)
        if after != before:
            plan.append({'path': rule['path'], 'before': before, 'after': after,
                         'from': ['site-preserving object-list insertion']})
    return plan


def prepare_driver(root, variants):
    """Patch the active driver, or the script template in an unpacked source tree."""
    relative = 'run_lapw'
    active = root / relative
    template = root / 'SRC/run_lapw'
    if active.is_symlink():
        if active.resolve() != template.resolve() or not template.is_file():
            raise BundleError(f'unsupported run_lapw symlink: {active}')
        relative = 'SRC/run_lapw'
    elif not active.exists():
        relative = 'SRC/run_lapw'
    before = read_regular(root / relative)
    if before is None:
        raise BundleError('--with-run-lapw requires run_lapw or SRC/run_lapw')
    digest = sha(before)
    if any(digest == v['target_sha256'] for v in variants):
        return []
    for v in variants:
        if digest == v['sha256']:
            with tempfile.TemporaryDirectory(prefix='wien-driver-') as temp:
                path = Path(temp) / 'run_lapw'
                path.write_bytes(before)
                apply_diff(path.parent, v['patch'])
                after = path.read_bytes()
            if sha(after) != v['target_sha256']:
                raise BundleError(f'unexpected patched {relative} digest')
            return [dict(path=relative, before=before, after=after, **{'from': [v['label']]})]
    raise BundleError(f'unknown {relative} (sha256={digest}); preserve site edits and review the driver patch')


def apply_diff(directory, relative):
    subprocess.run(['patch', '--batch', '--forward', '-F', '0', '-p1', '-d', str(directory),
                    '-i', str(BUNDLE / safe_relative(relative))], check=True)


def prepare_worker_sources(root, recipe, plan):
    """Check and instrument only temporary copies of the final ordinary sources."""
    pending = {item['path']: item['after'] for item in plan}
    exports = {}
    with tempfile.TemporaryDirectory(prefix='wien-export-check-') as temp:
        stage = Path(temp)
        for name, item in recipe['sources'].items():
            relative = safe_relative(name)
            before = pending[name] if name in pending else read_regular(root / relative)
            if sha(before) != item['base_sha256']:
                raise BundleError(f'unknown response-worker base: {name} (sha256={sha(before)})')
            path = stage / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(before)
        apply_diff(stage, recipe['patch'])
        for name, item in recipe['sources'].items():
            data = (stage / name).read_bytes()
            if sha(data) != item['target_sha256']:
                raise BundleError(f'unexpected response-worker source after patch: {name}')
            exports[name] = data
    return exports


def install_workers(root, destination, exports):
    """Compile four serial WIEN workers separately from the ordinary WIEN build."""
    products = [('SRC_lapw0', 'seq', 'lapw0'), ('SRC_lcore', 'lcore', 'lcore'),
                ('SRC_lapw1', 'real', 'lapw1'), ('SRC_lapw1', 'complex', 'lapw1c')]
    ignore = shutil.ignore_patterns('*.o', '*.mod', '*.smod', '*_tmp_.*', '.parallel',
                                    '.sequential', '.real', '.complex', '*.dSYM',
                                    'lapw0', 'lcore', 'lapw1', 'lapw1c')
    with tempfile.TemporaryDirectory(prefix='wien-response-build-') as temp:
        stage = Path(temp)
        for component in ('SRC_Globals', 'SRC_lapw0', 'SRC_lcore', 'SRC_lapw1'):
            shutil.copytree(root / component, stage / component, symlinks=True, ignore=ignore)
        # WIEN Makefiles may link site-provided libraries through -L../SRC_lib.
        if (root / 'SRC_lib').is_dir():
            (stage / 'SRC_lib').symlink_to(root / 'SRC_lib', target_is_directory=True)
        for name, data in exports.items():
            path = stage / name
            if path.is_symlink():
                raise BundleError(f'unexpected response source symlink: {name}')
            path.write_bytes(data)
        plan = []
        for component, target, executable in products:
            subprocess.run(['make', '-j1', '-C', str(stage / component),
                            'DESTDIR=./', target], check=True)
            product = stage / component / executable
            if not product.is_file() or not os.access(product, os.X_OK):
                raise BundleError(f'worker build did not produce {executable}')
            name = executable + '_response'
            after = product.read_bytes()
            before = read_regular(destination / name)
            if before != after or not os.access(destination / name, os.X_OK):
                plan.append({'path': name, 'before': before, 'after': after, 'install_mode': 0o755,
                             'from': ['dedicated WIEN response build']})
        if plan:
            destination.mkdir(parents=True, exist_ok=True)
            install(destination, plan)
    print(f'Four WIEN response workers verified in {destination}.')


def install(root, plan):
    """Back up a preflighted transaction and roll back any failed replacement."""
    stamp = dt.datetime.now().strftime('%Y%m%d_%H%M%S') + '-' + uuid.uuid4().hex[:8]
    parent = Path(os.environ.get('WIEN_STRESS_BACKUP_DIR', str(root)))
    backup = parent / ('stress_patch_backup_' + stamp)
    # Check for edits made while patch preflight was running.
    for item in plan:
        if read_regular(root / item['path']) != item['before']:
            raise BundleError('source changed during preflight: ' + item['path'])
    backup.mkdir(parents=True, exist_ok=False)
    records = []
    for item in plan:
        path = root / item['path']
        item['mode'] = item.get('install_mode', stat.S_IMODE(path.stat().st_mode) if path.exists() else 0o644)
        if item['before'] is not None:
            saved = backup / item['path']
            saved.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, saved)
        records.append({'path': item['path'], 'before_sha256': sha(item['before']),
                        'after_sha256': sha(item['after']), 'mode': item['mode'], 'from': item['from']})
    (backup / 'transaction.json').write_text(json.dumps(records, indent=2) + '\n')
    written = []
    try:
        for item in plan:
            path = root / item['path']
            path.parent.mkdir(parents=True, exist_ok=True)
            fd, tmp = tempfile.mkstemp(prefix='.wien-stress-', dir=path.parent)
            try:
                with os.fdopen(fd, 'wb') as stream:
                    stream.write(item['after'])
                    stream.flush()
                    os.fsync(stream.fileno())
                os.chmod(tmp, item['mode'])
                os.replace(tmp, path)
            finally:
                if os.path.exists(tmp):
                    os.unlink(tmp)
            written.append(item)
        for item in plan:
            if read_regular(root / item['path']) != item['after']:
                raise BundleError('installed content changed: ' + item['path'])
    except BaseException:
        for item in reversed(written):
            path = root / item['path']
            if item['before'] is None:
                path.unlink()
            else:
                shutil.copy2(backup / item['path'], path)
        raise
    print(f'Backup and transaction record: {backup}')
    return backup


def build(root, with_fftw, with_lapw0_mpi=False):
    # Serial make avoids upstream incomplete module dependencies; this is
    # independent of MPI execution of the resulting programs.
    steps = [('SRC_lapw0', []), ('SRC_lapw1', ['real']), ('SRC_lapw1', ['complex']),
             ('SRC_lcore', []), ('SRC_mixer', []), ('SRC_lapw2', ['real']), ('SRC_lapw2', ['complex'])]
    if with_lapw0_mpi:
        steps.insert(1, ('SRC_lapw0', ['para']))
    if with_fftw:
        steps += [('SRC_3ddens', []), ('SRC_hf', ['real']), ('SRC_hf', ['complex'])]
    for component, targets in steps:
        subprocess.run(['make', '-j1', '-C', str(root / component), *targets], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--check', action='store_true', help='preflight only (default)')
    modes.add_argument('--apply', action='store_true', help='apply all recognized source transitions')
    modes.add_argument('--verify', action='store_true', help='require exact current sources and object lists')
    parser.add_argument('--with-run-lapw', action='store_true', help='include the final screened DFT driver hook')
    parser.add_argument('--install-response-workers', type=Path, metavar='DBIN',
                        help='build four isolated response workers into DBIN; --check/--verify only preflight their sources')
    parser.add_argument('--with-debug', action='store_true', help='compatibility option; dormant diagnostics are included')
    parser.add_argument('--with-fftw-compat', action='store_true', help='also patch optional 3ddens/HF FFTW C interfaces')
    parser.add_argument('--with-lapw0-mpi', action='store_true', help='also build/install lapw0_mpi with the configured MPI compiler')
    parser.add_argument('--build', action='store_true', help='build affected programs with existing site settings')
    parser.add_argument('--install', action='store_true', help='build and install executables into WIENROOT')
    parser.add_argument('--install-helpers', type=Path, metavar='DIR')
    parser.add_argument('wienroot', type=Path)
    options = parser.parse_args()
    checked = verify_integrity()
    registry = json.loads((BUNDLE / 'sources.json').read_text())
    if registry.get('schema') != 1:
        raise BundleError('unsupported source registry schema')
    for item in registry['sources'].values():
        for variant in item['variants']:
            if variant['patch'] not in checked:
                raise BundleError('source transition missing from integrity manifest')
    for extra in [registry['response_workers'], *registry['run_lapw']]:
        if extra['patch'] not in checked:
            raise BundleError('driver or response-worker patch missing from integrity manifest')
    root = options.wienroot.expanduser().resolve()
    if not root.is_dir():
        raise BundleError(f'not a source directory: {root}')
    plan = prepare(root, registry, options.with_fftw_compat)
    if options.with_run_lapw:
        plan.extend(prepare_driver(root, registry['run_lapw']))
    if options.install_response_workers:
        exports = prepare_worker_sources(root, registry['response_workers'], plan)
    count = sum(options.with_fftw_compat or not x['optional'] for x in registry['sources'].values())
    print(f"{'upgradeable' if plan else 'applied'} {root}: {count} source files; {len(plan)} changes")
    for item in plan:
        print(f"  {'ADD' if item['before'] is None else 'UPDATE'} {item['path']}")
    if options.verify and plan:
        raise BundleError('tree has recognized older sources; apply the upgrade first')
    if options.verify or options.check:
        return 0
    if options.apply and plan:
        install(root, plan)
        if prepare(root, registry, options.with_fftw_compat):
            raise BundleError('post-install verification failed')
        if options.with_run_lapw and prepare_driver(root, registry['run_lapw']):
            raise BundleError('post-install driver verification failed')
        plan = []
        print('Complete current source hashes and both object lists verified.')
    if options.build or options.install:
        if plan:
            raise BundleError('use --apply before --build on older sources')
        build(root, options.with_fftw_compat, options.with_lapw0_mpi)
    if options.install:
        binaries = [('SRC_lapw0', 'lapw0'), ('SRC_lapw1', 'lapw1'), ('SRC_lapw1', 'lapw1c'),
                    ('SRC_lcore', 'lcore'), ('SRC_mixer', 'mixer'),
                    ('SRC_lapw2', 'lapw2'), ('SRC_lapw2', 'lapw2c')]
        if options.with_lapw0_mpi:
            binaries.append(('SRC_lapw0', 'lapw0_mpi'))
        if options.with_fftw_compat:
            binaries += [('SRC_3ddens', '3ddens'), ('SRC_hf', 'hf'), ('SRC_hf', 'hfc')]
        for directory, filename in binaries:
            shutil.copy2(root / directory / filename, root / filename)
        print('Installed rebuilt executables.')
    if options.install_response_workers:
        if plan:
            raise BundleError('use --apply before building workers on older sources')
        install_workers(root, options.install_response_workers.expanduser().resolve(), exports)
    if options.install_helpers:
        options.install_helpers.mkdir(parents=True, exist_ok=True)
        for name in ['rmt_stress_correction.py', 'radial_response_fd.py']:
            target = options.install_helpers / name
            shutil.copy2(BUNDLE / 'helpers' / name, target)
            target.chmod(0o755)
        print(f'Installed diagnostic helpers in {options.install_helpers}')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (BundleError, OSError, subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
