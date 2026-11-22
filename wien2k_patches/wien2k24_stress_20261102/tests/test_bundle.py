#!/usr/bin/env python3
"""Exercise deployment in disposable source trees; never modify --pristine-root.

Usage: python3 tests/test_bundle.py --pristine-root /path/to/WIEN2k24 -v
Requires the original source files of the seven covered components (including
optional 3ddens/HF), plus their Makefiles. No compiler or WIEN binaries needed.
"""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.dont_write_bytecode = True
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--pristine-root', required=True, type=Path)
parser.add_argument('--bundle', type=Path, default=Path(__file__).resolve().parents[1])
OPTIONS, EXTRA = parser.parse_known_args()
B = OPTIONS.bundle.resolve()
REG = json.loads((B / 'sources.json').read_text())
spec = importlib.util.spec_from_file_location('wien_bundle_installer', B / 'apply_wien2k24_stress.py')
installer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installer)


def snapshot(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob('*') if p.is_file()}


def patch(root, path, reverse=False):
    command = ['patch', '--batch', '-F', '2', '-p1', '-d', str(root), '-i', str(path)]
    if reverse:
        command.append('-R')
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        raise AssertionError(result.stdout + result.stderr)


class BundleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.storage = tempfile.TemporaryDirectory(prefix='wien-bundle-test-')
        cls.top = Path(cls.storage.name)
        cls.base = cls.top / 'original'
        for relative in [*REG['sources'], *[r['path'] for r in REG['makefiles']],
                         *REG['response_workers']['sources'], 'SRC/run_lapw', 'run_lapw',
                         'SRC_3ddens/Makefile', 'SRC_hf/Makefile']:
            source = OPTIONS.pristine_root / relative
            if source.is_file():
                dest = cls.base / relative
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, dest)
        if not (cls.base / 'SRC/run_lapw').is_file():
            (cls.base / 'SRC').mkdir(exist_ok=True)
            shutil.copy2(cls.base / 'run_lapw', cls.base / 'SRC/run_lapw')
        (cls.base / 'run_lapw').unlink(missing_ok=True)
        for line in (B / 'BASE_SOURCES.sha256').read_text().splitlines():
            digest, relative = line.split(None, 1)
            assert installer.sha((cls.base / relative).read_bytes()) == digest, relative
        cls.pristine_fingerprint = snapshot(OPTIONS.pristine_root)

    @classmethod
    def tearDownClass(cls):
        assert snapshot(OPTIONS.pristine_root) == cls.pristine_fingerprint
        cls.storage.cleanup()

    def setUp(self):
        self.storage = tempfile.TemporaryDirectory(dir=self.top)
        self.root = Path(self.storage.name) / 'case with spaces'
        shutil.copytree(self.base, self.root)

    def tearDown(self):
        self.storage.cleanup()

    def run_tool(self, *args, good=True, bundle=B):
        result = subprocess.run([sys.executable, str(bundle / 'apply_wien2k24_stress.py'),
                                 *args, str(self.root)], capture_output=True, text=True)
        self.assertEqual(result.returncode == 0, good, result.stdout + result.stderr)
        return result

    def assert_current(self, optional=False):
        for relative, item in REG['sources'].items():
            if optional or not item['optional']:
                self.assertEqual(installer.sha((self.root / relative).read_bytes()), item['target_sha256'])
        self.run_tool('--verify', *(['--with-fftw-compat'] if optional else []))

    def test_01_pristine_check_is_read_only(self):
        before = snapshot(self.root)
        self.run_tool('--check')
        self.assertEqual(snapshot(self.root), before)

    def test_02_pristine_apply_backup_and_idempotence(self):
        before = {p: (self.root / p).read_bytes() if (self.root / p).exists() else None
                  for p in REG['sources']}
        self.run_tool('--apply')
        self.assert_current()
        backups = list(self.root.glob('stress_patch_backup_*'))
        self.assertEqual(len(backups), 1)
        for p, data in before.items():
            if data is not None and not REG['sources'][p]['optional']:
                self.assertEqual((backups[0] / p).read_bytes(), data)
        first = snapshot(self.root)
        self.run_tool('--apply')
        self.assertEqual(snapshot(self.root), first)

    def historical_tree(self, state):
        shutil.rmtree(self.root)
        shutil.copytree(self.base, self.root)
        if state == 'pristine':
            return
        if state.startswith('legacy'):
            for p in sorted((B / 'patches/required').glob('*.patch')):
                patch(self.root, p)
            if state == 'legacy-full':
                for p in sorted((B / 'patches/optional').glob('010*.patch')):
                    patch(self.root, p)
        else:
            h = B / 'patches/history/20261105'
            patch(self.root, h / '0001-wien2k24-complete-stress.patch')
            if state == 'complete-pre-cache':
                patch(self.root, h / '0004-upgrade-complete-to-operator-cache.patch', True)

    def test_03_all_historical_whole_tree_upgrades(self):
        for state in ['legacy-required', 'legacy-full', 'complete-pre-cache', 'complete-20261105']:
            with self.subTest(state=state):
                self.historical_tree(state)
                self.run_tool('--apply')
                self.assert_current()

    def test_04_every_registered_variant(self):
        # Reverse each transition from its target to recover its exact baseline,
        # then test the public installer on that documented mixed-state tree.
        self.run_tool('--apply', '--with-fftw-compat')
        n = 0
        for relative, item in REG['sources'].items():
            for variant in item['variants']:
                with self.subTest(path=relative, sha=variant['sha256']):
                    patch(self.root, B / variant['patch'], True)
                    actual = installer.sha((self.root / relative).read_bytes()) if (self.root / relative).exists() else None
                    self.assertEqual(actual, variant['sha256'])
                    self.run_tool('--apply', '--with-fftw-compat')
                    self.assert_current(True)
                    n += 1
        self.assertGreaterEqual(n, 55)

    def test_05_site_flags_and_comments_preserved(self):
        original = {}
        for rule in REG['makefiles']:
            p = self.root / rule['path']
            data = p.read_bytes()
            data = b'# site example: ' + rule['object'].encode() + b' in a comment\nFC = mpif90-site\nLIBS = -L/site/lib -lsite\n' + data
            p.write_bytes(data)
            original[rule['path']] = data
        self.run_tool('--apply')
        for rule in REG['makefiles']:
            self.assertEqual((self.root / rule['path']).read_bytes(), installer.add_object(original[rule['path']], rule))

    def test_06_unknown_source_refused_without_partial_writes(self):
        p = self.root / 'SRC_mixer/scfana.f'
        p.write_bytes(p.read_bytes() + b'\n! independent local change\n')
        before = snapshot(self.root)
        self.run_tool('--apply', good=False)
        self.assertEqual(before, snapshot(self.root))

    def test_07_ambiguous_makefile_refused_without_partial_writes(self):
        p = self.root / 'SRC_lapw2/Makefile'
        p.write_bytes(p.read_bytes() + b'\nOBJS += local.o\n')
        before = snapshot(self.root)
        self.run_tool('--apply', good=False)
        self.assertEqual(before, snapshot(self.root))

    def test_08_missing_component_refused_without_partial_writes(self):
        (self.root / 'SRC_lcore/hfsd.f').unlink()
        before = snapshot(self.root)
        self.run_tool('--apply', good=False)
        self.assertEqual(before, snapshot(self.root))

    def test_09_optional_fftw_group(self):
        before = {p: (self.root / p).read_bytes() for p, x in REG['sources'].items() if x['optional']}
        self.run_tool('--apply')
        for p, data in before.items():
            self.assertEqual((self.root / p).read_bytes(), data)
        self.run_tool('--apply', '--with-fftw-compat')
        self.assert_current(True)

    def test_10_integrity_detects_corruption(self):
        corrupted = Path(self.storage.name) / 'damaged bundle'
        shutil.copytree(B, corrupted)
        p = corrupted / 'helpers/radial_response_fd.py'
        p.write_bytes(p.read_bytes() + b'\n# damage\n')
        before = snapshot(self.root)
        self.run_tool('--apply', good=False, bundle=corrupted)
        self.assertEqual(before, snapshot(self.root))

    def test_11_all_reviewable_cumulative_patches(self):
        for state in ('pristine', 'legacy-required', 'legacy-full', 'complete-pre-cache', 'complete-20261105'):
            with self.subTest(state=state):
                self.historical_tree(state)
                (upgrade,) = (B / 'patches/current').glob('*-' + state + '-to-current.patch')
                patch(self.root, upgrade)
                for name, record in REG['sources'].items():
                    if not record['optional']:
                        self.assertEqual(installer.sha((self.root / name).read_bytes()),
                                         record['target_sha256'], name)

    def test_12_failed_install_rolls_back(self):
        plan = installer.prepare(self.root, REG)
        originals = {p['path']: p['before'] for p in plan}
        real_replace = installer.os.replace
        calls = 0
        def failing_replace(*args):
            nonlocal calls
            calls += 1
            if calls == 5:
                raise OSError('injected replacement failure')
            return real_replace(*args)
        with mock.patch.object(installer.os, 'replace', side_effect=failing_replace):
            with self.assertRaises(OSError):
                installer.install(self.root, plan)
        for relative, data in originals.items():
            path = self.root / relative
            self.assertEqual(path.read_bytes() if path.exists() else None, data)

    def test_13_current_sources_missing_object_is_repaired(self):
        self.run_tool('--apply')
        rule = REG['makefiles'][0]
        p = self.root / rule['path']
        p.write_bytes(p.read_bytes().replace(b' ' + rule['object'].encode(), b'', 1))
        self.run_tool('--verify', good=False)
        self.run_tool('--apply')
        self.assert_current()

    def test_14_legacy_entry_points_are_current(self):
        r = subprocess.run(['bash', str(B / 'apply_wien2k24_stress_patches_legacy.sh'),
                            '--apply', str(self.root)], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assert_current()

    def test_15_final_driver_and_worker_preflight(self):
        destination = self.root.parent / 'workers'
        flags = ('--with-run-lapw', '--install-response-workers', str(destination))
        before = snapshot(self.root)
        self.run_tool('--check', '--install', '--install-helpers', str(destination), *flags)
        self.assertEqual(snapshot(self.root), before)
        self.assertFalse(destination.exists())
        script = self.root / 'SRC/run_lapw'
        mode = script.stat().st_mode
        (self.root / 'run_lapw').symlink_to('SRC/run_lapw')
        self.run_tool('--apply', '--with-run-lapw')
        self.assertTrue((self.root / 'run_lapw').is_symlink())
        self.assertEqual(script.stat().st_mode, mode)
        self.assertEqual(installer.sha(script.read_bytes()), REG['run_lapw'][0]['target_sha256'])
        after = snapshot(self.root)
        self.run_tool('--verify', *flags)
        self.run_tool('--apply', '--with-run-lapw')
        self.assertEqual(snapshot(self.root), after)
        self.assertFalse(destination.exists())

    def test_16_unknown_driver_or_worker_refused_before_apply(self):
        for name in ('run_lapw', 'SRC_lcore/resld.f'):
            with self.subTest(name=name):
                self.historical_tree('pristine')
                file = self.root / name
                file.write_bytes((file.read_bytes() if file.exists() else b'') + b'# site change')
                before = snapshot(self.root)
                self.run_tool('--apply', '--with-run-lapw', '--install-response-workers',
                              str(self.root.parent / 'workers'), good=False)
                self.assertEqual(snapshot(self.root), before)

    def test_17_failed_worker_build_preserves_binaries(self):
        self.run_tool('--apply')
        destination = self.root.parent / 'workers'
        destination.mkdir()
        (destination / 'lcore_response').write_bytes(b'previous build')
        (self.root / 'SRC_Globals').mkdir()
        before = snapshot(self.root)
        recipe = REG['response_workers']
        exports = installer.prepare_worker_sources(self.root, recipe, [])
        # Fail at the third target, after two new products have been built.
        def make(args, check):
            component = Path(args[args.index('-C') + 1])
            name = {'seq': 'lapw0', 'lcore': 'lcore'}.get(args[-1])
            if name is None:
                raise subprocess.CalledProcessError(1, args)
            path = component / name
            path.write_bytes(b'test executable')
            path.chmod(0o755)
        with mock.patch.object(installer.subprocess, 'run', side_effect=make):
            with self.assertRaises(subprocess.CalledProcessError):
                installer.install_workers(self.root, destination, exports)
        self.assertEqual(list(destination.iterdir()), [destination / 'lcore_response'])
        self.assertEqual((destination / 'lcore_response').read_bytes(), b'previous build')
        self.assertEqual(snapshot(self.root), before)

    def test_18_mpi_lapw0_is_opt_in(self):
        for enabled in (False, True):
            with mock.patch.object(installer.subprocess, 'run') as run:
                installer.build(self.root, False, enabled)
                targets = [call.args[0][-1] for call in run.call_args_list]
                self.assertEqual(targets.count('para'), int(enabled))


if __name__ == '__main__':
    unittest.main(argv=[sys.argv[0], *EXTRA])
