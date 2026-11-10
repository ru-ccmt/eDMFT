"""Check installed augmentation helpers and the rebuilt native radial archive."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import unittest


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--evidence', type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=False)
    source = args.root / 'src/python'
    tests = source / 'tests'
    sys.path[:0] = [str(source), str(tests)]
    import lapw_augmentation_response as augmentation
    assert Path(augmentation.__file__).resolve() == source / 'lapw_augmentation_response.py'
    suite = unittest.defaultTestLoader.discover(str(tests))
    with (args.evidence / 'installed_python.log').open('w') as log:
        result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit('Installed Python tests failed; see installed_python.log')
    print(f'{result.testsRun} installed-source Python tests passed', flush=True)

    # A new interpreter prevents previously imported source modules from
    # disguising an outdated or missing dbin installation.
    env = os.environ.copy()
    env['PYTHONPATH'] = '/Users/haule/dbin:' + str(tests)
    with (args.evidence / 'installed_dbin.log').open('w') as log:
        subprocess.run([sys.executable, '-m', 'unittest', '-v',
                        'test_lapw_augmentation_response'], cwd=args.work,
                       env=env, stdout=log, stderr=subprocess.STDOUT, check=True)

    native = args.root / 'src/dmft2/atpar'
    library = native / 'libatpar.a'
    executable = args.work / 'check_radial_potential'
    command = ['/opt/homebrew/bin/gfortran', '-O3', '-fopenmp',
               '-ffree-line-length-none', '-I' + str(native),
               str(args.root / 'src/dmft2/tests/test_radial_potential.f90'),
               str(library), '-o', str(executable)]
    with (args.evidence / 'installed_native_link.log').open('w') as log:
        subprocess.run(command, cwd=args.work, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    fixtures = args.root / 'src/doc/stress/native_hamiltonian_response_20261108/native_fixtures'
    staged = json.loads((args.evidence / 'evidence/native_checked_report.json').read_text())
    rows = []
    for pair in staged['rows']:
        run = subprocess.run([str(executable),
                              str(fixtures / (pair['left'] + '.input')),
                              str(fixtures / (pair['right'] + '.input'))],
                             capture_output=True, text=True, check=True)
        values = list(map(float, run.stdout.split()))
        if len(values) != 3:
            raise ValueError('unexpected native test output')
        rows.append(dict(left=pair['left'], right=pair['right'],
                         central=values[0], finite_difference=values[1],
                         adjoint=values[2]))
    report = dict(python_tests=result.testsRun, dbin_tests=9,
                  native_cases=len(rows), native_rows=rows,
                  maximum_errors={key: max(row[key] for row in rows)
                                  for key in ['central', 'finite_difference', 'adjoint']},
                  library_sha256=sha(library), compiler=command[0],
                  source_helper=str(augmentation.__file__),
                  source_helper_sha256=sha(Path(augmentation.__file__)),
                  dbin_helper_sha256=sha(Path('/Users/haule/dbin/lapw_augmentation_response.py')))
    assert report['source_helper_sha256'] == report['dbin_helper_sha256']
    (args.evidence / 'installed_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'native_rows'}, indent=2))


if __name__ == '__main__':
    main()
