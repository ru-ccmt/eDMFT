"""Validate the proposed source patch without writing to the user's tree."""
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import unittest

WORK = Path(__file__).resolve().parent
ROOT = Path('/Users/haule/projects/eDMFT_new')
STAGE = WORK / 'stage'

def main():
    sys.path[:0] = [str(STAGE), str(ROOT / 'src/python'),
                   str(ROOT / 'src/python/tests')]
    suite = unittest.defaultTestLoader.discover(str(ROOT / 'src/python/tests'))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromName(
        'test_lapw_augmentation_response'))
    with (WORK / 'python_regression.log').open('w') as log:
        result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit('Python tests failed; see python_regression.log')
    print(f'{result.testsRun} Python tests passed', flush=True)

    build = WORK / 'staged_build/src/dmft2/atpar'
    build.mkdir(parents=True, exist_ok=False)
    shutil.copy2(ROOT / 'src/Makefile.in', build.parents[1] / 'Makefile.in')
    for source in (ROOT / 'src/dmft2/atpar').glob('*.f90'):
        shutil.copy2(source, build)
    for name in ['Makefile', 'lapw_radial_potential.f90']:
        shutil.copy2(STAGE / name, build / name)
    with (WORK / 'staged_build.log').open('w') as log:
        subprocess.run(['make', '-j2', 'libatpar.a',
                        'PF90=/opt/homebrew/bin/mpif90'],
                       cwd=build, stdout=log, stderr=subprocess.STDOUT,
                       check=True)
    check = build / 'check_radial_potential'
    subprocess.run(['/opt/homebrew/bin/gfortran', '-O3', '-fopenmp',
                    '-ffree-line-length-none', '-I.',
                    str(STAGE / 'test_radial_potential.f90'),
                    str(build / 'libatpar.a'), '-o', str(check)],
                   cwd=build, check=True)
    fixtures = ROOT / 'src/doc/stress/native_hamiltonian_response_20261108/native_fixtures'
    prior = json.loads((WORK / 'native_pair_tests/report.json').read_text())
    rows = []
    for row in prior['rows']:
        run = subprocess.run([str(check),
                              str(fixtures / (row['left'] + '.input')),
                              str(fixtures / (row['right'] + '.input'))],
                             check=True, capture_output=True, text=True)
        errors = list(map(float, run.stdout.split()))
        if len(errors) != 3:
            raise ValueError('unexpected native output')
        rows.append(dict(left=row['left'], right=row['right'],
                         central=errors[0], finite_difference=errors[1],
                         adjoint=errors[2]))
    report = dict(python_tests=result.testsRun, native_cases=len(rows), rows=rows,
                  maximum_errors={name: max(r[name] for r in rows)
                                  for name in ('central', 'finite_difference', 'adjoint')})
    (WORK / 'staged_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'rows'}, indent=2))

if __name__ == '__main__':
    main()
