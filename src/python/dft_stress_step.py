#!/usr/bin/env python3
"""Compute final-only screened DFT stress from one converged WIEN iteration.

The last LAPW0/LAPW1/LAPW2 iteration used case.clmsum_old.  MIXER has since
replaced case.clmsum with the next density.  The isolated response snapshot
therefore restores clmsum_old, leaving the live calculation untouched.  It
reconstructs the native Mermin valence charge and differentiates the same
discrete charge and potential maps used for eDMFT, without an impurity
projector.  Only a converged driver invokes this program, once at the end.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
import uuid

import numpy as np

from dmft2_step import entries, local_definition
from lapw_dft_fermi import native_fermi_parameters
from response_mpi import resolve_launches

INPUT_SUFFIXES = (
    'struct', 'in0', 'in1', 'in1c', 'in2', 'in2c', 'inc', 'inm',
    'klist', 'kgen', 'rotlm', 'clmsum', 'vsp', 'vns',
)
THREADS = dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1', MKL_NUM_THREADS='1')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def in2_mode(path):
    return path.read_text().splitlines()[0].split()[0].upper()


def native_stress(text):
    rows = []
    for i in range(1, 4):
        found = re.findall(r'^:STRESS_GPa%03d:\s*(.*)$' % i, text, re.M)
        if not found:
            raise ValueError('native mixer stress row is missing')
        rows.append([float(x) for x in found[-1].split()[:3]])
    return np.array(rows)


def publication(report, reference_id, native):
    ry = np.asarray(report['stress_Ry_Bohr3'], float)
    gpa = np.asarray(report['stress_GPa'], float)
    if ry.shape != (3, 3) or gpa.shape != (3, 3) or not np.all(np.isfinite(gpa)):
        raise ValueError('response did not supply a finite 3x3 stress tensor')
    if not np.allclose(gpa, gpa.T, rtol=0, atol=1e-8):
        raise ValueError('response stress tensor is not symmetric')
    pressure = float(np.trace(gpa) / 3)
    lines = [f'\n:DFT_STRESS_RESPONSE: {reference_id} final screened DFT stress\n']
    lines.extend(':STRESS_RY%03d:' % i + ''.join('%20.10f' % x for x in row) + '\n'
                 for i, row in enumerate(ry, 1))
    lines.extend(':STRESS_GPa%03d:' % i + ''.join('%16.5f' % x for x in row)
                 + '   total\n' for i, row in enumerate(gpa, 1))
    lines.append(':PRESSURE: %18.5f GPa     total\n' % pressure)
    lines.append(':DFT_STRESS_RESPONSE_END: final-only publication\n')
    return ''.join(lines), dict(native_GPa=native.tolist(), screened_GPa=gpa.tolist(),
                               correction_GPa=(gpa-native).tolist(), pressure_GPa=pressure)


def run_step(command, cwd, log, *, env):
    with log.open('w') as stream:
        result = subprocess.run(command, cwd=cwd, env=env, stdout=stream,
                                stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'{shlex.join(command)} failed; see {log}')


def execute(case, *, root, wien_root, directory=Path('.')):
    directory = Path(directory).resolve()
    if not re.fullmatch(r'[A-Za-z0-9_.+-]+', case):
        raise ValueError('invalid WIEN case name')
    definition = directory/'lapw1.def'
    if not definition.is_file():
        raise ValueError('missing last-iteration lapw1.def')
    complex_case = entries(definition.read_text())[5].endswith('.in1c')
    in2 = directory / (case + ('.in2c' if complex_case else '.in2'))
    if in2_mode(in2) not in ('STRT', 'STRFT'):
        raise ValueError('final screened DFT stress requires STRT or STRFT in case.in2[c]')
    if in2_mode(directory / (case+'.in0')) != 'STR':
        raise ValueError('final screened DFT stress requires STR in case.in0')
    native_fermi_parameters(directory, case, complex_vectors=complex_case)  # Require Mermin TEMPS.
    old_charge = directory / (case+'.clmsum_old')
    if not old_charge.is_file():
        raise ValueError('missing pre-mixer case.clmsum_old for the final LAPW snapshot')
    scfm = directory / (case+'.scfm')
    scf = directory / (case+'.scf')
    if not scfm.is_file() or not scf.is_file():
        raise ValueError('MIXER must finish before the final response')
    if ':DFT_STRESS_RESPONSE:' in scfm.read_text():
        raise ValueError('this final DFT iteration already has a screened stress')
    native = native_stress(scfm.read_text())
    work = directory/(case+'.dft_response_runs')/(time.strftime('%Y%m%dT%H%M%S')+'-'+uuid.uuid4().hex[:12])
    snapshot = work/'reference'
    snapshot.mkdir(parents=True)
    fingerprints = {}
    for suffix in INPUT_SUFFIXES:
        source = old_charge if suffix == 'clmsum' else directory/(case+'.'+suffix)
        if source.is_file():
            target = snapshot/(case+'.'+suffix)
            shutil.copy2(source, target)
            fingerprints[source.name] = digest(source)
    (snapshot/'lapw1.def').write_text(local_definition(definition.read_text()))
    fingerprints['lapw1.def'] = digest(definition)
    reference_id = hashlib.sha256(json.dumps(fingerprints, sort_keys=True).encode()).hexdigest()
    prefix_file = directory/'mpi_prefix.dat'
    prefix = prefix_file.read_text().splitlines()[0].strip() if prefix_file.is_file() else ''
    mpi, native_mpi = resolve_launches(prefix, {}, analytic=True)
    workers = {name: root/name for name in ('lapw0_response', 'lcore_response',
                  'lapw1c_response' if complex_case else 'lapw1_response')}
    for path in [*workers.values(), root/'liblapw_response.so', wien_root/'x']:
        if not path.is_file():
            raise FileNotFoundError(path)
    config = dict(case=case, electronic_model='dft', reference=str(snapshot),
                  workspace=str(work/'central'), wien_root=str(wien_root),
                  lapw0_worker=str(workers['lapw0_response']),
                  lapw1_worker=str(workers['lapw1c_response' if complex_case else 'lapw1_response']),
                  lcore_worker=str(workers['lcore_response']),
                  lapw1_definition=str(snapshot/'lapw1.def'),
                  lapw1_mpi=mpi, native_mpi=native_mpi)
    (work/'prepare.json').write_text(json.dumps(config, indent=2)+'\n')
    env = os.environ | THREADS
    prepare = [sys.executable, str(root/'lapw_response_prepare.py'), str(work/'prepare.json')]
    run_step(prepare, directory, work/'prepare.log', env=env)
    command = [sys.executable, str(root/'lapw_total_stress_driver.py'),
               '--snapshot', str(work/'central'), '--case', case,
               '--potential-core', str(work/'central/potential_core.npz'),
               '--library', str(root/'liblapw_response.so'),
               '--output', str(work/'total_stress')]
    if mpi:
        command = mpi + command + ['--mpi']
    run_step(command, directory, work/'response.log', env=env)
    path = work/'total_stress/report.json'
    report = json.loads(path.read_text())
    if (report.get('converged') is not True or
        report.get('full_tensor_within_symmetry') is not True or
        report.get('diagnostic_only', False)):
        raise RuntimeError('screened DFT stress is incomplete; see '+str(path))
    for name, expected in fingerprints.items():
        if digest(directory/name) != expected:
            raise RuntimeError('live WIEN input changed during response: '+name)
    correction, comparison = publication(report, reference_id, native)
    summary = dict(reference_id=reference_id, report=str(path),
                   workspace=str(work), **comparison)
    (work/'comparison.json').write_text(json.dumps(summary, indent=2)+'\n')
    with scfm.open('a') as stream:
        stream.write(correction)
    with scf.open('a') as stream:
        stream.write(correction)
    reference = report.get('reference_fermi_response') or {}
    adjoint = report.get('adjoint') or {}
    (directory/(case+'.output2_response')).write_text(
        'Final screened DFT stress\n'
        f'Work area: {work}\nAnalytic report: {path}\n'
        f'Native mean tensile stress: {np.trace(native)/3:.8f} GPa\n'
        f'Screened mean tensile stress: {comparison["pressure_GPa"]:.8f} GPa\n'
        'Tensor correction in GPa:\n'
        + np.array2string(np.asarray(comparison['correction_GPa']), precision=8) + '\n'
        f'Hydrostatic explicit derivative: {report["explicit_Ry"][0]:.12g} Ry\n'
        f'Hydrostatic induced derivative: {report["induced_Ry"][0]:.12g} Ry\n'
        f'IN1 reference EF: {reference.get("in1_reference_Ry", float("nan")):.12g} Ry\n'
        f'EF-relative radial channels: {len(reference.get("radial_slots", []))}\n'
        f'Adjoint residual: {adjoint.get("residual_norm", float("nan")):.6g}\n'
        f'MPI ranks: {report.get("parallel", {}).get("ranks", 1)}\n'
        f'Response time: {report.get("timings", {}).get("total_seconds", float("nan")):.1f} s\n'
        'Preparation and response logs are in the work area.\n'
    )
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', required=True)
    p.add_argument('--root', type=Path, default=Path(os.environ.get('WIEN_DMFT_ROOT', Path(__file__).parent)))
    p.add_argument('--wien-root', type=Path, default=Path(os.environ.get('WIENROOT', '')))
    p.add_argument('--directory', type=Path, default=Path('.'))
    a = p.parse_args()
    result = execute(a.case, root=a.root.resolve(), wien_root=a.wien_root.resolve(), directory=a.directory)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
