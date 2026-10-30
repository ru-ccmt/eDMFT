#!/usr/bin/env python
"""Continue CTQMC/BCTQMC self-energies and optionally prepare A(k,w) and FS data."""
import argparse
import ast
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import utils


def ChangeMatsubaraFlag(fname, matsubara, x_range=None, nfrq=None):
    with open(fname) as fi:
        lines = fi.readlines()
    fields, sep, comment = lines[1].partition('#')
    dat = fields.split()
    dat[0] = str(matsubara)
    if x_range is not None:
        dat[4:6] = map(str, x_range)
    if nfrq is not None:
        dat[3] = str(nfrq)
    lines[1] = ' '.join(dat) + ('  #' + comment if sep else '\n')
    with open(fname, 'w') as fo:
        fo.writelines(lines)


def read_solver(fname):
    """Read the conventional solver='...' assignment without executing params.dat."""
    solver = 'CTQMC'  # Same default as run_dmft.py.
    for node in ast.parse(fname.read_text(), filename=str(fname)).body:
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AnnAssign):
            targets = [node.target]
        else:
            continue
        if any(isinstance(t, ast.Name) and t.id == 'solver' for t in targets):
            try:
                solver = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                raise ValueError("params.dat: solver must be a literal 'CTQMC' or 'BCTQMC'")
    if solver not in ('CTQMC', 'BCTQMC'):
        raise ValueError('Unsupported solver in params.dat: {!r}'.format(solver))
    return solver


def read_bctqmc_estimator(fname, impurity=0):
    """Read only estimator controls; do not execute arbitrary params.dat code."""
    name = 'iparams' + str(impurity)
    controls = {}
    for node in ast.parse(fname.read_text(), filename=str(fname)).body:
        targets = node.targets if isinstance(node, ast.Assign) else (
            [node.target] if isinstance(node, ast.AnnAssign) else [])
        if any(isinstance(t, ast.Name) and t.id == name for t in targets):
            if not isinstance(node.value, ast.Dict):
                raise ValueError('params.dat: {} must be a literal dictionary'.format(name))
            controls = {key.value: value for key, value in zip(node.value.keys, node.value.values)
                        if isinstance(key, ast.Constant)}
        for target in targets:
            if (isinstance(target, ast.Subscript) and isinstance(target.value, ast.Name)
                    and target.value.id == name and isinstance(target.slice, ast.Constant)):
                controls[target.slice.value] = node.value

    def value(key, default):
        node = controls.get(key)
        if node is None:
            return default
        if isinstance(node, (ast.List, ast.Tuple)):
            if not node.elts:
                raise ValueError('params.dat: {} must not be empty'.format(key))
            node = node.elts[0]
        try:
            return ast.literal_eval(node)
        except (ValueError, TypeError):
            raise ValueError('params.dat: {} must be a literal value'.format(key))

    estimator = value('self_energy_estimator', 'legacy')
    if not isinstance(estimator, str):
        raise ValueError('params.dat: self_energy_estimator must be a string')
    estimator = estimator.strip().lower()
    if estimator == 'legacy':
        blend = value('Sdyson_Bulla_blend', True)
        if not isinstance(blend, bool):
            raise ValueError('params.dat: Sdyson_Bulla_blend must be a boolean')
        return 'blended' if blend else 'bulla'
    if estimator == 'combined':
        estimator = 'blended'
    if estimator not in ('dyson', 'blended', 'bulla', 'kugler'):
        raise ValueError('Unknown self_energy_estimator: {!r}'.format(estimator))
    return estimator


def prepare_bulla_impurity(source, processed, trans, sigind, target):
    """Put selected Bulla components in the processed impurity sig.out convention."""
    with processed.open() as fi:
        headers = [fi.readline(), fi.readline()]
    edc = np.asarray(ast.literal_eval(headers[1].split('=', 1)[1]), dtype=float)
    previous_sinf = np.asarray(ast.literal_eval(headers[0].split('=', 1)[1]), dtype=float)
    raw = np.loadtxt(source, ndmin=2)
    old = np.loadtxt(processed, ndmin=2)
    if (raw.shape[1] != 1 + 2*len(edc) or old.shape[1] != raw.shape[1]
            or len(previous_sinf) != len(edc) or not np.all(np.isfinite(raw))
            or not np.all(np.isfinite(edc)) or np.any(np.diff(raw[:, 0]) <= 0)
            or raw[0, 0] <= 0):
        raise ValueError('{}: inconsistent impurity self-energy layout'.format(source))
    with trans.open() as fi:
        dimension = int(fi.readline().split()[0])
        fi.readline()  # Sigind heading.
        local = np.array([list(map(int, fi.readline().split()[:dimension]))
                          for _ in range(dimension)], dtype=int)
    if local.shape != sigind.shape:
        raise ValueError('{}: Trans.dat and case.indmfi dimensions differ'.format(source))
    labels = sorted(set(sigind.ravel()) - {0})
    labels = [label for label in labels if label > 0]
    columns = []
    diagonal = []
    for label in labels:
        candidates = set(local[sigind == label])
        if len(candidates) != 1:
            raise ValueError('{}: ambiguous impurity mapping for label {}'.format(source, label))
        column = candidates.pop() - 1
        if column < 0 or column >= len(edc):
            raise ValueError('{}: missing Bulla component for label {}'.format(source, label))
        columns.append(column)
        diagonal.append(label in np.diag(sigind))
    data = raw[:, [0] + [k for c in columns for k in (1+2*c, 2+2*c)]].copy()
    sinf = data[-1, 1::2].copy()
    dc = edc[columns].copy()
    for i, column in enumerate(columns):
        # Honor components suppressed by the existing impurity processing.
        suppressed = (not diagonal[i] and previous_sinf[column] == 0 and dc[i] == 0
                      and np.all(old[:, 1+2*column:3+2*column] == 0))
        projected = diagonal[i] and dc[i] <= -1000
        if suppressed or projected:
            sinf[i] = dc[i]
            data[:, 1+2*i] = dc[i]
            data[:, 2+2*i] = 0
    data[:, 1::2] -= dc  # ctqmc.HighFrequency writes Sigma - Edc.
    with target.open('w') as fo:
        fo.write('# s_oo= {!r}\n# Edc= {!r}\n'.format(sinf.tolist(), dc.tolist()))
        np.savetxt(fo, data, fmt='%.15e')


def prepare_continuation_sigma(job, directory, case, root, siginds, bulla_impurities):
    """Use assembled sig.inp normally; gather exceptional Bulla inputs in isolation."""
    target = directory / 'sig.inpx'
    if not bulla_impurities:
        shutil.copy2(job / 'sig.inp', target)
        print('Continuing assembled self-energy:', job / 'sig.inp', flush=True)
        return target.name
    with tempfile.TemporaryDirectory(prefix='.bulla-', dir=directory) as temporary:
        stage = Path(temporary)
        for suffix in ('.struct', '.indmfl', '.indmfldn', '.indmfi'):
            source = job / (case + suffix)
            if source.exists():
                (stage / source.name).symlink_to(source)
        # Preserve static-only components not produced by an impurity problem.
        shutil.copy2(job / 'sig.inp', stage / 'sig.inp')
        for impurity, sigind in siginds.items():
            original = job / ('imp.' + str(impurity))
            destination = stage / original.name
            destination.mkdir()
            if impurity in bulla_impurities:
                prepare_bulla_impurity(original / 'Sbulla.out', original / 'sig.out',
                                       original / 'Trans.dat', sigind, destination / 'sig.out')
            else:
                shutil.copy2(original / 'sig.out', destination / 'sig.out')
        # Reuse the real global component mapping and static subtraction.
        run([root / 'sgather.py', '-m', '1.0'], stage)
        shutil.copy2(stage / 'sig.inp', target)
    print('Continuing Bulla for impurities:', sorted(bulla_impurities), flush=True)
    return target.name


def run(command, directory, input_text=None):
    print('executing in {}: {}'.format(directory, shlex.join(map(str, command))), flush=True)
    subprocess.run(list(map(str, command)), cwd=str(directory), input=input_text,
                   text=True, check=True)


def prepare_maxent_params(job, directory, solver, x0):
    target = directory / 'maxent_params.dat'
    if not target.exists():
        source = job / 'maxent_params.dat'
        if source.exists():
            shutil.copy2(source, target)
        else:
            params = {
                'statistics': 'fermi', 'Ntau': 400, 'L': 30.0, 'x0': 0.01,
                'bwdth': 0.004, 'Nw': 450, 'gwidth': 2*15.0, 'idg': 1,
                'deltag': 0.01, 'Asteps': 4000, 'alpha0': 1000,
                'min_ratio': 0.001, 'iflat': 1,
                'Nitt': 500 if solver == 'BCTQMC' else 100, 'Nr': 0, 'Nf': 40,
            }
            if solver == 'BCTQMC':
                params['SymCum'] = True
            target.write_text('params = ' + repr(params) + '\n')
    if x0 is not None:
        # Append an override to preserve arbitrary Python expressions and comments.
        with target.open('a') as fo:
            fo.write("\nparams['x0'] = {!r}  # to_real_axis.py --x0\n".format(x0))


def read_real_sigma(fname):
    with fname.open() as fi:
        headers = [fi.readline(), fi.readline()]
    sizes = []
    for line, key in zip(headers, ('s_oo', 'Edc')):
        match = re.fullmatch(r'\s*#\s*' + key + r'\s*=\s*(.*?)\s*', line)
        if not match:
            raise ValueError('{}: expected two headers # s_oo=... and # Edc=...'.format(fname))
        values = np.asarray(ast.literal_eval(match.group(1)), dtype=float)
        if values.ndim != 1 or not np.all(np.isfinite(values)):
            raise ValueError('{}: invalid {} header'.format(fname, key))
        sizes.append(len(values))
    data = np.loadtxt(fname, ndmin=2)
    if not data.size or data.shape[1] < 3 or data.shape[1] % 2 != 1:
        raise ValueError('{}: expected frequency followed by Re/Im self-energy pairs'.format(fname))
    if any(n != (data.shape[1]-1)//2 for n in sizes):
        raise ValueError('{}: header lengths do not match self-energy columns'.format(fname))
    if not np.all(np.isfinite(data)) or np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError('{}: data must be finite with strictly increasing frequencies'.format(fname))
    if data[0, 0] >= 0:
        raise ValueError('{}: first data frequency must be negative (real-axis self-energy required)'.format(fname))
    return headers, data


def zero_frequency_sigma(fname):
    headers, data = read_real_sigma(fname)
    omega = data[:, 0]
    exact = np.flatnonzero(omega == 0.0)
    if len(exact):
        row = data[exact[0]].copy()
    else:
        below = np.flatnonzero(omega < 0)[-2:]
        above = np.flatnonzero(omega > 0)[:2]
        if len(below) != 2 or len(above) != 2:
            raise ValueError('{}: cubic interpolation needs two frequencies on each side of zero'.format(fname))
        points = data[np.concatenate((below, above))]
        # Lagrange evaluation of the cubic at zero, for every Re/Im column.
        weights = np.ones(4)
        for i in range(4):
            for j in range(4):
                if i != j:
                    weights[i] *= -points[j, 0] / (points[i, 0] - points[j, 0])
        row = np.concatenate(([0.0], weights @ points[:, 1:]))
    if not np.all(np.isfinite(row)):
        raise ValueError('{}: nonfinite zero-frequency interpolation'.format(fname))
    with fname.open('w') as fo:
        fo.writelines(headers)
        np.savetxt(fo, row[None, :], fmt='%.15e')


def prepare_mpi(source, target, ncores):
    if source.exists():
        text = source.read_text().strip()
        count = r'(?<!\S)(-n|-np|--np|--n|--ntasks)(?:\s+|=)(?:\d+|\$\{?NSLOTS\}?|[\"\']\$\{?NSLOTS\}?[\"\'])'
        text, replaced = re.subn(count, lambda m: m.group(1) + ' ' + str(ncores), text)
        if not replaced:
            text += ' -n {}'.format(ncores)
    else:
        text = 'mpirun -n {}'.format(ncores)
    target.write_text(text + '\n')


def cleanup_output(directory, case):
    """Remove generated executables, large vectors, and empty files (not directories)."""
    intermediates = {'dmft', 'lapw1', 'lapwso'}
    intermediates.update(case + suffix for suffix in ('.vector', '.vectorso', '.vectorsodn'))
    removed = []
    for path in sorted(directory.iterdir()):
        if path.name in intermediates and (path.is_file() or path.is_symlink()):
            path.unlink()  # For an executable symlink, remove only the local link.
            removed.append(path.name)
        elif path.is_file() and path.stat().st_size == 0:
            path.unlink()
            removed.append(path.name)
    if removed:
        print('Cleaned {}: {}'.format(directory, ', '.join(removed)), flush=True)


def positive_int(value):
    value = int(value)
    if value <= 0:
        raise argparse.ArgumentTypeError('must be positive')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fname', nargs='?', default='gc', help='legacy positional argument (unused)')
    parser.add_argument('-y', default=None, help='spectral range in eV, e.g. -y=-1.125:1.0')
    parser.add_argument('-a', default=3, type=positive_int, help='legacy option (unused; continue the current case-level sig.inp)')
    parser.add_argument('-o', default='onreal', help='real-axis directory (default: onreal)')
    parser.add_argument('-m', default='maxent', help='continuation directory (default: maxent)')
    parser.add_argument('-n', default=None, type=positive_int, help='spectral frequency point count')
    parser.add_argument('--x0', type=float, help='override maxent low-energy cutoff (new-file default: 0.01)')
    parser.add_argument('--ncores', default=10, type=positive_int, help='MPI process count (default: 10; in a job pass "$NSLOTS")')
    parser.add_argument('--Akw', action='store_true', help='prepare A(k,w); requires case.klist_band before execution')
    parser.add_argument('--FS', action='store_true', help='prepare FS/case.bxsf for local XCrySDen viewing')
    parser.add_argument('--FS-nk', default=5000, type=positive_int, help='requested FS k-point count, unshifted (default: 5000)')
    args = parser.parse_args()
    if args.x0 is not None and (not np.isfinite(args.x0) or args.x0 <= 0):
        parser.error('--x0 must be finite and positive')
    xrng = None
    if args.y is not None:
        try:
            lo, hi = args.y.split(':')
            xrng = [float(lo) if lo else -3.0, float(hi) if hi else 1.0]
            if not np.all(np.isfinite(xrng)) or xrng[0] >= xrng[1]:
                raise ValueError()
        except ValueError:
            parser.error('-y must specify increasing finite bounds, e.g. -y=-1.125:1.0')

    job = Path.cwd()
    maxent = (job / args.m).resolve()
    onreal = (job / args.o).resolve()
    fsdir = job / 'FS'
    directories = [job, maxent, onreal] + ([fsdir] if args.FS else [])
    if len(set(directories)) != len(directories) or any(
            a in b.parents or b in a.parents for i, a in enumerate(directories[1:])
            for b in directories[1:][i+1:]):
        raise ValueError('Job, continuation, real-axis, and FS directories must be distinct; output directories must not contain each other')
    solver = read_solver(job / 'params.dat')
    w2k = utils.W2kEnvironment()
    dmfe = utils.DmftEnvironment()
    root = Path(dmfe.ROOT)
    case = w2k.case
    band = job / (case + '.klist_band')
    if args.Akw and (not band.is_file() or not band.stat().st_size):
        raise ValueError('--Akw requires {}. Please supply it in the job directory before starting.'.format(band))
    for name in (case + '.struct', case + '.indmfl'):
        if not (job / name).is_file():
            raise ValueError('Missing input: ' + str(job / name))
    if not (job / 'sig.inp').is_file():
        raise ValueError('Missing assembled self-energy: ' + str(job / 'sig.inp'))
    siginds, bulla_impurities = {}, set()
    if solver == 'BCTQMC':
        from indmffile import ParsIndmfi
        siginds = ParsIndmfi(case)
        bulla_impurities = {i for i in siginds
                            if read_bctqmc_estimator(job / 'params.dat', i) in ('dyson', 'blended')}
        if bulla_impurities:
            for impurity in siginds:
                names = ('sig.out', 'Sbulla.out', 'Trans.dat') if impurity in bulla_impurities else ('sig.out',)
                for name in names:
                    source = job / ('imp.' + str(impurity)) / name
                    if not source.is_file():
                        raise ValueError('Missing BCTQMC input: ' + str(source))
    executables = [root / name for name in ('maxent_run.py', 'dmft_copy.py', 'x_dmft.py')]
    executables.append(Path(w2k.WIENROOT) / 'x_lapw')
    if bulla_impurities:
        executables.append(root / 'sgather.py')
    if args.FS:
        executables.append(root / 'FSconvert.py')
    for executable in executables:
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise ValueError('Missing or non-executable program: ' + str(executable))

    print('Postprocessing solver:', solver, flush=True)
    maxent.mkdir(parents=True, exist_ok=True)
    sigma_input = prepare_continuation_sigma(job, maxent, case, root, siginds, bulla_impurities)
    prepare_maxent_params(job, maxent, solver, args.x0)
    run([root / 'maxent_run.py', sigma_input], maxent)
    read_real_sigma(maxent / 'Sig.out')

    onreal.mkdir(parents=True, exist_ok=True)
    run([root / 'dmft_copy.py', os.path.relpath(job, onreal)], onreal)
    shutil.copy2(maxent / 'Sig.out', onreal / 'sig.inp')
    read_real_sigma(onreal / 'sig.inp')
    prepare_mpi(job / 'mpi_prefix.dat', onreal / 'mpi_prefix.dat', args.ncores)
    magnetic = (onreal / (case + '.indmfldn')).is_file()
    for suffix in ('.indmfl', '.indmfldn') if magnetic else ('.indmfl',):
        ChangeMatsubaraFlag(onreal / (case + suffix), 0, xrng, args.n)
    spin_orbit = (onreal / (case + '.inso')).is_file()
    xlapw = Path(w2k.WIENROOT) / 'x_lapw'
    xdmft = root / 'x_dmft.py'
    run([xlapw, 'lapw0', '-f', case], onreal)
    run([xdmft, 'lapw1'], onreal)
    if spin_orbit:
        run([xdmft, 'lapwso'], onreal)
    run([xdmft, 'dmft1'], onreal)
    if magnetic:
        run([xdmft, 'dmft1', '-l', 'dn'], onreal)

    if args.Akw:
        shutil.copy2(band, onreal / band.name)
        run([xdmft, 'lapw1', '--band'], onreal)
        if spin_orbit:
            run([xdmft, 'lapwso', '--band'], onreal)
        run([xdmft, 'dmftp'], onreal)
        if magnetic:
            run([xdmft, 'dmftp', '-l', 'dn'], onreal)

    if args.FS:
        if not (onreal / 'SOlapm.dat').is_file():
            raise ValueError('FS requires onreal/SOlapm.dat produced by dmft1')
        fsdir.mkdir(exist_ok=True)
        run([root / 'dmft_copy.py', os.path.relpath(onreal, fsdir)], fsdir)
        if not (fsdir / 'SOlapm.dat').is_file():
            raise ValueError('dmft_copy.py did not copy SOlapm.dat into FS/')
        run([xlapw, 'kgen', '-f', case], fsdir, '{}\n0\n'.format(args.FS_nk))
        shutil.copy2(onreal / 'mpi_prefix.dat', fsdir / 'mpi_prefix.dat')
        run([xlapw, 'lapw0', '-f', case], fsdir)
        run([xdmft, 'lapw1'], fsdir)
        if spin_orbit:
            run([xdmft, 'lapwso'], fsdir)
        zero_frequency_sigma(fsdir / 'sig.inp')
        run([xdmft, 'dmft0'], fsdir)
        run([root / 'FSconvert.py'], fsdir)
        bxsf = fsdir / (case + '.bxsf')
        if not bxsf.is_file():
            raise ValueError('FSconvert.py did not produce ' + str(bxsf))
        print('Fermi surface file: {}\nCopy it locally and open with: xcrysden --bxsf {}'.format(bxsf, bxsf.name))

    # Keep intermediates until all requested stages have completed successfully.
    cleanup_output(onreal, case)
    if args.FS:
        cleanup_output(fsdir, case)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, SyntaxError, subprocess.CalledProcessError) as error:
        sys.exit('to_real_axis.py: ' + str(error))
