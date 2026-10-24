#!/usr/bin/env python
"""Continue CTQMC/BCTQMC self-energies and optionally prepare A(k,w) and FS data."""
import argparse
import ast
import glob
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

import numpy as np
import utils
from saverage import saverage_sig


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
    parser.add_argument('-a', default=3, type=positive_int, help='recent CTQMC self-energies to average (default: 3, capped at floor(N/10)+1)')
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
    if solver == 'BCTQMC':
        for name in ('Sbulla.out', 'sig.out'):
            if not (job / 'imp.0' / name).is_file():
                raise ValueError('Missing BCTQMC input: ' + str(job / 'imp.0' / name))
    sigmas = []
    if solver == 'CTQMC':
        for name in glob.glob('sig.inp.*.*'):
            match = re.fullmatch(r'sig\.inp\.(\d+)\.(\d+)', name)
            if match:
                sigmas.append((tuple(map(int, match.groups())), name))
        sigmas.sort(reverse=True)
        if not sigmas:
            raise ValueError('CTQMC requires sig.inp.<iteration>.<step> files to average')
    executables = [root / name for name in ('maxent_run.py', 'dmft_copy.py', 'x_dmft.py')]
    executables.append(Path(w2k.WIENROOT) / 'x_lapw')
    if solver == 'BCTQMC':
        executables.append(root / 'Sig2header.py')
    if args.FS:
        executables.append(root / 'FSconvert.py')
    for executable in executables:
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise ValueError('Missing or non-executable program: ' + str(executable))

    print('Postprocessing solver:', solver, flush=True)
    maxent.mkdir(parents=True, exist_ok=True)
    if solver == 'CTQMC':
        selected = [name for _, name in sigmas[:min(args.a, len(sigmas)//10 + 1)]]
        print('averaging over last few sigmas:', selected, flush=True)
        saverage_sig(selected, str(maxent / 'sig.inpx'))
        sigma_input = 'sig.inpx'
    else:
        run([root / 'Sig2header.py', os.path.relpath(job / 'imp.0' / 'Sbulla.out', maxent)], maxent)
        sigma_input = 'Sbulla.inp'
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
