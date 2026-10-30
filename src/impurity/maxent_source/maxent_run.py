#!/usr/bin/env python
"""Continue Sigma through matrix or scalar auxiliary cumulants.

Usage: python maxent_run.py sig.inp --indmfi UTe2.indmfi --outdir continued
The input with # s_oo / # Edc headers contains dynamic Sigma; output Sig.out
retains that convention and the headers. Headerless input is split using the
last Matsubara point's real part, as in the original script.
Only maxent_run.py changes; the scalar maxentropy.py backend is reused.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import argparse
import ast
from concurrent.futures import ProcessPoolExecutor
import importlib
import json
from pathlib import Path
import runpy
import sys
import numpy as np
from scipy import interpolate, optimize


def read_sigma(path):
    headers = []
    metadata = {}
    with Path(path).open() as stream:
        for line in stream:
            if not line.lstrip().startswith('#'):
                break
            headers.append(line.rstrip())
            text = line.lstrip()[1:].strip()
            if '=' in text:
                key, value = text.split('=', 1)
                if key.strip() in ('s_oo', 'Edc'):
                    metadata[key.strip()] = np.asarray(ast.literal_eval(value.strip()), float)
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] < 3 or data.shape[1] % 2 != 1 or not np.isfinite(data).all():
        raise ValueError('Self-energy requires finite frequency, real/imaginary column pairs')
    omega = data[:, 0]
    if len(omega) < 8 or omega[0] <= 0 or not np.all(np.diff(omega) > 0):
        raise ValueError('At least eight increasing positive Matsubara frequencies are required')
    beta = np.pi / omega[0]
    if not np.allclose(omega, (2*np.arange(len(omega))+1)*np.pi/beta, rtol=2e-9, atol=1e-10):
        raise ValueError('InverseFourier requires a complete positive fermionic Matsubara mesh')
    sigma = data[:, 1::2] + 1j*data[:, 2::2]
    dynamic_input = 's_oo' in metadata
    if dynamic_input:
        for name, values in metadata.items():
            if values.shape != (sigma.shape[1],) or not np.isfinite(values).all():
                raise ValueError(name+' header does not match self-energy columns')
        static = metadata['s_oo']
    else:
        if metadata:
            raise ValueError('An Edc header without s_oo has an ambiguous static convention')
        static = sigma[-1].real.copy()
        sigma = sigma-static
    return omega, sigma, static, headers, dynamic_input, beta


def find_indmfi(input_path, explicit=None):
    """Prefer a single mapping beside input, then in its parent directory.

    Keep the directory of an input symlink: maxent/sig.inp may point to
    ../sig.inp while a local maxent/*.indmfi should still take precedence.
    """
    if explicit is not None:
        path = Path(explicit).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError('Explicit indmfi file does not exist: '+str(path))
        return str(path)
    directory = Path(input_path).expanduser().absolute().parent.resolve()
    for folder in (directory, directory.parent):
        candidates = sorted(p for p in folder.glob('*.indmfi') if p.is_file())
        if len(candidates) > 1:
            raise ValueError('Multiple indmfi files in '+str(folder)+'; choose one with --indmfi: '
                             +', '.join(p.name for p in candidates))
        if candidates:
            return str(candidates[0].resolve())
    return None


def read_indmfi(path, ncolumns):
    if path is None:
        # Scalar legacy input, one independent diagonal orbital per column.
        return [np.diag(np.arange(1, ncolumns+1))]
    lines = []
    for line in Path(path).read_text().splitlines():
        text = line.split('#', 1)[0].strip()
        if text:
            lines.append(text)
    count = int(lines[0]); cursor = 1; blocks = []
    for _ in range(count):
        dim = int(lines[cursor]); cursor += 1
        matrix = np.array([[int(v) for v in row.split()] for row in lines[cursor:cursor+dim]], int)
        cursor += dim
        if matrix.shape != (dim, dim) or np.any(matrix < 0) or np.any(matrix > ncolumns):
            raise ValueError('indmfi needs square blocks with labels 0..number of Sigma columns')
        if np.any((np.diag(matrix) == 0) & (np.any(matrix != 0, axis=1) | np.any(matrix != 0, axis=0))):
            raise ValueError('A coupled orbital needs a nonzero diagonal indmfi label')
        blocks.append(matrix)
    if cursor != len(lines):
        raise ValueError('Unexpected trailing indmfi data')
    used = {int(x) for block in blocks for x in block.flat if x > 0}
    if used != set(range(1, ncolumns+1)):
        raise ValueError('indmfi must account for every input component, including measured zeros')
    return blocks


def connected_blocks(mappings):
    """Factor disconnected orbital groups and reuse identically labelled blocks."""
    unique = []; seen = {}; physical = []
    for block_id, mapping in enumerate(mappings):
        unseen = set(np.flatnonzero(np.diag(mapping) > 0).tolist())
        while unseen:
            start = min(unseen); unseen.remove(start); group = [start]; todo = [start]
            while todo:
                i = todo.pop()
                neighbors = {j for j in unseen if mapping[i, j] or mapping[j, i]}
                unseen -= neighbors; group.extend(sorted(neighbors)); todo.extend(sorted(neighbors))
            group.sort(key=lambda i: (int(mapping[i, i]), i))
            local = mapping[np.ix_(group, group)]
            key = (local.shape, tuple(local.flat))
            if key not in seen:
                seen[key] = len(unique); unique.append(local)
            physical.append({'indmfi_block': block_id, 'orbitals': group, 'unique_block': seen[key]})
    return unique, physical


def expand(values, mapping):
    matrix = np.zeros((len(values), len(mapping), len(mapping)), complex)
    for i, j in zip(*np.nonzero(mapping)):
        matrix[:, i, j] = values[:, mapping[i, j]-1]
    return matrix


def choose_shift(omega, sigma, symmetrize):
    """Keep the old scalar SymCum shifts, then use ONE diagonal shift matrix.

    Auxiliary projections must never be symmetrized independently: they must
    all be projections of the same matrix cumulant.
    """
    n = sigma.shape[1]
    if not symmetrize:
        return np.zeros(n), {'enabled': False}
    stop = min(len(omega), max(8, int(np.searchsorted(omega, 5.0))))
    w = omega[:stop]; shift = []; residuals = []
    for i in range(n):
        diagonal = sigma[:stop, i, i]
        def objective(c):
            m = 1./(1j*w+c-diagonal)
            value = interpolate.UnivariateSpline(w, m.real, s=0)(0.)
            derivative = interpolate.UnivariateSpline(w, (-m*m).real, s=0)(0.)
            return value, derivative
        initial = float(diagonal[0].real)
        result = optimize.root_scalar(objective, x0=initial, fprime=True, method='newton')
        root = float(result.root)
        # The objective also approaches zero at infinity; reject runaway roots.
        bound = max(10., 2.*float(np.max(abs(diagonal))))
        error = float(abs(objective(root)[0]))
        if not result.converged or not np.isfinite(root) or abs(root-initial) > bound or error > 1e-6:
            raise ValueError('Scalar SymCum shift optimization failed; disable SymCum to use zero shift')
        shift.append(root); residuals.append(error)
    return np.array(shift), {'enabled': True, 'scheme': 'legacy scalar orbital shifts, shared diagonal C for every auxiliary',
                             'maximum_scalar_residual': max(residuals)}


def projection_vectors(size, symmetric):
    """Minimal positive projections: diagonal, real sum, and phase sum.

    A plain minus projection only repeats information about Fij+Fji.
    The relative phase -i supplies the independent Fij-Fji combination.
    """
    identity = np.eye(size, dtype=complex)
    for i in range(size):
        yield ('d', i, i), identity[i]
    for i in range(size):
        for j in range(i+1, size):
            yield ('rp', i, j), (identity[i]+identity[j])/np.sqrt(2.)
            if not symmetric:
                yield ('ip', i, j), (identity[i]-1j*identity[j])/np.sqrt(2.)


def assemble_projected(results, size):
    """Polarization identity for spectra OR complex analytic functions.

    At the same complex frequency Fji need not equal conj(Fij).
    Only the spectral matrix is Hermitian; do not conjugate the functions.
    """
    first = next(iter(results.values()))
    matrix = np.zeros((len(first), size, size), complex)
    for i in range(size):
        matrix[:, i, i] = results[('d', i, i)]
    for i in range(size):
        for j in range(i+1, size):
            diagonal = .5*(results[('d', i, i)]+results[('d', j, j)])
            real_part = results[('rp', i, j)]-diagonal
            phase_part = results[('ip', i, j)]-diagonal if ('ip', i, j) in results else 0.
            matrix[:, i, j] = real_part+1j*phase_part
            matrix[:, j, i] = real_part-1j*phase_part
    return matrix


def matrix_spectrum(results, size):
    return assemble_projected(results, size)


def projected_input(omega, sigma, shift, vector, method, cumulant=None):
    projected_shift = float(np.real(np.vdot(vector, shift*vector)))
    if method == 'matrix':
        if cumulant is None:
            cumulant = np.linalg.inv(1j*omega[:, None, None]*np.eye(len(shift))+np.diag(shift)-sigma)
        projected = np.einsum('i,wij,j->w', vector.conj(), cumulant, vector)
    elif method == 'scalar':
        projected_sigma = np.einsum('i,wij,j->w', vector.conj(), sigma, vector)
        projected = 1./(1j*omega+projected_shift-projected_sigma)
    else:
        raise ValueError('cumulant_method must be matrix or scalar')
    return projected, projected_shift


def integration_weights(omega):
    return np.r_[.5*(omega[1]-omega[0]), .5*(omega[2:]-omega[:-2]), .5*(omega[-1]-omega[-2])]


def constrain_spectrum(spectrum, omega, mode):
    """Enforce matrix positivity and the cumulant's unit zeroth moment.

    This is an explicit post-fit correction, not a joint matrix-MaxEnt fit.
    Its magnitude and the resulting back-fit errors are recorded.
    """
    weights = integration_weights(omega)
    eigenvalues, vectors = np.linalg.eigh(spectrum)
    minimum = float(eigenvalues.min())
    if mode == 'check' and minimum < -1e-10:
        raise ValueError('Independent auxiliary fits give a nonpositive matrix spectrum; use --matrix-psd project and inspect the back-fit errors')
    positive = np.einsum('wij,wj,wkj->wik', vectors, np.maximum(eigenvalues, 0.), vectors.conj())
    moment = np.einsum('w,wij->ij', weights, positive)
    norm_eigenvalues, norm_vectors = np.linalg.eigh(moment)
    if norm_eigenvalues.min() <= 1e-12:
        raise ValueError('Matrix spectrum has a singular integrated weight')
    whitening = (norm_vectors/np.sqrt(norm_eigenvalues)) @ norm_vectors.conj().T
    corrected = whitening @ positive @ whitening
    relative_change = np.sqrt(np.sum(weights[:, None, None]*abs(corrected-spectrum)**2)/max(np.sum(weights[:, None, None]*abs(spectrum)**2), 1e-300))
    return corrected, {'minimum_raw_eigenvalue': minimum, 'relative_spectral_correction': float(relative_change),
                       'moment_eigenvalues_before_normalization': norm_eigenvalues.tolist(),
                       'unit_moment_error': float(np.max(abs(np.einsum('w,wij->ij', weights, corrected)-np.eye(spectrum.shape[1]))))}


def retarded_cumulant(omega, spectrum, backend):
    hilbert = np.empty_like(spectrum, dtype=complex)
    for i in range(spectrum.shape[1]):
        for j in range(spectrum.shape[2]):
            a = spectrum[:, i, j]
            hilbert[:, i, j] = -np.pi*(backend.KramarsKronig(omega, a.real)+1j*backend.KramarsKronig(omega, a.imag))
    return hilbert-1j*np.pi*spectrum


def continue_task(payload):
    # Every task has a separate directory: maxentropy.py writes intermediate files.
    task, params, directory = payload
    old = Path.cwd(); directory = Path(directory); directory.mkdir(exist_ok=True)
    try:
        os.chdir(directory)
        backend = importlib.import_module('maxentropy')
        with open('maxent.log', 'w') as log:
            spectrum, omega = backend.MaximumEntropy(params, task['tau'], task['gt'], log=log)
        return task['block'], task['key'], np.asarray(spectrum), np.asarray(omega)
    finally:
        os.chdir(old)


def tau_backfit(omega, scalar_spectrum, task):
    weights = integration_weights(omega)
    tau = task['tau']; beta = tau[-1]
    kernel = np.exp(-tau[:, None]*omega[None, :]-np.logaddexp(0., -beta*omega)[None, :])
    fit = -kernel @ (scalar_spectrum*weights)
    error = fit-task['gt']
    return fit, {'tau_rms_error': float(np.sqrt(np.mean(error**2))), 'tau_max_error': float(np.max(abs(error)))}


def finish_block(omega_real, spectra, block, block_tasks, om_check, backend, method, psd_mode):
    n = len(block['mapping']); weights = integration_weights(omega_real)
    kernel = weights/(1j*om_check[:, None]-omega_real[None, :])
    info = {}; artifacts = {}
    if method == 'matrix':
        raw = matrix_spectrum(spectra, n)
        a, info = constrain_spectrum(raw, omega_real, psd_mode)
        cumulant = retarded_cumulant(omega_real, a, backend)
        continued = omega_real[:, None, None]*np.eye(n)+np.diag(block['shift'])-np.linalg.inv(cumulant)
        back = np.einsum('nw,wij->nij', kernel, a)
        sigma_back = 1j*om_check[:, None, None]*np.eye(n)+np.diag(block['shift'])-np.linalg.inv(back)
        first_moment_offset = np.einsum('w,wij->ij', weights*omega_real, a)+np.diag(block['shift'])
        fit_spectra = {task['key']: np.einsum('i,wij,j->w', task['vector'].conj(), a, task['vector']).real for task in block_tasks}
        artifacts.update(spectral_matrix=a, raw_spectral_matrix=raw, cumulant=cumulant)
    else:
        # Each scalar cumulant has unit mass. No joint matrix-PSD repair is
        # applied after nonlinear scalar inversions; report Sigma causality.
        fit_spectra = {}; functions = {}; back_functions = {}; offsets = {}; masses = []
        for task in block_tasks:
            key = task['key']; raw = spectra[key]; mass = float(raw @ weights)
            if mass <= 0:
                raise ValueError('A scalar cumulant spectrum has nonpositive mass')
            a = raw/mass; fit_spectra[key] = a; masses.append(mass)
            retarded = retarded_cumulant(omega_real, a[:, None, None], backend)[:, 0, 0]
            functions[key] = omega_real+task['projected_shift']-1./retarded
            back_functions[key] = 1j*om_check+task['projected_shift']-1./(kernel @ a)
            offsets[key] = np.array([task['projected_shift']+np.dot(weights*omega_real, a)])
        continued = assemble_projected(functions, n)
        sigma_back = assemble_projected(back_functions, n)
        first_moment_offset = assemble_projected(offsets, n)[0]
        info = {'relative_spectral_correction': float(max(abs(1./np.array(masses)-1.))),
                'scalar_masses_before_normalization': masses,
                'matrix_positivity_enforced': False}
        artifacts.update(scalar_spectra=np.array([fit_spectra[t['key']] for t in block_tasks]),
                         projection_keys=np.array([t['key'] for t in block_tasks], dtype=str))
    return continued, sigma_back, first_moment_offset, fit_spectra, info, artifacts


def write_sigma(path, omega, sigma, headers):
    data = np.empty((len(omega), 1+2*sigma.shape[1]))
    data[:, 0] = omega; data[:, 1::2] = sigma.real; data[:, 2::2] = sigma.imag
    with Path(path).open('w') as stream:
        if headers:
            stream.write('\n'.join(headers)+'\n')
        np.savetxt(stream, data, fmt='%.15g')


def format_diagnostics(report):
    """Human-readable companion to the complete machine-readable JSON."""
    def number(value):
        return "n/a" if value is None else f"{value:.5g}"
    blocks = report.get('continuation', [])
    projections = [p for block in blocks for p in block.get('projections', [])]
    method = report.get('cumulant_method', 'unknown')
    minimum = report.get('final_sigma_minimum_decay_eigenvalue')
    causal = report.get('final_sigma_causality_passed', minimum >= -1e-8 if minimum is not None else None)
    status = 'PASS' if causal else ('FAIL' if causal is not None else 'not checked')
    representation = report.get('label_representation_max_error')
    mapping_status = 'not checked' if representation is None else ('PASS' if representation <= 1e-6 else 'FAIL')
    lines = [
        'ANALYTIC CONTINUATION DIAGNOSTICS',
        '',
        f"Method:       {method} cumulant",
        f"Input:        {report.get('input', 'unknown')}",
        f"Mapping:      {report.get('indmfi') or 'none (independent diagonal columns)'}",
        f"Beta:         {number(report.get('beta'))} eV^-1",
        f"Components:   {report.get('ncomponents', '?')}    Unique blocks: {report.get('nunique_blocks', '?')}    Scalar continuations: {report.get('nprojections', '?')}",
        '',
        f"Matrix causality: {status}",
        f"  Smallest eigenvalue of -(Sigma-Sigma^dagger)/(2i): {number(minimum)} eV",
        '  Values below -1e-8 eV fail this check.',
        f"Component mapping: {mapping_status} (maximum representation error: {number(representation)} eV)",
    ]
    if method == 'matrix':
        lines.append('Matrix spectral positivity and unit normalization: '+str(report.get('matrix_psd_mode', 'unknown')))
    else:
        lines.append('Joint matrix positivity is checked, but not repaired, in scalar mode.')
    correction = max((b.get('relative_spectral_correction', 0.) for b in blocks), default=0.)
    label = 'PSD/normalization spectral correction' if method == 'matrix' else 'scalar normalization correction'
    lines.append(f'Largest {label}: {100*correction:.5g}%')
    ratio = max((p.get('rms_over_deltag', 0.) for p in projections), default=0.)
    lines.extend([
        f'Largest imaginary-time fit RMS / deltag: {ratio:.5g}',
        '',
        'SELF-ENERGY BACK-CONTINUATION ERRORS (eV)',
        'Block  Components             First iw    iw < 1 eV   iw < 10 eV  Tail offset',
        '-----  --------------------  ----------  ----------  ----------  -----------',
    ])
    mappings = report.get('blocks', [])
    for block in blocks:
        bid = block['block']
        mapping = mappings[bid].get('mapping', []) if bid < len(mappings) else []
        labels = ','.join(str(v) for v in sorted({v for row in mapping for v in row if v > 0}))
        lines.append(f"{bid:5d}  {labels:20s}  {number(block.get('self_energy_back_error_first_frequency')):>10s}  "
                     f"{number(block.get('self_energy_back_error_below_1eV')):>10s}  "
                     f"{number(block.get('self_energy_back_error_below_10eV')):>10s}  "
                     f"{number(block.get('inferred_dynamic_static_offset_max')):>11s}")
    lines.extend([
        '',
        'Errors are maxima over checked Matsubara points in each window.',
        'They compare the input self-energy with the back-continued fit.',
        'They measure reconstruction accuracy; they are not statistical error bars.',
        'Tail offset is the inferred constant part of dynamic Sigma at infinity.',
        'The first spectral moment is not constrained, and this offset is not removed.',
        '',
        'Full details: diagnostics.json',
        'Individual imaginary-time fits: b*/gt_backfit.dat',
    ])
    return '\n'.join(lines)+'\n'


def run(args, comm=None):
    rank = comm.Get_rank() if comm is not None else 0
    size = comm.Get_size() if comm is not None else 1
    backend = importlib.import_module('maxentropy')
    outdir = Path(args.outdir).resolve()
    if rank == 0:
        outdir.mkdir(parents=True, exist_ok=True)
        omega, sigma, static, headers, dynamic_input, beta = read_sigma(args.input)
        params = runpy.run_path(args.params)['params']
        method = args.cumulant or params.get('cumulant_method', 'matrix')
        if method not in ('matrix', 'scalar'):
            raise ValueError('cumulant_method must be matrix or scalar')
        if params['statistics'] != 'fermi':
            raise ValueError('Matrix cumulants need fermionic statistics')
        indmfi = find_indmfi(args.input, args.indmfi)
        if indmfi is not None:
            print('Using indmfi mapping: '+indmfi, flush=True)
        mappings = read_indmfi(indmfi, sigma.shape[1])
        unique, physical = connected_blocks(mappings)
        tau = np.linspace(0., beta, int(params['Ntau'])+1)
        tasks = []; blocks = []
        for b, mapping in enumerate(unique):
            s = expand(sigma, mapping)
            decay = -(s-s.conj().transpose(0, 2, 1))/(2j)
            causal_min = float(np.linalg.eigvalsh(decay).min())
            if causal_min < -1e-8:
                raise ValueError(f'Input Sigma block {b} violates matrix causality: minimum decay eigenvalue {causal_min}')
            shift, shift_info = choose_shift(omega, s, params.get('SymCum', False))
            cumulant = None
            inputs = {'omega': omega, 'sigma_dynamic': s, 'mapping': mapping, 'shift': shift, 'method': method}
            if method == 'matrix':
                cumulant = np.linalg.inv(1j*omega[:, None, None]*np.eye(len(mapping))+np.diag(shift)-s)
                inputs['cumulant'] = cumulant
            blocks.append({'mapping': mapping, 'shift': shift, 'shift_info': shift_info, 'input_decay_min': causal_min})
            np.savez(outdir/f'cumulant_input.block{b}.npz', **inputs)
            for key, vector in projection_vectors(len(mapping), np.array_equal(mapping, mapping.T)):
                projected, projected_shift = projected_input(omega, s, shift, vector, method, cumulant)
                if np.max(projected.imag) > 1e-9:
                    raise ValueError('A projected cumulant is not causal')
                gt = backend.InverseFourier(projected, omega, tau, beta, min(int(params['Nf']), len(omega)))
                normalization = float(-gt[0]-gt[-1])
                if abs(normalization-1.) > .02 or np.max(gt) > max(1e-6, float(params['deltag'])):
                    raise ValueError('Inverse Fourier cumulant is not a positive-spectrum fermionic function')
                name = f'b{b}_{key[0]}_{key[1]}_{key[2]}'
                folder = outdir/name; folder.mkdir(exist_ok=True)
                np.savetxt(folder/'gw0.dat', np.c_[omega, projected.real, projected.imag])
                np.savetxt(folder/'gt0.dat', np.c_[tau, gt])
                tasks.append({'block': b, 'key': key, 'vector': vector, 'gt': gt, 'tau': tau, 'name': name,
                              'inverse_fourier_normalization': normalization, 'projected_shift': projected_shift})
        manifest = {'input': str(Path(args.input).resolve()), 'indmfi': indmfi, 'beta': beta,
                    'dynamic_input': dynamic_input, 'cumulant_method': method, 'shift_rule': 'C_v = v^dagger diag(C_i) v; no independent auxiliary optimization', 'ncomponents': sigma.shape[1], 'nphysical_blocks': len(physical),
                    'nunique_blocks': len(blocks), 'nprojections': len(tasks), 'physical_blocks': physical,
                    'blocks': [{'mapping': b['mapping'].tolist(), 'shift': b['shift'].tolist(), 'shift_info': b['shift_info'],
                                'input_decay_min': b['input_decay_min']} for b in blocks]}
        (outdir/'preparation.json').write_text(json.dumps(manifest, indent=2)+'\n')
        print(f'Prepared {method} cumulants: {len(blocks)} unique blocks and {len(tasks)} positive projections; beta={beta:.12g}', flush=True)
        bundle = (params, tasks)
    else:
        bundle = None
    if comm is not None:
        bundle = comm.bcast(bundle, root=0)
    params, tasks = bundle
    if args.prepare_only:
        return
    payloads = [(task, params, str(outdir/task['name'])) for task in tasks[rank::size]]
    if args.workers > 1:
        if size > 1:
            raise ValueError('Use MPI ranks or --workers, not both')
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(continue_task, payloads))
    else:
        results = [continue_task(payload) for payload in payloads]
    if comm is not None:
        gathered = comm.gather(results, root=0)
        if rank == 0:
            results = [item for chunk in gathered for item in chunk]
    if rank != 0:
        return
    omega_real = results[0][3]
    for _, _, a, w in results:
        if not np.array_equal(w, omega_real) or not np.isfinite(a).all() or np.min(a) < -1e-12:
            raise ValueError('Scalar MaxEnt spectra must be finite, nonnegative, and on a common mesh')
    # Match the original driver's convention of omitting the mesh's zero point.
    keep = np.abs(omega_real) >= 1e-6
    omega_real = omega_real[keep]
    component_sum = np.zeros((len(omega_real), sigma.shape[1]), complex); counts = np.zeros(sigma.shape[1])
    diagnostics = dict(manifest)
    diagnostics['matrix_psd_mode'] = args.matrix_psd if method == 'matrix' else 'not_enforced_for_scalar'
    diagnostics['continuation'] = []
    for b, block in enumerate(blocks):
        spectra = {key: a[keep] for block_id, key, a, _ in results if block_id == b}
        block_tasks = [t for t in tasks if t['block'] == b]
        pick = np.unique(np.r_[np.arange(min(200, len(omega))),
                               np.geomspace(min(200, len(omega)-1), len(omega)-1, 100).astype(int)])
        om_check = omega[pick]
        continued, sigma_back, first_moment_offset, fit_spectra, constraint_info, artifacts = finish_block(
            omega_real, spectra, block, block_tasks, om_check, backend, method, args.matrix_psd)
        entry = {'block': b, **constraint_info, 'projections': []}
        for task in block_tasks:
            fit, info = tau_backfit(omega_real, fit_spectra[task['key']], task)
            np.savetxt(outdir/task['name']/'gt_backfit.dat', np.c_[task['tau'], task['gt'], fit])
            entry['projections'].append({'name': task['name'], **info, 'rms_over_deltag': info['tau_rms_error']/float(params['deltag']),
                                         'inverse_fourier_normalization': task['inverse_fourier_normalization'],
                                         'projected_shift': task['projected_shift']})
        target = expand(sigma[pick], block['mapping'])
        error = np.abs(sigma_back-target)
        entry['self_energy_back_error_first_frequency'] = float(error[0].max())
        entry['self_energy_back_error_below_1eV'] = float(error[om_check < 1.].max()) if np.any(om_check < 1.) else None
        entry['self_energy_back_error_below_10eV'] = float(error[om_check < 10.].max()) if np.any(om_check < 10.) else None
        entry['inferred_dynamic_static_offset_max'] = float(np.abs(first_moment_offset).max())
        np.savez(outdir/f'continued.block{b}.npz', omega=omega_real,
                 sigma_dynamic=continued, mapping=block['mapping'], shift=block['shift'], method=method, **artifacts)
        for i, j in zip(*np.nonzero(block['mapping'])):
            label = block['mapping'][i, j]-1
            component_sum[:, label] += continued[:, i, j]; counts[label] += 1
        diagnostics['continuation'].append(entry)
    reduced = component_sum/counts
    if not dynamic_input:
        reduced = reduced+static
    write_sigma(outdir/'Sig.out', omega_real, reduced, headers if dynamic_input else [])
    # Recheck the final labelled representation after symmetry averaging.
    minimum = np.inf; representation_error = 0.
    for b, block in enumerate(blocks):
        reconstructed = expand(reduced, block['mapping'])
        minimum = min(minimum, float(np.linalg.eigvalsh(-(reconstructed-reconstructed.conj().transpose(0, 2, 1))/(2j)).min()))
        saved = np.load(outdir/f'continued.block{b}.npz')['sigma_dynamic']
        if not dynamic_input:
            saved = saved+expand(static[None, :], block['mapping'])
        representation_error = max(representation_error, float(np.max(abs(saved-reconstructed))))
    diagnostics['final_sigma_minimum_decay_eigenvalue'] = minimum
    diagnostics['final_sigma_causality_passed'] = bool(minimum >= -1e-8)
    diagnostics['label_representation_max_error'] = representation_error
    (outdir/'diagnostics.json').write_text(json.dumps(diagnostics, indent=2)+'\n')
    (outdir/'diagnostics.txt').write_text(format_diagnostics(diagnostics))
    worst_fit = max(x['rms_over_deltag'] for b in diagnostics['continuation'] for x in b['projections'])
    max_correction = max(b['relative_spectral_correction'] for b in diagnostics['continuation'])
    print(f'Wrote {outdir / "Sig.out"}; largest spectral correction={max_correction:.4g}, worst tau RMS/deltag={worst_fit:.4g}', flush=True)
    worst_first = max(b['self_energy_back_error_first_frequency'] for b in diagnostics['continuation'])
    print(f'Largest first-Matsubara Sigma back-error={worst_first:.5g} eV; see diagnostics.txt for the summary and diagnostics.json for full details.', flush=True)
    if minimum < -1e-8 or representation_error > 1e-6 or worst_fit > 2.:
        print('Continuation needs review: inspect diagnostics.txt and gt_backfit.dat before using Sig.out.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input')
    parser.add_argument('--indmfi', help='default: find one *.indmfi beside input, then in its parent directory')
    parser.add_argument('--params', default='maxent_params.dat')
    parser.add_argument('--outdir', default='.')
    parser.add_argument('--workers', type=int, default=1, help='independent local MaxEnt processes')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--cumulant', choices=['matrix', 'scalar'], help='default: params[cumulant_method], or matrix')
    parser.add_argument('--matrix-psd', choices=['project', 'check'], default='project',
                        help='matrix method only; scalar method reports Sigma causality without repairing it')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be positive')
    comm = None
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
    except ImportError:
        if int(os.environ.get('OMPI_COMM_WORLD_SIZE', os.environ.get('PMI_SIZE', '1'))) > 1:
            raise RuntimeError('mpi4py is required for an MPI launch; use ordinary Python for serial/local workers')
    try:
        run(args, comm)
    except Exception:
        if comm is not None and comm.Get_size() > 1:
            import traceback
            traceback.print_exc(); comm.Abort(1)
        raise


if __name__ == '__main__':
    main()
