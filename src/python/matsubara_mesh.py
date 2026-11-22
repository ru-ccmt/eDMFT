"""Select a common lattice quadrature for the sampled impurity self-energy."""
from math import exp, log
from pathlib import Path

import numpy as np

DEFAULT_TOL = 1.e-6
HEADER = '#EDMFT_AUTO_TAIL_V1'


def impurity_nom(requested):
    """A frozen impurity's input can outlive the params.dat used to generate it."""
    counts = [] if requested is None else [int(requested)]
    for path in Path('.').glob('imp.*/PARAMS'):
        params = {}
        for line in path.read_text().splitlines():
            fields = line.partition('#')[0].split()
            if len(fields) >= 2 and fields[0] in ('nom', 'beta'):
                params[fields[0]] = fields[1]
        if 'nom' in params:
            nom = int(params['nom'])
            # The native CTQMC solver resolves nom=0 to floor(2 beta).
            if nom == 0 and 'beta' in params:
                nom = int(2 * float(params['beta']))
            counts.append(nom)
    if not counts or max(counts) < 1:
        raise ValueError('automatic tail requires --nom or a positive nom in imp.*/PARAMS')
    return max(counts)


def log_indices(n, nom, ntail):
    if nom < 1 or ntail < 1:
        raise ValueError('nom and integer ntail must be positive')
    if n <= nom + 1:
        return np.arange(n)
    ntail = min(n - nom, ntail)
    if ntail == 1:
        return np.r_[np.arange(nom), n-1]
    step = log((n-1.)/nom)/(ntail-1.)
    return np.unique(np.r_[np.arange(nom), [int(nom*exp(i*step)+.5) for i in range(ntail)]])


def endpoint(x, y):
    # Match dmft2's spline3 boundary condition: second derivative of the
    # cubic through the four nearest nodes. Rescale to avoid small powers.
    t = x - x[0]
    scale = np.max(abs(t))
    coeff = np.linalg.solve(np.vander(t / scale, 4, increasing=True), y)
    return 2 * coeff[2] / scale**2


def surrogate(omega, dynamic, broadening):
    """Scalar probes of the tail. A static level is measured relative to mu."""
    z = 1j * omega
    probes = [dynamic.real, dynamic.imag]
    # Include the historical damping floor for the diagonal probes, and
    # also monitor the original component (including off-diagonal entries).
    sigma = dynamic.real - 1j * np.maximum(abs(dynamic.imag), broadening)
    for xi in (-omega[0], 0., omega[0]):
        g0 = 1 / (z - xi)
        dg = g0**2 * sigma / (1 - g0*sigma)
        probes.extend((dg.real, (z*dg).real, -np.log1p(-g0*sigma).real))
    values = np.column_stack(probes)
    norm = np.abs(values).sum(axis=0)
    return values[:, norm > 1.e-20] / norm[norm > 1.e-20]


def refine(x, values, tolerance, seed):
    """Control the discrete sum of absolute interpolation errors per probe."""
    from scipy.interpolate import CubicSpline

    nodes = seed.copy()
    if values.shape[1] == 0:
        return nodes
    while len(nodes) < len(x):
        u, v = x[nodes][::-1], values[nodes][::-1]
        bc = ((2, endpoint(u[:4],v[:4])), (2, endpoint(u[:-5:-1],v[:-5:-1])))
        residual = abs(CubicSpline(u,v,bc_type=bc)(x) - values)
        errors = residual.sum(axis=0)
        if max(errors) <= tolerance:
            break
        # Split the interval making the largest contribution. After many
        # splits (e.g. a noisy tail) refine several disjoint intervals per
        # pass to avoid quadratic selection time. Always retain exact data.
        y = residual[:, errors.argmax()]
        prefix = np.r_[0., np.cumsum(y)]
        local = np.diff(prefix[nodes])
        intervals = np.flatnonzero(np.diff(nodes) > 1)
        batch = max(1, len(nodes)//32)
        worst = intervals[np.argsort(local[intervals])[-batch:]]
        left, right = nodes[worst], nodes[worst+1]
        med = np.searchsorted(prefix, (prefix[left]+prefix[right])/2)
        nodes = np.sort(np.r_[nodes, np.clip(med, left+1, right-1)])
    return nodes


def auto_mesh(omega, components, broadening, nom, tolerance):
    """Return lattice rows and the 1-based rows for the native moment fit.

    The 24 fit frequencies are selected independently of tolerance, from
    the same 300-point logarithmic reference used in the Ce validation.
    They travel inside the split file, including for execute-only runs.
    """
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('tail tolerance must be finite and positive')
    n = len(omega)
    if nom < 1:
        raise ValueError('nom must be positive')
    if n == 0 or omega[0] <= 0 or not np.allclose(omega / omega[0], 2*np.arange(n)+1, atol=2.e-5, rtol=1.e-9):
        raise ValueError('automatic selection needs a contiguous positive Matsubara input in sig.inp')
    nom = min(nom,n)
    fit = log_indices(n,nom,300)[-24:]
    if n-nom < 8:
        return np.arange(n), fit+1
    w = omega[nom:]
    x = 1/w
    seed = np.unique(np.rint(np.geomspace(nom,n-1,8)).astype(int) - nom)
    seed = np.union1d(seed, fit[fit>=nom]-nom)
    if len(seed) < 4:
        seed = np.arange(len(w))
    # All channels and both spins use one common mesh. Recheck every proxy
    # after forming their union: inserting knots changes a cubic spline.
    unique = {}
    for dynamic in components:
        a = np.ascontiguousarray(dynamic[nom:])
        if not np.all(np.isfinite(a)):
            raise ValueError('non-finite self-energy in Matsubara tail')
        unique.setdefault(a.tobytes(),a)
    nodes = seed
    while len(nodes)<len(w):
        before = len(nodes)
        for dynamic in unique.values():
            probes = surrogate(w,dynamic,broadening)
            nodes = refine(x,probes,tolerance,nodes)
        if len(nodes)==before:
            break
    rows = np.r_[np.arange(nom), nom+nodes]
    return rows, np.searchsorted(rows,fit)+1


def cached_auto_mesh(omega, components, broadening, nom, tolerance, cache):
    """Reuse an identical mesh during charge iterations; the fit uses live sigma."""
    import hashlib
    import os
    import tempfile
    from zipfile import BadZipFile

    # A file timestamp is insufficient: ssplit can receive different inputs
    # at the same path or the two spins can be prepared in either order.
    key = hashlib.sha256(f'v1:{nom}:{tolerance!r}:{broadening!r}'.encode())
    for a in (omega, *components):
        key.update(np.ascontiguousarray(a).view(np.uint8))
    signature = key.hexdigest()
    try:
        with cache.open('rb') as source, np.load(source, allow_pickle=False) as stored:
            if str(stored['key']) == signature:
                rows, fit = stored['rows'], stored['fit']
                if (rows.ndim == fit.ndim == 1 and rows.size and fit.size
                        and np.issubdtype(rows.dtype,np.integer) and np.issubdtype(fit.dtype,np.integer)
                        and rows[0]==0 and rows[-1]==len(omega)-1 and np.all(np.diff(rows)>0)
                        and np.array_equal(rows[:nom], np.arange(min(nom, len(omega))))
                        and fit[0]>=1 and fit[-1]==len(rows) and np.all(np.diff(fit)>0)
                        and np.array_equal(rows[fit-1], log_indices(len(omega), nom, 300)[-24:])):
                    return rows, fit
    except (OSError, KeyError, ValueError, EOFError, BadZipFile):
        pass
    rows, fit = auto_mesh(omega, components, broadening, nom, tolerance)
    # Concurrent spin preparations may share a directory. Readers always
    # see a complete old or new file. An unwritable cache only costs time.
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=cache.parent, delete=False) as out:
            temporary = Path(out.name)
            np.savez(out, key=signature, rows=rows, fit=fit)
        os.replace(temporary,cache)
    except OSError:
        pass
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return rows, fit
