"""Analytic derivative of the converged LAPW1 radial-energy search roots.

For a searched channel SELECT takes E=(E_bottom+E_top)/2, where the
physical radial slope vanishes at bottom and radial value vanishes at top.
At fixed RMT each root obeys dE=-F_V.d(rV)/F_E. This differentiates the
root continuation, not the final 1e-8-Ha staircase of native bisection.
Successful two-root searches and stable CONT branches are supported. A
missing top gives max(E_bottom,E_input); the fixed-input branch has zero
derivative, and the bottom branch has the full root derivative. Branch
crossings, HDLO or changing reference Fermi energy are not inferred.
"""

import re
from pathlib import Path
import numpy as np
from scipy.optimize import root_scalar
from lapw_radial_native import NativeRadialContext
from lapw_valence_state import _snapshot_inputs


def search_root_terms(energy, bottom, top):
    """Roots/weights of SELECT's continuous branch, using its Ry sentinels.

    SELECT writes -520 Ry for an unfound boundary. These are not eigenvalues
    and must never be sent to the radial root finder. Its CONT branch is
    exactly constant when the saved input energy is above the found bottom.
    A 0.002-Ry ambiguity margin exceeds SCF1's three-decimal root rounding.
    """
    if bottom < -260 and top < -260:
        return []
    if bottom < -260:
        raise ValueError("native SELECT has no valid bottom root")
    if top < -260:
        if energy > bottom + 0.002:
            return []
        if abs(energy - bottom) < 0.0006:
            # At the max kink this branch cannot be established from the
            # rounded SCF line alone. Require explicit policy support later.
            raise NotImplementedError(
                "CONT bottom/max branch needs the saved input energy"
            )
        raise ValueError("inconsistent native CONT branch")
    return [(bottom, 2, 0.5), (top, 0, 0.5)]


def searched_energy_derivatives(directory, case, *, complex_vectors, library):
    directory = Path(directory)
    m = _snapshot_inputs(case, directory, complex_vectors)
    pattern = r":E(\d+)_(\d+):\s*E\(\s*\d+\)=\s*(\S+)\s+E\(BOTTOM\)=\s*(\S+)\s+E\(TOP\)=\s*(\S+)"
    rows = re.findall(pattern, (directory / (case + ".scf1")).read_text())
    if not rows:
        raise ValueError("no successful native radial searches found")
    result = {}
    details = []
    for ls, as_, es, bs, ts in rows:
        l, a = int(ls), int(as_) - 1
        e0, b0, t0 = map(float, (es, bs, ts))
        root_terms = search_root_terms(e0, b0, t0)
        if not root_terms:
            continue
        s = m.s
        r = np.geomspace(s.r0[a], s.rmt[a], s.jrj[a])
        v = m.v[a][0, 0]
        kwargs = dict(
            r0=r[0],
            dx=np.log(r[-1] / r[0]) / (len(r) - 1),
            rmt=r[-1],
            lo_energies=np.array([]),
            ell=l,
            charge=s.Znuc[a],
            relativistic=s.mode.strip() == "RELA",
            lapw=True,
            quadrature="wien2k24",
            library=library,
        )
        derivatives = []
        roots = []
        for guess, component, weight in root_terms:

            def residual(energy):
                with NativeRadialContext(v, energy=energy, **kwargs) as c:
                    return float(c.primitive[-1, component, 0])

            width = 0.002
            while (
                residual(guess - width) * residual(guess + width) > 0 and width < 0.05
            ):
                width *= 2
            root = root_scalar(
                residual, bracket=[guess - width, guess + width], xtol=2e-13, rtol=1e-14
            ).root
            with NativeRadialContext(v, energy=root, **kwargs) as c:
                bar = np.zeros_like(c.primitive)
                bar[-1, component, 0] = 1
                gv, ge, _ = c.pullback(bar, np.zeros_like(c.lo_coefficients))
                if abs(ge) < 1e-12:
                    raise ValueError("singular radial root derivative")
                derivatives.append(-weight * gv / ge)
            roots.append(root)
        main = m.header[a][l]
        main = main - 200 if main > 150 else main
        local = m.lo[a][l, : m.ilo[a][l]] if l < 4 else np.zeros(0)
        energies = np.r_[main, local]
        slots = np.flatnonzero(abs(energies - e0) < 6e-5)
        if not len(slots):
            raise ValueError("searched energy does not match central vector header")
        root_energy = float(np.mean(roots))
        if max(abs(energies[slots] - root_energy)) > 2e-6:
            raise ValueError("native search root differs from vector energy")
        matrix = result.setdefault((a, l), np.zeros((len(energies), len(r))))
        for slot in slots:
            if np.any(matrix[slot]):
                raise ValueError("ambiguous repeated radial search")
            matrix[slot] = sum(derivatives)
        details.append(
            dict(
                atom=a,
                l=l,
                slots=slots.tolist(),
                roots_Ry=roots,
                native_energy_Ry=float(energies[slots[0]]),
                root_error_Ry=root_energy - float(energies[slots[0]]),
            )
        )
    return result, details
