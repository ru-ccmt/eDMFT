#!/usr/bin/env python3
# @Copyright 2025 Jia-Ming Wang
"""
init_dmft_kspace.py -- Interactive initializer for K-space-native cluster DMFT.

Generates case.indmfl + case.indmfi for a 2×2 plaquette cluster in the
K-space-native (Fourier) basis, where one block of dimension 4*N_orb spans
all cluster sites, and the cftrans encodes the plaquette Fourier × d-shell T2C.

Designed to produce the same indmfl/indmfi as the CCOC reference:
    bold_ctqmc/example/CCOC-OCA/OCAfull/CCOC_eDMFT_result/CCOC.indmfl

Usage (interactive):
    python init_dmft_kspace.py

Usage (non-interactive, for testing):
    python init_dmft_kspace.py --nointeract \
        --iatoms 17,18,19,20 --L 2 --projector 5 \
        --emin 61 --emax 204 --case CCOC

Key functions (importable):
    build_plaquette_cftrans(l=2)   -> 20×20 complex array
    build_plaquette_sigind(n_orbs, sig_K)  -> 20×20 int array
    make_plaquette_indmfl(case, iatoms, L, ...)
"""

import argparse
import os
import sys
import numpy as np

# ---------------------------------------------------------------------------
# Helper: locate the eDMFT python source so the module can be imported both
# from src/python/ and from the installed bin/ directory.
# ---------------------------------------------------------------------------
def _ensure_imports():
    here = os.path.dirname(os.path.abspath(__file__))
    for candidate in [here, os.path.join(here, '..', 'python')]:
        if os.path.isfile(os.path.join(candidate, 'indmffile.py')):
            if candidate not in sys.path:
                sys.path.insert(0, candidate)
            return
    # fall back to PYTHONPATH
_ensure_imports()

import indmffile
from cubic_harmonics import Spheric2Cubic
from inorg_edmft_io import plaquette_unitary
from gen_cix import Get_momenta_Sigind


# ---------------------------------------------------------------------------
# K-point definitions for a 2×2 plaquette
# ---------------------------------------------------------------------------
# Fractional reciprocal-space coordinates (units of 2π/a) for the 4 cluster
# momenta of a 2×2 plaquette:
#   K0 = Γ  = (0,   0)
#   K1 = X  = (0.5, 0)   [≡ (π,0) in 1st-BZ units]
#   K2 = Y  = (0,   0.5) [≡ (0,π)]
#   K3 = M  = (0.5, 0.5) [≡ (π,π)]
# The ordering matches plaquette_unitary() in inorg_edmft_io.py:
#   U[K,site] = 0.5 * exp(i K·R_site)
PLAQUETTE_Ks = np.array([[0.0, 0.0],
                          [0.5, 0.0],
                          [0.0, 0.5],
                          [0.5, 0.5]])


def build_plaquette_cftrans(l=2):
    """Build the 20×20 complex crystal-field / Fourier transformation matrix
    for a 2×2 plaquette with *l*-shell d-orbitals (l=2).

    The matrix maps the site×spherical-harmonic basis onto the K-point×
    real-harmonic basis, matching the CCOC reference (CCOC.indmfl L38-58).

    Layout in the returned matrix
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Rows (output K-real-harmonic basis):
        row r = custom_orb * N_K + K_idx
        custom orb ordering: 0=dx2y2, 1=dz2, 2=dxz, 3=dyz, 4=dxy
        K_idx: 0=Γ, 1=X(π,0), 2=Y(0,π), 3=M(π,π)

    Columns (input site×spherical-harmonic basis):
        col c = site_idx * N_sph + sph_m_idx
        site_idx: 0..3 (Cu atoms)
        sph_m_idx: 0..4 correspond to m = -2,-1,0,+1,+2

    The construction is:
        cf[K_idx, custom_orb; site, sph_m] =
            U[K_idx, site] * T2C[T2C_orb, sph_m]
    where T2C_orb = T2C_orb_for_custom[custom_orb] permutes from T2C row order
    (dz2=0, dx2y2=1, dxz=2, dyz=3, dxy=4) to the CCOC custom order that puts
    dx2y2 first.

    Parameters
    ----------
    l : int
        Angular momentum (default 2 = d-shell).

    Returns
    -------
    cf : np.ndarray, shape (N_K*(2*l+1), N_sites*(2*l+1)), dtype complex
        Crystal-field / plaquette-Fourier transformation matrix.
    """
    N_K = 4         # 2×2 plaquette
    N_sites = 4
    N_sph = 2 * l + 1   # 5 for d-shell

    # T2C: transforms spherical -> real-harmonic, shape (N_sph, N_sph)
    # Row order in Spheric2Cubic(2): [dz2, dx2y2, dxz, dyz, dxy]
    T2C = Spheric2Cubic(l)

    # Plaquette Fourier matrix U[K, site] = 0.5 * exp(i K·R_site)
    # For the 2×2 plaquette all phases are real (±1) so U is purely real.
    U = plaquette_unitary()   # shape (4, 4)

    # Standard Kronecker product: cf_kron[K*N_sph + T2C_orb, site*N_sph + sph_m]
    cf_kron = np.kron(U, T2C)   # shape (20, 20)

    # Re-order rows to match the CCOC reference convention:
    #   custom orb ordering: dx2y2 (T2C row 1), dz2 (0), dxz (2), dyz (3), dxy (4)
    #   ref row index r = custom_orb * N_K + K_idx
    # => ref[r] = cf_kron[K_idx * N_sph + T2C_orb]
    T2C_orb_for_custom = [1, 0, 2, 3, 4]   # dx2y2 first
    perm = []
    for r in range(N_K * N_sph):
        custom_orb = r // N_K
        K_idx = r % N_K
        T2C_orb = T2C_orb_for_custom[custom_orb]
        perm.append(K_idx * N_sph + T2C_orb)

    return cf_kron[perm, :]


def build_plaquette_sigind(n_orbs=20, corr_orb_idx=0, lift_xy=False):
    """Build the 20×20 diagonal Sigind matrix for the 2×2 plaquette CCOC setup.

    In the K-space basis with row ordering as in build_plaquette_cftrans:
        rows 0..3  = (dx2y2, K=Γ/X/Y/M)  <- correlated
        rows 4..19 = (dz2,dxz,dyz,dxy) × K <- not correlated (sigind=0)

    The 4 K-point components are grouped by the plaquette symmetry:
        Γ(0,0)   -> sigind = 1  (Γ-Γ component, x^2-y^2_00)
        X(π,0)   -> sigind = 2  (Γ-M component, x^2-y^2_0p)
        Y(0,π)   -> sigind = 2  (same as X by 90° rotation)
        M(π,π)   -> sigind = 3  (M-M component, x^2-y^2_pp)

    Parameters
    ----------
    n_orbs : int
        Total dimension of the sigind matrix (default 20 = 4 K-pts × 5 d-orbs).
    corr_orb_idx : int
        Index of the correlated orbital in the custom row ordering (0 = dx2y2).

    Returns
    -------
    sigind : np.ndarray, shape (n_orbs, n_orbs), dtype int
        Diagonal Sigind matrix.
    legends : list of str
        Labels for the independent self-energy components.
    """
    N_K = 4

    if lift_xy:
        # Nambu / d-wave: LIFT the X=Y degeneracy so Phi_X and Phi_Y become
        # INDEPENDENT components.  The d-wave gap requires Phi_X = -Phi_Y, which
        # collapses to zero if X and Y share a single sigind index.  K=4 distinct
        # Gamma/X/Y/M -> [1,2,3,4].  (The d-wave SIGN is imposed downstream by
        # gen_cix superc=1 / the SC cix, not by these index values.)
        sig_K = [1, 2, 3, 4]
        legends = ['x^2-y^2_G', 'x^2-y^2_X', 'x^2-y^2_Y', 'x^2-y^2_M']
    else:
        # Normal state: X=Y merged by C4 symmetry -> K=3.
        # Get symmetry grouping of K-points (same logic as gen_cix.Get_momenta_Sigind)
        sig_K_raw = Get_momenta_Sigind(PLAQUETTE_Ks)   # returns 0-indexed [0,1,1,2]
        sig_K = [s + 1 for s in sig_K_raw]             # 1-indexed: [1,2,2,3]
        legends = ['x^2-y^2_00', 'x^2-y^2_0p', 'x^2-y^2_pp']

    sigind = np.zeros((n_orbs, n_orbs), dtype=int)
    # corr_orb_idx=0 (dx2y2) occupies rows 0..3 in the custom ordering
    base_row = corr_orb_idx * N_K
    for k_idx in range(N_K):
        row = base_row + k_idx
        sigind[row, row] = sig_K[k_idx]

    return sigind, legends


def make_plaquette_indmfl(case, iatoms, L, emin, emax, projector=5,
                          matsubara=1, broadc=0.025, broadnc=0.025,
                          om_npts=200, om_emin=-3.0, om_emax=1.0, lift_xy=False):
    """Build an indmffile.Indmfl for a 2×2 plaquette cluster and write it.

    All four cluster atoms (iatoms) are gathered into ONE cix block of
    dimension 4*(2L+1).  The cftrans encodes the plaquette Fourier ⊗ T2C
    as returned by build_plaquette_cftrans().

    Parameters
    ----------
    case     : str   Wien2k case name (file prefix).
    iatoms   : list  Four atom indices (e.g. [17,18,19,20]).
    L        : int   Angular-momentum shell (2 = d-shell).
    emin,emax: int   Band indices passed to projector 5 (hybridization window).
    projector: int   Projector type (default 5).
    matsubara: int   1 = imaginary axis, 0 = real axis.
    broadc,broadnc : float   Broadenings (eV).
    om_npts  : int   Number of ω points.
    om_emin,om_emax : float  ω window (eV) for spectral plots.

    Returns
    -------
    inl : indmffile.Indmfl   The populated (and written) Indmfl object.
    """
    if len(iatoms) != 4:
        raise ValueError("build_plaquette_indmfl: need exactly 4 iatoms for 2×2 plaquette")
    qsplit = 2    # real-harmonic basis, spin-degenerate

    inl = indmffile.Indmfl(case)

    # Head parameters
    inl.emin = emin
    inl.emax = emax
    inl.Qrenorm = 1
    inl.projector = projector
    inl.matsubara = matsubara
    inl.broadc = broadc
    inl.broadnc = broadnc
    inl.om_npts = om_npts
    inl.om_emin = om_emin
    inl.om_emax = om_emax

    # Store hybr_emin/hybr_emax so write_head uses the integer band-index branch
    inl.hybr_emin = emin
    inl.hybr_emax = emax

    # Atom list: all four atoms, no local rotation, no shift
    for iatom in iatoms:
        inl.atoms[iatom] = (0, [], None, 0)   # (locrot_shift, new_xyz, shift_vec, Rmt2)
        inl.locrot[iatom] = (0, 0)

    # cix: all four atoms belong to ONE cix block (index 1)
    icix = 1
    inl.cix[icix] = [(iatom, L, qsplit) for iatom in iatoms]

    # cixgrp: one unique correlated problem
    inl.cixgrp[1] = [icix]

    # Build cftrans (20×20 complex)
    cf = build_plaquette_cftrans(l=L)

    # Build sigind + legends (20×20 diagonal)
    N_orbs = 4 * (2 * L + 1)    # 4 K-points × 5 d-orbs = 20
    sigind, legends = build_plaquette_sigind(n_orbs=N_orbs, corr_orb_idx=0, lift_xy=lift_xy)

    inl.siginds[icix] = sigind
    inl.cftrans[icix] = cf
    inl.legends[icix] = legends

    # Lsa / icpsa (needed by write_atomlist internal logic)
    for iatom in iatoms:
        inl.Lsa.append(np.array([L], dtype=int))
        inl.icpsa.append(np.array([icix], dtype=int))

    inl.write()
    return inl


def make_plaquette_indmfi(inl):
    """Write case.indmfi from the populated Indmfl object.

    Parameters
    ----------
    inl : indmffile.Indmfl   Must be already populated (cixgrp, siginds set).

    Returns
    -------
    indmfi : indmffile.Indmfi
    """
    indmfi = indmffile.Indmfi(inl)
    indmfi.write()
    return indmfi


# ---------------------------------------------------------------------------
# Validation: compare generated cftrans / sigind against CCOC reference
# ---------------------------------------------------------------------------
def validate_against_ccoc_reference(ref_indmfl_path):
    """Numeric diff against the CCOC reference indmfl.

    Returns
    -------
    max_cf_diff  : float   Max |cftrans_generated - cftrans_ref|
    sigind_match : bool    True if sigind diagonals match exactly
    legend_match : bool    True if legends match exactly
    """
    # Read reference
    ref_case = ref_indmfl_path[:-7]   # strip '.indmfl'
    ref_inl = indmffile.Indmfl(os.path.basename(ref_case))
    ref_inl.read(ref_indmfl_path)

    # Generated
    cf_gen = build_plaquette_cftrans(l=2)
    sigind_gen, legends_gen = build_plaquette_sigind(n_orbs=20, corr_orb_idx=0)

    # cftrans diff
    cf_ref = ref_inl.cftrans[1]
    max_cf_diff = np.max(np.abs(cf_gen - cf_ref))

    # Sigind
    sig_ref_diag = list(np.diag(ref_inl.siginds[1]))
    sig_gen_diag = list(np.diag(sigind_gen))
    sigind_match = (sig_gen_diag == sig_ref_diag)

    # Legends
    legend_match = (legends_gen == ref_inl.legends[1])

    return max_cf_diff, sigind_match, legend_match, sig_gen_diag, sig_ref_diag


# ---------------------------------------------------------------------------
# Interactive / non-interactive main
# ---------------------------------------------------------------------------
def _detect_case():
    """Return Wien2k case name from current directory."""
    import glob
    structs = glob.glob('*.struct')
    if structs:
        return structs[0][:-7]
    # fall back to directory name
    return os.path.basename(os.getcwd())


def _ask(prompt, default):
    """Ask user a question; return default if empty input."""
    try:
        ans = input(prompt + f" [{default}]: ").strip()
    except (EOFError, KeyboardInterrupt):
        return default
    return ans if ans else default


def main():
    parser = argparse.ArgumentParser(
        description='Generate K-space plaquette indmfl/indmfi for 2×2 cluster DMFT')
    parser.add_argument('--nointeract', action='store_true',
                        help='Non-interactive mode: use defaults / argv')
    parser.add_argument('--case', default=None,
                        help='Wien2k case name (default: auto-detect)')
    parser.add_argument('--iatoms', default='17,18,19,20',
                        help='Four Cu iatom indices, comma-separated')
    parser.add_argument('--L', type=int, default=2,
                        help='Angular momentum (2=d-shell)')
    parser.add_argument('--projector', type=int, default=5,
                        help='Projector type')
    parser.add_argument('--emin', type=int, default=61,
                        help='Lower band index for projector 5')
    parser.add_argument('--emax', type=int, default=204,
                        help='Upper band index for projector 5')
    parser.add_argument('--matsubara', type=int, default=1,
                        help='1=imaginary axis, 0=real axis')
    parser.add_argument('--broadc', type=float, default=0.025)
    parser.add_argument('--broadnc', type=float, default=0.025)
    parser.add_argument('--om_npts', type=int, default=200)
    parser.add_argument('--om_emin', type=float, default=-3.0)
    parser.add_argument('--om_emax', type=float, default=1.0)
    parser.add_argument('--validate', default=None,
                        help='Path to reference .indmfl to validate against')
    parser.add_argument('--nambu', action='store_true',
                        help='Generate a K=4 (X≠Y lifted) normal indmfl for the '
                             'Nambu/d-wave workflow; feed to nambu_indmfl_gen to double it')
    args = parser.parse_args()

    nointeract = args.nointeract

    # Detect case
    case = args.case if args.case else _detect_case()
    if not nointeract:
        case = _ask('Wien2k case name', case)

    print(f"\nInitializing K-space 2×2 plaquette DMFT for case: {case}")

    # Iatoms
    iatoms_str = args.iatoms
    if not nointeract:
        iatoms_str = _ask('Four plaquette iatom indices (comma-separated)', iatoms_str)
    iatoms = [int(x.strip()) for x in iatoms_str.split(',')]
    if len(iatoms) != 4:
        print(f"ERROR: need exactly 4 iatoms, got {len(iatoms)}: {iatoms}")
        sys.exit(1)
    print(f"  Cluster atoms (iatoms): {iatoms}")

    # L
    L = args.L
    if not nointeract:
        L = int(_ask('Angular momentum shell (2=d)', str(L)))

    # Projector
    projector = args.projector
    if not nointeract:
        projector = int(_ask('Projector type', str(projector)))

    # Band range (for projector 5)
    emin = args.emin
    emax = args.emax
    if not nointeract:
        emin = int(_ask('Band index emin', str(emin)))
        emax = int(_ask('Band index emax', str(emax)))

    # Matsubara
    matsubara = args.matsubara
    if not nointeract:
        a_str = _ask('Imaginary (i) or real (r) axis', 'i' if matsubara else 'r')
        matsubara = 1 if a_str.strip().lower() != 'r' else 0

    broadc   = args.broadc
    broadnc  = args.broadnc
    om_npts  = args.om_npts
    om_emin  = args.om_emin
    om_emax  = args.om_emax

    if not nointeract:
        broadc  = float(_ask('Broadening correlated (eV)', str(broadc)))
        broadnc = float(_ask('Broadening noncorr (eV)', str(broadnc)))
        om_npts = int(_ask('Number of ω points', str(om_npts)))
        om_emin = float(_ask('ω_min (eV)', str(om_emin)))
        om_emax = float(_ask('ω_max (eV)', str(om_emax)))

    # -------------------------------------------------------------------
    # Generate and write case.indmfl
    # -------------------------------------------------------------------
    print(f"\nBuilding 20×20 plaquette cftrans (4 K-pts × 5 d-orbs)...")
    if args.nambu:
        print("  --nambu: lifting X=Y degeneracy -> K=4 [1,2,3,4] (for d-wave Nambu)")
    inl = make_plaquette_indmfl(
        case, iatoms, L, emin, emax,
        projector=projector, matsubara=matsubara,
        broadc=broadc, broadnc=broadnc,
        om_npts=om_npts, om_emin=om_emin, om_emax=om_emax, lift_xy=args.nambu)
    print(f"  Written: {case}.indmfl")

    # -------------------------------------------------------------------
    # Generate and write case.indmfi
    # -------------------------------------------------------------------
    indmfi = make_plaquette_indmfi(inl)
    print(f"  Written: {case}.indmfi")

    print(f"\nSigind diagonal (first 5): {list(np.diag(inl.siginds[1])[:5])}")
    print(f"Legends: {inl.legends[1]}")
    max_cf = np.max(np.abs(inl.cftrans[1]))
    print(f"cftrans max |entry|: {max_cf:.8f}  (expected ~0.5 for dyz, ~0.35355339 for dx2y2)")

    # -------------------------------------------------------------------
    # Optional validation against reference
    # -------------------------------------------------------------------
    if args.validate:
        ref_path = args.validate
        if not os.path.isfile(ref_path):
            print(f"\nWARNING: reference file not found: {ref_path}")
        else:
            print(f"\nValidating against reference: {ref_path}")
            max_cf_diff, sigind_match, legend_match, sig_gen_diag, sig_ref_diag = \
                validate_against_ccoc_reference(ref_path)
            print(f"  cftrans max |diff|: {max_cf_diff:.2e}")
            print(f"  Sigind diag match:  {sigind_match}")
            print(f"    generated: {sig_gen_diag}")
            print(f"    reference: {sig_ref_diag}")
            print(f"  Legends match:      {legend_match}")
            if max_cf_diff <= 1e-8 and sigind_match and legend_match:
                print("  VALIDATION PASSED (all checks ≤ 1e-8)")
            else:
                print("  VALIDATION FAILED")
                sys.exit(1)

    print("\nDone.")
    return inl, indmfi


if __name__ == '__main__':
    main()
