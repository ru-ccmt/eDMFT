#!/usr/bin/env python3
"""
norg_edmft_io.py — Format conversion bridge between eDMFT framework and iNORG solver.

coded by Jia-Ming Wang (jiaming.w@rutgers.edu, Rutgers University, U.S.) date 2025 - now

Usage:
    # Preprocessing: convert eDMFT files → iNORG input files
    python norg_edmft_io.py preprocess [--edinput input.tsv] [--params PARAMS.norg] [--output-dir .]

    # Postprocessing: convert iNORG output (K-basis) → eDMFT format (site-basis)
    python norg_edmft_io.py postprocess [--sigind Trans.dat] [--nband 4] [--output-dir .]

Functions can also be imported and called from inorg.py.
"""

import numpy as np
import os
import sys
import re
import argparse


# =============================================================================
# iNORG Green function file format
# =============================================================================
#
# Header line:
#   w  output_Re11  output_Re12 ... output_ReNN  w  output_Im11  output_Im12 ... output_ImNN
# Data lines (one per Matsubara frequency):
#   omega  Re(g[0][0]) Re(g[0][1]) ... Re(g[N-1][N-1])  omega  Im(g[0][0]) Im(g[0][1]) ... Im(g[N-1][N-1])
#
# - First column: Matsubara frequency omega_n = (2n+1)*pi/beta
# - Next norbs^2 columns: real parts of g[m][n] (row-major)
# - Repeat: frequency + norbs^2 imaginary parts
# - Scientific notation, space-separated


# =============================================================================
# input.tsv format (written by eDMFT's prepare_EDinput)
# =============================================================================
#
# Tij
# <norb x norb real matrix — impurity energy levels>
#
# Delta
# <one line per Matsubara frequency>
#   (omega+0j) (Delta[0,0]) (Delta[0,1]) ... (Delta[norb-1,norb-1])
#   Each value is Python complex format: (real+imagj)
#
# Vij_abba
# <nband_half x nband_half real matrix — density-density interaction>
#
# Vij_abab
# <nband_half x nband_half real matrix — pair-hopping/spin-flip interaction>


# =============================================================================
# K-basis transformation utilities (2x2 plaquette)
# =============================================================================

def plaquette_unitary():
    """Build 4x4 unitary matrix for 2x2 plaquette site->K transform.

    Site positions: r0=(0,0), r1=(1,0), r2=(0,1), r3=(1,1)
    K-points: K0=(0,0), K1=(pi,0), K2=(0,pi), K3=(pi,pi)

    U[K,r] = (1/2) * exp(i K.r)
    """
    return 0.5 * np.array([
        [1,  1,  1,  1],
        [1, -1,  1, -1],
        [1,  1, -1, -1],
        [1, -1, -1,  1],
    ], dtype=float)


def transform_matrix_to_K(mat, U):
    """T^K = U @ T^site @ U^dagger"""
    return U @ mat @ U.conj().T


def transform_green_to_K(green, U):
    """Delta^K(iw) = U @ Delta^site(iw) @ U^dagger for each frequency.

    Args:
        green: (nomg, norb, norb) complex array
        U: (norb, norb) unitary matrix

    Returns:
        (nomg, norb, norb) complex array in K-basis
    """
    Ud = U.conj().T
    return np.einsum('ij,njk,kl->nil', U, green, Ud)


def transform_matrix_to_site(mat_K, U):
    """T^site = U^dagger @ T^K @ U (inverse transform)"""
    return U.conj().T @ mat_K @ U


def transform_green_to_site(green_K, U):
    """G^site(iw) = U^dagger @ G^K(iw) @ U for each frequency.

    Args:
        green_K: (nomg, norb, norb) complex array in K-basis
        U: (norb, norb) unitary matrix

    Returns:
        (nomg, norb, norb) complex array in site-basis
    """
    Ud = U.conj().T
    return np.einsum('ij,njk,kl->nil', Ud, green_K, U)


def transform_interaction_to_K_vijkl(vij_abba, vij_abab, U, nband, nambu=False):
    """Transform Ising interaction from site-basis to K-basis Vijkl format.

    For Nambu representation (nambu=True):
      After Nambu transformation (c+_{b,dn} -> psi_{b+nband}, c_{b,dn} -> psi+_{b+nband}):
      - abba term: (a, b+nband, b+nband, a) with value -V_abba[a][b]
      - abab term: (a, a+nband, b+nband, b) with value -V_abab[a][b]

    For normal state (nambu=False):
      Interaction +U n_up n_dn in two-spin representation:
      - abba term: (a, b+nband, b+nband, a) with value +V_abba[a][b]
      - abab term: (a, a+nband, b+nband, b) with value +V_abab[a][b]

    K-basis transformation of the full 4-index tensor:
    V^K_{IJKL} = sum_{ijkl} U*[I,i] U*[J,j] U[K,k] U[L,l] V^site_{ijkl}

    Returns list of (I, J, K, L, value) tuples in K-basis.
    """
    tol = 1e-12
    vijkl = []
    nband_half = vij_abba.shape[0]

    # Sign factor: -1 for Nambu (hole convention), +1 for normal state
    sign = -1.0 if nambu else +1.0

    norbs = 2 * nband
    vijkl_tensor = np.zeros((norbs, norbs, norbs, norbs))

    for a in range(nband_half):
        for b in range(nband_half):
            # abba contribution
            if abs(vij_abba[a][b]) > tol:
                val = sign * vij_abba[a][b]
                for I in range(nband):
                    for L in range(nband):
                        coeff_IL = U[I, a].conj() * U[L, a]
                        if abs(coeff_IL) < tol:
                            continue
                        for J in range(nband):
                            for K in range(nband):
                                coeff_JK = U[J, b].conj() * U[K, b]
                                if abs(coeff_JK) < tol:
                                    continue
                                vijkl_tensor[I][J + nband][K + nband][L] += val * coeff_IL * coeff_JK

            # abab contribution (a != b)
            if a != b and abs(vij_abab[a][b]) > tol:
                val = sign * vij_abab[a][b]
                for I in range(nband):
                    for J in range(nband):
                        coeff_IJ = U[I, a].conj() * U[J, a].conj()
                        if abs(coeff_IJ) < tol:
                            continue
                        for K in range(nband):
                            for L in range(nband):
                                coeff_KL = U[K, b] * U[L, b]
                                if abs(coeff_KL) < tol:
                                    continue
                                vijkl_tensor[I][J + nband][K + nband][L] += val * coeff_IJ * coeff_KL

    # Collect non-zero entries
    for I in range(norbs):
        for J in range(norbs):
            for K in range(norbs):
                for L in range(norbs):
                    if abs(vijkl_tensor[I][J][K][L]) > tol:
                        vijkl.append((I, J, K, L, vijkl_tensor[I][J][K][L]))

    return vijkl


def compute_nambu_shift(vij_abba):
    """Compute per-orbital half-fill shift from Vij_abba.

    shift[a] = (1/2) * sum_b Vij_abba[a][b]

    For uniform U case (diagonal Vij_abba): shift[a] = U/2
    """
    return 0.5 * np.sum(vij_abba, axis=1)


def parse_complex(s):
    """Parse a Python-style complex number string like '(1.5+2.3j)' or '(1.5-2.3j)'."""
    s = s.strip()
    if s.startswith('(') and s.endswith(')'):
        s = s[1:-1]
    return complex(s)


def read_edinput(filename):
    """Read an input.tsv file written by eDMFT's prepare_EDinput().

    Returns:
        tij: numpy array (norb, norb) — impurity energy levels
        omega: numpy array (nomg,) — Matsubara frequencies (real)
        delta: numpy array (nomg, norb, norb) — complex hybridization function
        vij_abba: numpy array (nband_half, nband_half) — density-density interaction
        vij_abab: numpy array (nband_half, nband_half) — pair-hopping interaction
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find section boundaries
    sections = {}
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped in ('Tij', 'Delta', 'Vij_abba', 'Vij_abab', 'Vijkl'):
            sections[stripped] = i

    # --- Read Tij ---
    tij_start = sections['Tij'] + 1
    tij_lines = []
    for i in range(tij_start, len(lines)):
        stripped = lines[i].strip()
        if not stripped or stripped in sections:
            break
        tij_lines.append(stripped)
    tij = np.array([[float(x) for x in row.split()] for row in tij_lines])
    norb = tij.shape[0]

    # --- Read Delta ---
    delta_start = sections['Delta'] + 1
    delta_lines = []
    for i in range(delta_start, len(lines)):
        stripped = lines[i].strip()
        if not stripped or stripped in sections:
            break
        delta_lines.append(stripped)

    nomg = len(delta_lines)
    omega = np.zeros(nomg)
    delta = np.zeros((nomg, norb, norb), dtype=complex)

    for i, line in enumerate(delta_lines):
        # Split by ')' then find complex numbers
        parts = re.findall(r'\([^)]+\)', line)
        if len(parts) != norb * norb + 1:
            raise ValueError(
                f"Line {i}: expected {norb*norb+1} complex values, got {len(parts)}")
        # First value is (omega + 0j)
        omega[i] = parse_complex(parts[0]).imag  # omega_n is the imaginary part
        # if omega is real-valued (stored as real+0j), extract real part
        omega_cmplx = parse_complex(parts[0])
        if omega_cmplx.imag == 0:
            omega[i] = omega_cmplx.real
        else:
            omega[i] = omega_cmplx.imag

        # Remaining norb*norb values are Delta matrix elements (row-major)
        idx = 1
        for m in range(norb):
            for n in range(norb):
                delta[i, m, n] = parse_complex(parts[idx])
                idx += 1

    # --- Read Vij_abba ---
    vij_abba = None
    vij_abab = None
    vijkl = None

    if 'Vij_abba' in sections:
        vij_start = sections['Vij_abba'] + 1
        vij_lines = []
        for i in range(vij_start, len(lines)):
            stripped = lines[i].strip()
            if not stripped or stripped in sections:
                break
            vij_lines.append(stripped)
        vij_abba = np.array([[float(x) for x in row.split()] for row in vij_lines])

    if 'Vij_abab' in sections:
        vij_start = sections['Vij_abab'] + 1
        vij_lines = []
        for i in range(vij_start, len(lines)):
            stripped = lines[i].strip()
            if not stripped or stripped in sections:
                break
            vij_lines.append(stripped)
        vij_abab = np.array([[float(x) for x in row.split()] for row in vij_lines])

    if 'Vijkl' in sections:
        vijkl_start = sections['Vijkl'] + 1
        vijkl = []
        for i in range(vijkl_start, len(lines)):
            stripped = lines[i].strip()
            if not stripped or stripped in sections:
                break
            parts = stripped.split()
            if len(parts) == 5:
                vijkl.append((int(parts[0]), int(parts[1]),
                              int(parts[2]), int(parts[3]),
                              float(parts[4])))

    return tij, omega, delta, vij_abba, vij_abab, vijkl


def write_inorg_green(filename, omega, g_complex, norbs, label="output"):
    """Write a Green function file in iNORG format.

    Args:
        filename: output file path
        omega: array of Matsubara frequencies (real values, i.e., Im(z))
        g_complex: array (nomg, norbs, norbs) of complex Green function values
        norbs: number of orbitals
        label: prefix for column headers (e.g., "output" or "999")
    """
    nomg = len(omega)
    with open(filename, 'w') as f:
        # Header
        header_parts = [f"{'w':>24s}"]
        for m in range(norbs):
            for n in range(norbs):
                header_parts.append(f"{label+'_Re'+str(m+1)+str(n+1):>24s}")
        header_parts.append(f"{'w':>24s}")
        for m in range(norbs):
            for n in range(norbs):
                header_parts.append(f"{label+'_Im'+str(m+1)+str(n+1):>24s}")
        f.write("  ".join(header_parts) + "\n")

        # Data
        for i in range(nomg):
            parts = [f"{omega[i]:24.16e}"]
            for m in range(norbs):
                for n in range(norbs):
                    parts.append(f"{g_complex[i, m, n].real:24.16e}")
            parts.append(f"{omega[i]:24.16e}")
            for m in range(norbs):
                for n in range(norbs):
                    parts.append(f"{g_complex[i, m, n].imag:24.16e}")
            f.write("  ".join(parts) + "\n")


def read_inorg_green(filename, norbs=None):
    """Read a Green function file in iNORG format.

    Args:
        filename: input file path
        norbs: number of orbitals (auto-detected from header if None)

    Returns:
        omega: array of Matsubara frequencies
        g_complex: array (nomg, norbs, norbs) of complex Green function values
    """
    with open(filename, 'r') as f:
        # Read and parse header to determine norbs
        header = f.readline()
        if norbs is None:
            # Count columns: 1 (w) + norbs^2 (Re) + 1 (w) + norbs^2 (Im)
            ncols = len(header.split())
            # ncols = 2 + 2*norbs^2  =>  norbs^2 = (ncols - 2) / 2
            norbs_sq = (ncols - 2) // 2
            norbs = int(np.sqrt(norbs_sq))
            assert norbs * norbs == norbs_sq, \
                f"Cannot determine norbs from {ncols} columns"

        # Read data
        data = np.loadtxt(f)
        nomg = data.shape[0]
        ncols_expected = 2 + 2 * norbs * norbs

        omega = data[:, 0]
        g_complex = np.zeros((nomg, norbs, norbs), dtype=complex)

        for i in range(nomg):
            idx_re = 1
            idx_im = 1 + norbs * norbs + 1
            for m in range(norbs):
                for n in range(norbs):
                    re_val = data[i, idx_re]
                    im_val = data[i, idx_im]
                    g_complex[i, m, n] = complex(re_val, im_val)
                    idx_re += 1
                    idx_im += 1

    return omega, g_complex


def write_interaction_file(filename, vij_abba=None, vij_abab=None, vijkl=None):
    """Write the interaction tensor file for iNORG.

    Args:
        filename: output file path
        vij_abba: (nband_half, nband_half) array for Ising format
        vij_abab: (nband_half, nband_half) array for Ising format
        vijkl: list of (i, j, k, l, value) tuples for Full Coulomb format
    """
    with open(filename, 'w') as f:
        if vijkl is not None:
            # Full Coulomb format
            f.write("Vijkl\n")
            for (i, j, k, l, value) in vijkl:
                f.write(f"{i} {j} {k} {l} {value:.16e}\n")
        elif vij_abba is not None:
            # Ising format
            f.write("Vij_abba\n")
            for row in vij_abba:
                f.write(" ".join(f"{v:.16e}" for v in row) + "\n")
            f.write("Vij_abab\n")
            if vij_abab is not None:
                for row in vij_abab:
                    f.write(" ".join(f"{v:.16e}" for v in row) + "\n")
            else:
                nband_half = vij_abba.shape[0]
                zeros = np.zeros((nband_half, nband_half))
                for row in zeros:
                    f.write(" ".join(f"{v:.16e}" for v in row) + "\n")


def read_sigind_matrix(filename):
    """Read sigind_matrix from Trans.dat file.

    Args:
        filename: path to Trans.dat

    Returns:
        sigind: numpy array (ntot, ntot) — sigind matrix
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find the Sigind section
    sigind_start = None
    matrix_size = None
    for i, line in enumerate(lines):
        if '#---- Sigind follows' in line:
            sigind_start = i + 1
        elif line.strip() and sigind_start is None:
            # First non-empty line may have matrix size
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                matrix_size = int(parts[0])

    if sigind_start is None:
        # No marker found, try reading from beginning (old format)
        sigind_start = 0
        if matrix_size is not None:
            sigind_start = 1  # skip size line

    # Read matrix
    rows = []
    for i in range(sigind_start, len(lines)):
        line = lines[i].strip()
        if not line or line.startswith('#'):
            if rows:
                break
            continue
        try:
            row = [int(x) for x in line.split()]
            rows.append(row)
        except ValueError:
            break

    if not rows:
        raise ValueError(f"Could not read sigind matrix from {filename}")

    return np.array(rows)


def build_sigind_for_test(nband):
    """Build a simple diagonal sigind_matrix for testing.

    For nband orbitals without spin polarization:
    sigind_matrix is (2*nband, 2*nband) with diagonal elements
    [1, 2, ..., nband, 1, 2, ..., nband] (spin-up then spin-down share same indices).

    For iNORG's nband (which equals norb in eDMFT), the sigind is norb x norb
    with diagonal elements [1, 2, ..., nband].
    """
    # Simple diagonal sigind for nband orbitals
    norb = nband  # In eDMFT, norb = nband (correlated orbitals per spin if spin_dep=False)
    sigind = np.zeros((norb, norb), dtype=int)
    for i in range(norb):
        sigind[i, i] = i + 1
    return sigind


def expand_delta_from_edinput(delta_compact, norb):
    """Expand eDMFT compact Delta (norb x norb) to iNORG norbs (= 2*norb, up+dn).

    In Nambu representation for iNORG:
    - Indices 0..norb-1: spin-up orbitals (particle sector)
    - Indices norb..2*norb-1: spin-down orbitals (hole sector)

    For non-spin-dependent eDMFT, Delta_up = Delta_dn = delta_compact.
    In Nambu: Delta_full is block-diagonal: diag(Delta_up, -Delta_dn^*)

    Actually for the eDMFT interface, we pass norbs = 2*nband to iNORG where nband=norb.
    The input.tsv provides Delta for the correlated orbitals (norb x norb).
    For if_mat_type=1, norg_sets=2 (up/dn blocks), so norbs=2*nband and
    the hybridization is stored as a norbs x norbs = (2*norb) x (2*norb) matrix.

    In this case, we place delta_compact in both spin blocks:
    Delta_full[0:norb, 0:norb] = delta_compact   (spin-up)
    Delta_full[norb:, norb:]   = delta_compact    (spin-down, same as up for non-polarized)
    """
    norbs = 2 * norb
    nomg = delta_compact.shape[0]
    delta_full = np.zeros((nomg, norbs, norbs), dtype=complex)
    # Spin-up block
    delta_full[:, :norb, :norb] = delta_compact
    # Spin-down block (same as up for non-spin-polarized)
    delta_full[:, norb:, norb:] = delta_compact
    return delta_full


def expand_tij_from_edinput(tij_compact):
    """Expand eDMFT compact Tij (norb x norb) to iNORG norbs (= 2*norb).

    Places Tij in both spin blocks:
    Tij_full[0:norb, 0:norb] = tij_compact   (spin-up)
    Tij_full[norb:, norb:]   = tij_compact    (spin-down)
    """
    norb = tij_compact.shape[0]
    norbs = 2 * norb
    tij_full = np.zeros((norbs, norbs))
    tij_full[:norb, :norb] = tij_compact
    tij_full[norb:, norb:] = tij_compact
    return tij_full


def generate_prmtr_txt(output_dir, norbs, beta, unit_omg, num_omg, sym="psc"):
    """Generate a minimal prmtr.txt for the iNORG solver.

    This provides default values that will be overridden by APIedmft::update().
    The key parameters (nband, norbs, norg_sets, etc.) come from PARAMS.norg.
    """
    prmtr_path = os.path.join(output_dir, "prmtr.txt")
    with open(prmtr_path, 'w') as f:
        f.write(f"{'nband':20s} {norbs}\n")
        f.write(f"{'unit_omg':20s} {unit_omg:.16e}\n")
        f.write(f"{'max_omg':20s} {unit_omg * num_omg * 2:.16e}\n")
        f.write(f"{'hubbU':20s} 0.0\n")
        f.write(f"{'mu':20s} 0.0\n")
        f.write(f"{'norg_sets':20s} {norbs}\n")
        f.write(f"{'sym':20s} {sym}\n")
        f.write(f"{'ndiv':20s} 1\n")
        f.write(f"{'iter_max':20s} 1\n")
        f.write(f"{'iter_min':20s} 1\n")
    return prmtr_path


def preprocess_edinput(edinput_file, params_file=None, output_dir=".", sym="psc"):
    """Preprocess eDMFT input.tsv file to generate iNORG input files.

    This is the main preprocessing function when input.tsv is available.
    It reads the combined eDMFT input file and generates:
      - input_hb.txt           (hybridization in iNORG Green format, K-basis)
      - input_tij.txt          (on-site energy matrix, K-basis)
      - input_interaction.txt  (interaction tensor, K-basis Vijkl format)

    For Nambu (sym='psc'/etc.): includes Nambu half-fill shift in Tij, -V in interaction.
    For normal state (sym='normal-pm'/etc.): no shift, +V in interaction.

    The site->K basis transformation is performed here in Python preprocessing.
    The C++ solver works entirely in K-basis.

    Args:
        edinput_file: path to input.tsv
        params_file: path to PARAMS.norg (optional, for additional settings)
        output_dir: directory to write output files
        sym: symmetry type for determining Nambu vs normal state
    """
    print(f"Reading input.tsv from: {edinput_file}")

    tij, omega, delta_compact, vij_abba, vij_abab, vijkl = read_edinput(edinput_file)
    norb = tij.shape[0]
    norbs = 2 * norb  # iNORG uses norbs = 2*nband (Nambu / up+dn)
    nomg = len(omega)

    is_nambu = not str(sym).startswith('normal')

    print(f"  norb (per spin) = {norb}")
    print(f"  norbs (total)   = {norbs}")
    print(f"  nomg            = {nomg}")
    print(f"  sym             = {sym} (nambu={is_nambu})")
    print(f"  omega range     = [{omega[0]:.6f}, {omega[-1]:.6f}]")
    if vij_abba is not None:
        print(f"  Vij_abba shape  = {vij_abba.shape}")
        print(f"  Vij_abba diag   = {np.diag(vij_abba)}")
    if vijkl is not None:
        print(f"  Vijkl entries   = {len(vijkl)}")

    # Build unitary matrix for site->K transform
    U = plaquette_unitary()
    print(f"  Using 2x2 plaquette unitary (4x4)")

    # 1. Compute Nambu half-fill shift and add to Tij diagonal (Nambu only)
    tij_shifted = tij.copy()
    if is_nambu and vij_abba is not None:
        shift = compute_nambu_shift(vij_abba)
        print(f"  Nambu shift     = {shift}")
        for a in range(norb):
            tij_shifted[a, a] += shift[a]
    elif not is_nambu:
        print(f"  Normal state: no Nambu shift applied")

    # 2. Transform Tij to K-basis
    tij_K = transform_matrix_to_K(tij_shifted, U)
    print(f"  Tij (K-basis) diagonal = {np.diag(tij_K)}")

    # 3. Transform Delta to K-basis
    delta_K = transform_green_to_K(delta_compact, U)
    print(f"  Delta transformed to K-basis")

    # 4. Expand Delta to full norbs x norbs (spin-up + spin-down blocks)
    delta_full = expand_delta_from_edinput(delta_K, norb)
    print(f"  Delta expanded to shape: {delta_full.shape}")

    # 5. Transform interaction to K-basis Vijkl format
    if vij_abba is not None:
        vijkl_K = transform_interaction_to_K_vijkl(
            vij_abba, vij_abab if vij_abab is not None else np.zeros_like(vij_abba),
            U, norb, nambu=is_nambu)
        print(f"  Interaction: {len(vijkl_K)} Vijkl entries in K-basis (nambu={is_nambu})")
    elif vijkl is not None:
        # Already in Vijkl format, pass through (assumed already K-basis or no transform needed)
        vijkl_K = vijkl
        print(f"  Interaction: {len(vijkl_K)} Vijkl entries (pass-through)")
    else:
        vijkl_K = None
        print(f"  Warning: no interaction data found")

    # Write input_tij.txt (K-basis, nband x nband, includes Nambu shift)
    tij_path = os.path.join(output_dir, "input_tij.txt")
    np.savetxt(tij_path, tij_K, fmt='%.16e')
    print(f"  Written: {tij_path}")

    # Write input_hb.txt (iNORG Green format, K-basis)
    hb_path = os.path.join(output_dir, "input_hb.txt")
    write_inorg_green(hb_path, omega, delta_full, norbs, label="output")
    print(f"  Written: {hb_path}")

    # Write input_interaction.txt (K-basis Vijkl format)
    interaction_path = os.path.join(output_dir, "input_interaction.txt")
    write_interaction_file(interaction_path, vijkl=vijkl_K)
    print(f"  Written: {interaction_path}")

    # Compute beta and unit_omg from omega spacing
    if nomg >= 2:
        domg = omega[1] - omega[0]
        # omega_n = (2n+1)*pi/beta, so domg = 2*pi/beta
        beta = 2 * np.pi / domg
        unit_omg = np.pi / beta
        print(f"  Derived: beta = {beta:.6f}, unit_omg = {unit_omg:.10e}")

    # Generate minimal prmtr.txt if needed
    prmtr_path = os.path.join(output_dir, "prmtr.txt")
    if not os.path.exists(prmtr_path):
        generate_prmtr_txt(output_dir, norbs, beta, unit_omg, nomg)
        print(f"  Generated: {prmtr_path}")

    print("Preprocessing complete.")
    return {
        'norb': norb,
        'norbs': norbs,
        'nomg': nomg,
        'beta': beta,
        'omega': omega,
        'tij': tij,
        'tij_K': tij_K,
        'delta': delta_full,
    }


def write_trans_dat(filename, sigind):
    """Write a Trans.dat file with the given sigind matrix."""
    n = sigind.shape[0]
    with open(filename, 'w') as f:
        f.write(f"{n} {n} #  size of Sigind and CF\n")
        f.write("#---- Sigind follows\n")
        for i in range(n):
            row = "  ".join(f"{sigind[i,j]:4d}" for j in range(n))
            f.write(f"  {row}\n")
        f.write("#---- CF follows\n")
        # Write identity CF (crystal field = identity)
        for i in range(n):
            row = []
            for j in range(n):
                if i == j:
                    row.append("  1.00000000   0.00000000")
                else:
                    row.append("  0.00000000   0.00000000")
            f.write("".join(row) + "\n")


def postprocess_from_norg(output_dir=".", sigind=None, nband=None,
                          seimp_file=None, gfimp_file=None, sym="psc",
                          transform_to_site=True):
    """Convert iNORG output to eDMFT sigind-compressed format.

    iNORG output is in K-basis (norbs x norbs). This function:
    1. Reads the norbs x norbs output
    2. Extracts the nband x nband spin-up block (indices 0:nband)
    3. Inverse K->site transforms: Sigma^site(iw) = U^dagger @ Sigma^K(iw) @ U
    4. Applies sigind compression and writes Sig.out/Gf.out

    Reads:
      - io/seimp.txt or seimp_file  (self-energy, K-basis)
      - io/gfimp.txt or gfimp_file  (Green's function, K-basis)

    Writes:
      - Sig.out  (eDMFT sigind format, site-basis)
      - Gf.out   (eDMFT sigind format, site-basis)

    Args:
        output_dir: directory containing iNORG output
        sigind: sigind_matrix (numpy array) or None (will build simple diagonal)
        nband: number of bands (per spin), required if sigind is None
        seimp_file: explicit path to self-energy file
        gfimp_file: explicit path to Green's function file
    """
    # Build inverse transform matrix
    U = plaquette_unitary()

    # Find self-energy file
    if seimp_file is None:
        candidates = [
            os.path.join(output_dir, "io", "seimp.txt"),
            os.path.join(output_dir, "seimp.txt"),
        ]
        # Also check for zicXXX.mb.seimp.txt pattern
        io_dir = os.path.join(output_dir, "io")
        if os.path.isdir(io_dir):
            for fn in sorted(os.listdir(io_dir), reverse=True):
                if fn.endswith(".mb.seimp.txt"):
                    candidates.insert(0, os.path.join(io_dir, fn))
        for c in candidates:
            if os.path.exists(c):
                seimp_file = c
                break

    if gfimp_file is None:
        candidates = [
            os.path.join(output_dir, "io", "gfimp.txt"),
            os.path.join(output_dir, "gfimp.txt"),
        ]
        io_dir = os.path.join(output_dir, "io")
        if os.path.isdir(io_dir):
            for fn in sorted(os.listdir(io_dir), reverse=True):
                if fn.endswith(".mb.gfimp.txt"):
                    candidates.insert(0, os.path.join(io_dir, fn))
        for c in candidates:
            if os.path.exists(c):
                gfimp_file = c
                break

    # Read the Green functions (K-basis, norbs x norbs)
    if seimp_file and os.path.exists(seimp_file):
        print(f"Reading self-energy from: {seimp_file}")
        omega, seimp = read_inorg_green(seimp_file)
        norbs = seimp.shape[1]
        norb = norbs // 2  # iNORG norbs = 2*norb (Nambu)

        if nband is None:
            nband = norb
        if sigind is None:
            sigind = build_sigind_for_test(norb)

        # Extract spin block(s) and inverse K->site transform
        is_normal = str(sym).startswith('normal')
        if is_normal:
            # Normal state: average spin-up and spin-down blocks
            seimp_up = seimp[:, :norb, :norb]
            seimp_dn = seimp[:, norb:, norb:]
            seimp_K_block = (seimp_up + seimp_dn) / 2
            print(f"  Averaged spin-up and spin-down blocks (normal state)")
        else:
            # Nambu: only use spin-up (particle sector)
            seimp_K_block = seimp[:, :norb, :norb]
        if transform_to_site:
            seimp_site = transform_green_to_site(seimp_K_block, U)
            print(f"  Inverse K->site transform applied to self-energy")
        else:
            seimp_site = seimp_K_block
            print(f"  K-native: self-energy kept in K-basis (no K->site rotation)")

        # Write Sig.out (site-basis)
        sig_out_path = os.path.join(output_dir, "Sig.out")
        write_edmft_green(sig_out_path, omega, seimp_site, sigind, norb,
                          header_type="sigma")
        print(f"  Written: {sig_out_path}")

        # Extract and write anomalous self-energy block (Nambu mode only)
        if not is_normal and seimp.shape[1] > norb:
            # Nambu mode detected: norbs = 2*norb, anomalous block exists
            seimp_anomal_K = seimp[:, :norb, norb:]  # Anomalous block (off-diagonal)
            if transform_to_site:
                seimp_anomal_site = transform_green_to_site(seimp_anomal_K, U)
                print(f"  Inverse K->site transform applied to anomalous self-energy")
            else:
                seimp_anomal_site = seimp_anomal_K
                print(f"  K-native: anomalous self-energy kept in K-basis")

            # Write Phi.out (site-basis, same format as Sig.out)
            phi_out_path = os.path.join(output_dir, "Phi.out")
            write_edmft_green(phi_out_path, omega, seimp_anomal_site, sigind, norb,
                              header_type="sigma")
            print(f"  Written: {phi_out_path}")
    else:
        print("Warning: self-energy file not found")

    if gfimp_file and os.path.exists(gfimp_file):
        print(f"Reading Green's function from: {gfimp_file}")
        omega, gfimp = read_inorg_green(gfimp_file)
        norbs = gfimp.shape[1]
        norb = norbs // 2

        if nband is None:
            nband = norb
        if sigind is None:
            sigind = build_sigind_for_test(norb)

        # Extract spin block(s) and inverse K->site transform
        is_normal = str(sym).startswith('normal')
        if is_normal:
            # Normal state: average spin-up and spin-down blocks
            gfimp_up = gfimp[:, :norb, :norb]
            gfimp_dn = gfimp[:, norb:, norb:]
            gfimp_K_block = (gfimp_up + gfimp_dn) / 2
            print(f"  Averaged spin-up and spin-down blocks (normal state)")
        else:
            # Nambu: only use spin-up (particle sector)
            gfimp_K_block = gfimp[:, :norb, :norb]
        if transform_to_site:
            gfimp_site = transform_green_to_site(gfimp_K_block, U)
            print(f"  Inverse K->site transform applied to Green's function")
        else:
            gfimp_site = gfimp_K_block
            print(f"  K-native: Green's function kept in K-basis (no K->site rotation)")

        # Write Gf.out (site-basis)
        gf_out_path = os.path.join(output_dir, "Gf.out")
        write_edmft_green(gf_out_path, omega, gfimp_site, sigind, norb,
                          header_type="green")
        print(f"  Written: {gf_out_path}")

        # Extract and write anomalous Green function F (Nambu mode only)
        if not is_normal and gfimp.shape[1] > norb:
            gfimp_anomal_K = gfimp[:, :norb, norb:]  # off-diagonal block
            if transform_to_site:
                gfimp_anomal_site = transform_green_to_site(gfimp_anomal_K, U)
                print(f"  Inverse K->site transform applied to anomalous Green function")
            else:
                gfimp_anomal_site = gfimp_anomal_K
                print(f"  K-native: anomalous Green function kept in K-basis")

            f_out_path = os.path.join(output_dir, "F.out")
            write_edmft_green(f_out_path, omega, gfimp_anomal_site, sigind, norb,
                              header_type="green")
            print(f"  Written: {f_out_path}")
    else:
        print("Warning: Green's function file not found")

    print("Postprocessing complete.")


def write_edmft_green(filename, omega, g_complex, sigind, norb,
                      header_type="sigma"):
    """Write Green function in eDMFT sigind-compressed format.

    Format:
        # s_oo= [val1, val2, ..., valN]     (for Sig.out only)
        # Edc= [val1, val2, ..., valN]      (for Sig.out only)
        omega_0  re_1 im_1  re_2 im_2  ...  re_M im_M
        ...

    The sigind matrix maps (i,j) → index k. For each unique k, we average
    all g[i][j] where sigind[i][j] == k, then output as (re_k, im_k).

    For iNORG's Nambu representation, norbs = 2*norb. We extract the
    spin-up block g[0:norb, 0:norb] for the eDMFT output.

    Args:
        filename: output file path
        omega: Matsubara frequencies
        g_complex: (nomg, norbs, norbs) complex Green function
        sigind: (norb, norb) sigind matrix
        norb: number of orbitals per spin
        header_type: "sigma" or "green"
    """
    nomg = len(omega)
    # Slice sigind to match data dimensions — in Nambu mode, the full sigind
    # is 40×40 (max=12) but the data is only norb×norb (4×4, pp block).
    # Without slicing, max_idx=12 creates zero-padded components 4-12 that
    # break the merge→sgather pipeline (sgather reads 12 pairs from a 48-pair file).
    sigind = np.array(sigind)[:norb, :norb]
    max_idx = int(np.max(sigind))

    with open(filename, 'w') as f:
        if header_type == "sigma":
            # Write s_oo and Edc headers (high-frequency limit)
            # s_oo = Re(Sigma(i*inf)) for each sigind component
            # For now, extract from the last Matsubara frequency
            s_oo = np.zeros(max_idx)
            edc = np.zeros(max_idx)
            for k in range(1, max_idx + 1):
                vals = []
                for i in range(norb):
                    for j in range(norb):
                        if sigind[i, j] == k:
                            vals.append(g_complex[-1, i, j].real)
                if vals:
                    s_oo[k - 1] = np.mean(vals)
            f.write("# s_oo= " + str(list(s_oo)) + "\n")
            f.write("# Edc= " + str(list(edc)) + "\n")

        # Write data: omega  re_1 im_1  re_2 im_2 ...
        for n in range(nomg):
            parts = [f"{omega[n]:.17e}"]
            for k in range(1, max_idx + 1):
                vals = []
                for i in range(norb):
                    for j in range(norb):
                        if sigind[i, j] == k:
                            vals.append(g_complex[n, i, j])
                if vals:
                    avg = np.mean(vals)
                    parts.append(f"{avg.real:.17e}")
                    parts.append(f"{avg.imag:.17e}")
                else:
                    parts.append("0.0")
                    parts.append("0.0")
            f.write("  ".join(parts) + "\n")


def read_edmft_green(filename, sigind, norb):
    """Read eDMFT sigind-compressed Green function and expand to full matrix.

    Args:
        filename: input file path (e.g., Sig.out, Delta.inp)
        sigind: (norb, norb) sigind matrix
        norb: number of orbitals

    Returns:
        omega: array of frequencies
        g_complex: (nomg, norb, norb) complex Green function
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Skip comment/header lines
    data_lines = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        data_lines.append(stripped)

    nomg = len(data_lines)
    max_idx = int(np.max(sigind))
    omega = np.zeros(nomg)
    g_complex = np.zeros((nomg, norb, norb), dtype=complex)

    for n, line in enumerate(data_lines):
        vals = [float(x) for x in line.split()]
        omega[n] = vals[0]
        # Read (re, im) pairs for each sigind component
        sigind_vals = {}
        for k in range(max_idx):
            re_val = vals[1 + 2 * k]
            im_val = vals[2 + 2 * k]
            sigind_vals[k + 1] = complex(re_val, im_val)

        # Expand to full matrix
        for i in range(norb):
            for j in range(norb):
                idx = sigind[i, j]
                if idx > 0 and idx in sigind_vals:
                    g_complex[n, i, j] = sigind_vals[idx]

    return omega, g_complex


def preprocess_for_norg(params_norg_file, trans_dat_file=None,
                        delta_inp_file=None, sig_out_file=None,
                        Uc=None, CoulombF='Ising', spin_dep=False,
                        output_dir="."):
    """Convert eDMFT files to iNORG format (alternative to preprocess_edinput).

    This function works with separate eDMFT files instead of a combined input.tsv.

    Generates:
      - input_hb.txt            (hybridization in iNORG Green format)
      - input_seloc.txt         (previous self-energy, optional)
      - input_interaction.txt   (interaction tensor from eDMFT Uc)

    Args:
        params_norg_file: path to PARAMS.norg
        trans_dat_file: path to Trans.dat
        delta_inp_file: path to Delta.inp (eDMFT sigind format)
        sig_out_file: path to Sig.out (previous iteration, optional)
        Uc: 4D numpy array of Coulomb interaction (optional)
        CoulombF: interaction type ('Ising', 'Full', etc.)
        spin_dep: whether spin-dependent
        output_dir: directory to write output files
    """
    # Read PARAMS.norg
    params = read_params_norg(params_norg_file)
    nband = params.get('orbital', 4)
    norb = nband  # per spin
    norbs = 2 * norb
    beta = params.get('beta', 50.0)
    unit_omg = np.pi / beta

    # Read sigind
    if trans_dat_file and os.path.exists(trans_dat_file):
        sigind = read_sigind_matrix(trans_dat_file)
        # Truncate to norb x norb if needed
        if sigind.shape[0] > norb:
            sigind = sigind[:norb, :norb]
    else:
        sigind = build_sigind_for_test(norb)

    # Read Delta.inp and convert
    if delta_inp_file and os.path.exists(delta_inp_file):
        omega, delta_compact = read_edmft_green(delta_inp_file, sigind, norb)
        delta_full = expand_delta_from_edinput(delta_compact, norb)
        hb_path = os.path.join(output_dir, "input_hb.txt")
        write_inorg_green(hb_path, omega, delta_full, norbs)
        print(f"Written: {hb_path}")

    # Read Sig.out (previous iteration) and convert
    if sig_out_file and os.path.exists(sig_out_file):
        omega_sig, sig_compact = read_edmft_green(sig_out_file, sigind, norb)
        sig_full = expand_delta_from_edinput(sig_compact, norb)
        seloc_path = os.path.join(output_dir, "input_seloc.txt")
        write_inorg_green(seloc_path, omega_sig, sig_full, norbs)
        print(f"Written: {seloc_path}")

    # Write interaction from Uc tensor
    if Uc is not None:
        interaction_path = os.path.join(output_dir, "input_interaction.txt")
        write_interaction_from_Uc(interaction_path, Uc, sigind,
                                  CoulombF, spin_dep)
        print(f"Written: {interaction_path}")

    print("Preprocessing complete.")


def write_interaction_from_Uc(filename, Uc, Sigind, CoulombF='Ising',
                               spin_dep=False):
    """Convert eDMFT Uc tensor to iNORG interaction file.

    Following inorg.py:1226-1342 (prepare_EDinput) logic.

    Args:
        filename: output file path
        Uc: 4D numpy array [m1,m2,m3,m4] — Coulomb interaction tensor
        Sigind: sigind_matrix
        CoulombF: interaction type ('Ising', 'IsingS', 'IsingK', 'Full')
        spin_dep: whether spin-dependent
    """
    ntot = len(Sigind)
    divider = 2 if spin_dep else 1
    norb_half = ntot // divider
    # Count correlated orbitals
    norb = np.count_nonzero(np.diag(Sigind[:norb_half, :norb_half]))

    if CoulombF in ['Ising', 'IsingS', 'IsingK']:
        vij_abba = np.zeros((norb, norb))
        vij_abab = np.zeros((norb, norb))

        i1 = 0
        for m1 in range(norb_half):
            if Sigind[m1, m1] == 0:
                continue
            i2 = 0
            for m2 in range(norb_half):
                if Sigind[m2, m2] == 0:
                    continue
                for m3 in range(norb_half):
                    if Sigind[m3, m3] == 0:
                        continue
                    for m4 in range(norb_half):
                        if Sigind[m4, m4] == 0:
                            continue
                        if m1 == m4 and m2 == m3:
                            if abs(Uc[m1, m2, m3, m4]) > 1e-8:
                                vij_abba[i1, i2] = Uc[m1, m2, m3, m4].real
                        if m1 == m3 and m2 == m4 and m1 != m2:
                            if abs(Uc[m1, m2, m3, m4]) > 1e-8:
                                vij_abab[i1, i2] = Uc[m1, m2, m3, m4].real
                i2 += 1
            i1 += 1

        write_interaction_file(filename, vij_abba=vij_abba, vij_abab=vij_abab)
    else:
        # Full Coulomb
        vijkl = []
        i1 = 0
        for m1 in range(norb_half):
            if Sigind[m1, m1] == 0:
                continue
            i2 = 0
            for m2 in range(norb_half):
                if Sigind[m2, m2] == 0:
                    continue
                i3 = 0
                for m3 in range(norb_half):
                    if Sigind[m3, m3] == 0:
                        continue
                    i4 = 0
                    for m4 in range(norb_half):
                        if Sigind[m4, m4] == 0:
                            continue
                        if abs(Uc[m1, m2, m3, m4]) > 1e-8:
                            vijkl.append((i1, i2, i3, i4,
                                          Uc[m1, m2, m3, m4].real))
                        i4 += 1
                    i3 += 1
                i2 += 1
            i1 += 1

        write_interaction_file(filename, vijkl=vijkl)


def read_params_norg(filename):
    """Read PARAMS.norg file (key-value format with optional comments).

    Returns:
        dict of parameters
    """
    params = {}
    with open(filename, 'r') as f:
        for line in f:
            # Remove comments
            line = line.split('#')[0].strip()
            if not line:
                continue
            parts = line.split(None, 1)
            if len(parts) < 2:
                continue
            key = parts[0]
            val_str = parts[1].strip()

            # Try to parse as JSON-like value
            if val_str.startswith('['):
                # Parse array
                val_str = val_str.strip('[]')
                try:
                    vals = [float(x.strip()) for x in val_str.split(',') if x.strip()]
                    params[key] = vals
                except ValueError:
                    params[key] = val_str
            elif val_str.startswith("'") and val_str.endswith("'"):
                params[key] = val_str.strip("'")
            else:
                try:
                    if '.' in val_str or 'e' in val_str.lower():
                        params[key] = float(val_str)
                    else:
                        params[key] = int(val_str)
                except ValueError:
                    params[key] = val_str

    return params


# =============================================================================
# input.tsv -> solver_input.tsv conversion
# =============================================================================

def write_solver_input_tsv(filename, tij, omega, delta, vijkl):
    """Write solver_input.tsv in the format that C++ api_edmft.cpp expects.

    Format:
      Tij section: nband x nband real matrix (space-separated)
      Delta section: Python complex format per frequency line
        (omega+0j) (Delta[0,0]) (Delta[0,1]) ... (Delta[nband-1,nband-1])
      Vijkl section: sparse 4-index format (0-based)
        i j k l value

    Args:
        filename: output file path
        tij: (nband, nband) real array -- impurity energy levels
        omega: (nomg,) real array -- Matsubara frequencies
        delta: (nomg, nband, nband) complex array -- hybridization
        vijkl: list of (I, J, K, L, value) tuples or None
    """
    nband = tij.shape[0]

    with open(filename, 'w') as f:
        # --- Tij ---
        f.write("Tij\n")
        for i in range(nband):
            row = " ".join(f"{tij[i, j]:.18e}" for j in range(nband))
            f.write(row + "\n")
        f.write("\n")

        # --- Delta (Python complex format) ---
        f.write("Delta\n")
        for n in range(len(omega)):
            parts = [f" ({omega[n]:.18e}+{0.0:.18e}j)"]
            for i in range(nband):
                for j in range(nband):
                    val = delta[n, i, j]
                    re, im = val.real, val.imag
                    if im >= 0:
                        parts.append(f" ({re:.18e}+{im:.18e}j)")
                    else:
                        parts.append(f" ({re:.18e}{im:.18e}j)")
            f.write("".join(parts) + "\n")
        f.write("\n")

        # --- Vijkl (sparse 4-index, 0-based) ---
        if vijkl is not None:
            f.write("Vijkl\n")
            for (I, J, K, L, value) in vijkl:
                f.write(f"{I} {J} {K} {L} {value:.18e}\n")
        f.write("\n")


def convert_edinput_to_solver_input(edinput_file, output_file,
                                     sigind=None, transform_to_K=True,
                                     sym='psc', transform_kinetic=True):
    """Convert input.tsv to solver_input.tsv for C++ api_edmft.cpp.

    input.tsv (from prepare_EDinput in inorg.py) contains:
      - Tij: norb x norb (may include both spins if spin_dep)
      - Delta: norb^2 complex values per frequency (sigind-expanded)
      - Vij_abba/Vij_abab or Vijkl: interaction

    solver_input.tsv (for C++) needs:
      - Tij: nband x nband (one spin block)
      - Delta: nband x nband Python complex format
      - Vijkl: sparse 4-index format (0-based)

    The C++ internally handles:
      - Spin expansion: block_diag(up, -dn) for Nambu, block_diag(up, dn) for normal
      - Hartree shift computation from Vijkl (Nambu only)

    Args:
        edinput_file: path to input.tsv
        output_file: path to write solver_input.tsv
        sigind: sigind matrix (numpy array) for determining nband
        transform_to_K: whether to apply site->K basis transform
        sym: symmetry type ('psc'/'pfm'/'dsc'/'afm' for Nambu,
             'normal-pm'/'normal-af' for normal state)
        transform_kinetic: whether to apply the site->K rotation to Tij/Delta.
             Set False for the K-native path (dmft1 already emits a K-diagonal
             Delta.imp via the indmfl cftrans, so re-rotating is the redundant
             round-trip we eliminate). The interaction (Vijkl) is still
             transformed to K regardless, because the Coulomb matrix Uc is
             always built site-basis.
    """
    tij, omega, delta, vij_abba, vij_abab, vijkl = read_edinput(edinput_file)
    norb = tij.shape[0]

    is_nambu = not str(sym).startswith('normal')

    # nband = output matrix dimension, nband_phys = physical cluster sites
    # Nambu: nband=8 (full pp+Φ+Φ†+hh), nband_phys=4 (for unitary/interaction)
    # Normal: nband=4 (one spin block extracted), nband_phys=4
    if is_nambu and norb == 8:
        nband = norb
        nband_phys = norb // 2
    elif sigind is not None:
        ntot = sigind.shape[0]
        half = ntot // 2
        nband_up = np.count_nonzero(np.diag(sigind[:half, :half]))
        nband_dn = np.count_nonzero(np.diag(sigind[half:, half:]))
        if nband_up == nband_dn and norb == nband_up + nband_dn:
            nband = nband_up
            tij = tij[:nband, :nband]
            delta = delta[:, :nband, :nband]
        elif norb == nband_up:
            nband = norb
        else:
            nband = norb
        nband_phys = nband
    else:
        nband = norb
        nband_phys = norb

    print(f"convert_edinput_to_solver_input: norb={norb}, nband={nband}, nband_phys={nband_phys}, sym={sym}, nambu={is_nambu}")
    print(f"  omega range: [{omega[0]:.6f}, {omega[-1]:.6f}], {len(omega)} points")

    if transform_to_K and transform_kinetic and nband_phys == 4:
        U = plaquette_unitary()
        if is_nambu and nband == 8:
            import scipy.linalg
            U_nambu = scipy.linalg.block_diag(U, U)
            tij = np.real(transform_matrix_to_K(tij, U_nambu))
            delta = transform_green_to_K(delta, U_nambu)
            print(f"  Applied site->K basis transform (8x8 Nambu block-diag)")
        else:
            tij = np.real(transform_matrix_to_K(tij, U))
            delta = transform_green_to_K(delta, U)
            print(f"  Applied site->K basis transform (2x2 plaquette)")
    elif transform_to_K and not transform_kinetic and nband_phys == 4:
        print(f"  K-native: Tij/Delta already K-diagonal; skipping site->K rotation "
              f"(interaction Vijkl still transformed to K)")
    elif transform_to_K and nband_phys != 4:
        print(f"  WARNING: K-basis transform skipped (nband_phys={nband_phys} != 4)")

    vijkl_out = vijkl
    if vijkl is None and vij_abba is not None:
        if transform_to_K and nband_phys == 4:
            U = plaquette_unitary()
            vij_abab_use = vij_abab if vij_abab is not None else np.zeros_like(vij_abba)
            vijkl_out = transform_interaction_to_K_vijkl(
                vij_abba, vij_abab_use, U, nband_phys, nambu=is_nambu)
            print(f"  Interaction: Ising -> {len(vijkl_out)} Vijkl entries (K-basis, nambu={is_nambu})")
        else:
            vijkl_out = _ising_to_vijkl_site(vij_abba, vij_abab, nband_phys, nambu=is_nambu)
            print(f"  Interaction: Ising -> {len(vijkl_out)} Vijkl entries (site-basis, nambu={is_nambu})")

    write_solver_input_tsv(output_file, tij, omega, delta, vijkl_out)
    print(f"  Written: {output_file}")


def _ising_to_vijkl_site(vij_abba, vij_abab, nband, nambu=False):
    """Convert Ising interaction matrices to Vijkl sparse format (site-basis).

    For Nambu (nambu=True): sign = -1 (hole convention)
      abba: (a, b+nband, b+nband, a) with -V_abba
      abab: (a, a+nband, b+nband, b) with -V_abab (a!=b)

    For normal state (nambu=False): sign = +1
      abba: (a, b+nband, b+nband, a) with +V_abba
      abab: (a, a+nband, b+nband, b) with +V_abab (a!=b)

    Args:
        vij_abba: (nband_half, nband_half) density-density
        vij_abab: (nband_half, nband_half) pair-hopping (or None)
        nband: number of orbitals per spin
        nambu: if True, use Nambu sign convention (-V); if False, normal (+V)

    Returns:
        list of (i, j, k, l, value) tuples (0-based indices)
    """
    tol = 1e-12
    vijkl = []
    nband_half = vij_abba.shape[0]
    sign = -1.0 if nambu else +1.0

    for a in range(nband_half):
        for b in range(nband_half):
            # abba contribution
            if abs(vij_abba[a][b]) > tol:
                vijkl.append((a, b + nband, b + nband, a, sign * vij_abba[a][b]))

            # abab contribution (a != b)
            if vij_abab is not None and a != b and abs(vij_abab[a][b]) > tol:
                vijkl.append((a, a + nband, b + nband, b, sign * vij_abab[a][b]))

    return vijkl


# =============================================================================
# CLI interface
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Format conversion between eDMFT and iNORG solver")
    subparsers = parser.add_subparsers(dest='command', help='Sub-command')

    # Preprocess sub-command
    pre_parser = subparsers.add_parser('preprocess',
                                        help='Convert eDMFT files to iNORG format')
    pre_parser.add_argument('--edinput', type=str, default=None,
                            help='Path to input.tsv (combined format)')
    pre_parser.add_argument('--params', type=str, default='PARAMS.norg',
                            help='Path to PARAMS.norg')
    pre_parser.add_argument('--trans', type=str, default=None,
                            help='Path to Trans.dat')
    pre_parser.add_argument('--delta', type=str, default=None,
                            help='Path to Delta.inp (eDMFT sigind format)')
    pre_parser.add_argument('--sig', type=str, default=None,
                            help='Path to Sig.out (previous iteration)')
    pre_parser.add_argument('--output-dir', type=str, default='.',
                            help='Output directory')

    # Postprocess sub-command
    post_parser = subparsers.add_parser('postprocess',
                                         help='Convert iNORG output to eDMFT format')
    post_parser.add_argument('--sigind', type=str, default=None,
                             help='Path to Trans.dat')
    post_parser.add_argument('--nband', type=int, default=None,
                             help='Number of bands per spin')
    post_parser.add_argument('--seimp', type=str, default=None,
                             help='Path to self-energy file')
    post_parser.add_argument('--gfimp', type=str, default=None,
                             help='Path to Green function file')
    post_parser.add_argument('--output-dir', type=str, default='.',
                             help='Output directory')

    args = parser.parse_args()

    if args.command == 'preprocess':
        if args.edinput:
            preprocess_edinput(args.edinput, args.params, args.output_dir)
        elif args.delta:
            preprocess_for_norg(args.params, args.trans, args.delta,
                                args.sig, output_dir=args.output_dir)
        else:
            print("Error: specify --edinput or --delta for preprocessing")
            sys.exit(1)

    elif args.command == 'postprocess':
        sigind = None
        if args.sigind and os.path.exists(args.sigind):
            sigind = read_sigind_matrix(args.sigind)
        postprocess_from_norg(args.output_dir, sigind, args.nband,
                              args.seimp, args.gfimp)

    else:
        parser.print_help()


if __name__ == '__main__':
    main()
