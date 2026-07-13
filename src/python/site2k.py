#!/usr/bin/env python3
"""
site2k.py — Transform solver_input.tsv from site-basis to K-basis.

Reads solver_input.tsv (Tij, Delta, Vij_abba, Vij_abab in site-basis)
and writes solver_input_k.tsv with all quantities transformed to K-basis
using the 2x2 plaquette unitary:

    U[K,r] = (1/2) * exp(i K.r)

    Site positions: r0=(0,0), r1=(1,0), r2=(0,1), r3=(1,1)
    K-points:       K0=(0,0), K1=(pi,0), K2=(0,pi), K3=(pi,pi)

Transforms:
    T^K       = U @ T^site @ U^dagger
    Delta^K   = U @ Delta^site @ U^dagger   (per frequency)
    V^K_{IJKL} = sum_{ijkl} U*[I,i] U*[J,j] U[K,k] U[L,l] V^site_{ijkl}

Usage:
    python site2k.py [input_file] [output_file]
"""

import sys
import os

# Add parent directory so we can import from bin-cmt.rutgers
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin-cmt.rutgers'))

import numpy as np
from inorg_edmft_io import (
    read_edinput,
    plaquette_unitary,
    transform_matrix_to_K,
    transform_green_to_K,
    transform_interaction_to_K_vijkl,
)


def write_solver_input_tsv(filename, tij, omega, delta, vijkl):
    """Write solver_input_k.tsv format.

    Args:
        filename: output file path
        tij: (norb, norb) real array
        omega: (nomg,) real array — Matsubara frequencies
        delta: (nomg, norb, norb) complex array
        vijkl: list of (I, J, K, L, value) tuples — K-basis interaction
    """
    norb = tij.shape[0]

    with open(filename, 'w') as f:
        # --- Tij ---
        f.write("Tij\n")
        for i in range(norb):
            row = " ".join(f"{tij[i, j]:.18e}" for j in range(norb))
            f.write(row + "\n")

        # blank line
        f.write("\n")

        # --- Delta ---
        f.write("Delta\n")
        for n in range(len(omega)):
            parts = [f" ({omega[n]:.18e}+{0.0:.18e}j)"]
            for i in range(norb):
                for j in range(norb):
                    val = delta[n, i, j]
                    re, im = val.real, val.imag
                    if im >= 0:
                        parts.append(f" ({re:.18e}+{im:.18e}j)")
                    else:
                        parts.append(f" ({re:.18e}{im:.18e}j)")
            f.write("".join(parts) + "\n")

        # blank line
        f.write("\n")

        # --- Vijkl ---
        if vijkl is not None:
            f.write("Vijkl\n")
            for (I, J, K, L, value) in vijkl:
                f.write(f"{I} {J} {K} {L} {value:.18e}\n")

        f.write("\n")


def main():
    # Default file paths
    input_file = os.path.join(os.path.dirname(__file__), "solver_input.tsv")
    output_file = os.path.join(os.path.dirname(__file__), "solver_input_k.tsv")

    if len(sys.argv) >= 2:
        input_file = sys.argv[1]
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]

    print(f"Reading: {input_file}")
    tij, omega, delta, vij_abba, vij_abab, vijkl = read_edinput(input_file)

    norb = tij.shape[0]
    nomg = len(omega)
    print(f"  norb = {norb}")
    print(f"  nomg = {nomg}")
    print(f"  omega range = [{omega[0]:.6f}, {omega[-1]:.6f}]")

    # Build 2x2 plaquette unitary
    U = plaquette_unitary()
    print(f"  Plaquette unitary (4x4):")
    print(U)

    # Transform Tij to K-basis and zero out numerical noise
    tij_K = transform_matrix_to_K(tij, U)
    tij_K[np.abs(tij_K) < 1e-12] = 0.0
    print(f"\n  Tij (site-basis) diagonal = {np.diag(tij)}")
    print(f"  Tij (K-basis) diagonal    = {np.diag(tij_K)}")

    # Transform Delta to K-basis
    delta_K = transform_green_to_K(delta, U)
    print(f"  Delta transformed to K-basis: shape {delta_K.shape}")

    # Transform interaction to K-basis via 4-index tensor transform
    vijkl_K = None
    if vij_abba is not None:
        vij_abab_use = vij_abab if vij_abab is not None else np.zeros_like(vij_abba)
        vijkl_K = transform_interaction_to_K_vijkl(vij_abba, vij_abab_use, U, norb)
        print(f"  Vij_abba (site) diagonal = {np.diag(vij_abba)}")
        if vij_abab is not None:
            print(f"  Vij_abab (site) diagonal = {np.diag(vij_abab)}")
        print(f"  Vijkl (K-basis): {len(vijkl_K)} non-zero entries")

    # Write output
    write_solver_input_tsv(output_file, tij_K, omega, delta_K, vijkl_K)
    print(f"\nWritten: {output_file}")


if __name__ == '__main__':
    main()
