#!/usr/bin/env python3
"""Convert normal-state solver_input.tsv to Nambu representation.

Usage:
    python3 convert_to_nambu.py input_normal.tsv output_nambu.tsv

Conversion rules (from doc/iNORG_Nambu_Interface.md):
  Tij:   T(4×4) → [T, 0; 0, -T] (8×8)
  Delta: Δ(4×4) → [Δ, 0; 0, -conj(Δ)] (8×8), anomalous blocks = 0
  Vijkl: V → -V (sign flip for Nambu hole convention)
"""
import sys
import numpy as np


def parse_complex(s):
    s = s.strip()
    if s.startswith('(') and s.endswith(')'):
        s = s[1:-1]
    return complex(s)


def read_solver_input(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Parse Tij
    idx = 0
    while idx < len(lines) and lines[idx].strip() != 'Tij':
        idx += 1
    idx += 1
    tij_rows = []
    while idx < len(lines) and lines[idx].strip():
        tij_rows.append([float(x) for x in lines[idx].split()])
        idx += 1
    tij = np.array(tij_rows)
    nband = tij.shape[0]

    # Parse Delta
    while idx < len(lines) and lines[idx].strip() != 'Delta':
        idx += 1
    idx += 1
    omega = []
    delta = []
    while idx < len(lines) and lines[idx].strip():
        parts = lines[idx].split()
        if not parts:
            idx += 1
            continue
        # First element is omega (real or complex with Im=0)
        om = parse_complex(parts[0]).real
        omega.append(om)
        # Remaining elements are nband*nband complex values
        vals = [parse_complex(p) for p in parts[1:]]
        delta.append(np.array(vals).reshape(nband, nband))
        idx += 1
    omega = np.array(omega)
    delta = np.array(delta)  # (nomg, nband, nband)

    # Parse Vijkl
    while idx < len(lines) and lines[idx].strip() != 'Vijkl':
        idx += 1
    idx += 1
    vijkl = []
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            idx += 1
            continue
        parts = line.split()
        if len(parts) >= 5:
            i, j, k, l = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])
            v = float(parts[4])
            vijkl.append((i, j, k, l, v))
        idx += 1

    return tij, omega, delta, vijkl, nband


def write_solver_input(filename, tij, omega, delta, vijkl):
    nband = tij.shape[0]
    with open(filename, 'w') as f:
        # Tij
        f.write("Tij\n")
        for i in range(nband):
            f.write(" ".join(f"{tij[i,j]:.18e}" for j in range(nband)) + "\n")
        f.write("\n")

        # Delta
        f.write("Delta\n")
        for n in range(len(omega)):
            parts = [f"({omega[n]:.18e}+0.000000000000000000e+00j)"]
            for i in range(nband):
                for j in range(nband):
                    z = delta[n, i, j]
                    parts.append(f"({z.real:.18e}{z.imag:+.18e}j)")
            f.write(" " + " ".join(parts) + "\n")
        f.write("\n")

        # Vijkl
        f.write("Vijkl\n")
        for (i, j, k, l, v) in vijkl:
            f.write(f"{i} {j} {k} {l} {v:.18e}\n")


def convert_to_nambu(tij_pp, omega, delta_pp, vijkl_normal):
    nband = tij_pp.shape[0]
    nband_nambu = 2 * nband

    # Tij: [T_pp, 0; 0, -T_pp]
    tij_nambu = np.zeros((nband_nambu, nband_nambu))
    tij_nambu[:nband, :nband] = tij_pp
    tij_nambu[nband:, nband:] = -tij_pp

    # Delta: [Δ_pp, 0; 0, -conj(Δ_pp)], anomalous = 0
    nomg = len(omega)
    delta_nambu = np.zeros((nomg, nband_nambu, nband_nambu), dtype=complex)
    delta_nambu[:, :nband, :nband] = delta_pp
    delta_nambu[:, nband:, nband:] = -np.conj(delta_pp)
    # Anomalous blocks (off-diagonal) = 0 (initial iteration)

    # Vijkl: sign flip V → -V
    vijkl_nambu = [(i, j, k, l, -v) for (i, j, k, l, v) in vijkl_normal]

    return tij_nambu, delta_nambu, vijkl_nambu


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python3 convert_to_nambu.py input_normal.tsv output_nambu.tsv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    tij, omega, delta, vijkl, nband = read_solver_input(input_file)
    print(f"Read normal-state: nband={nband}, nomg={len(omega)}, nvijkl={len(vijkl)}")
    print(f"  Tij diag: {np.diag(tij)}")
    print(f"  Vijkl[0] value: {vijkl[0][4]:.6f} (U/4 = {4*vijkl[0][4]:.3f})")

    tij_n, delta_n, vijkl_n = convert_to_nambu(tij, omega, delta, vijkl)
    print(f"\nNambu conversion:")
    print(f"  Tij_pp diag: {np.diag(tij_n)[:nband]}")
    print(f"  Tij_hh diag: {np.diag(tij_n)[nband:]}")
    print(f"  T_pp + T_hh = {np.diag(tij_n)[:nband] + np.diag(tij_n)[nband:]}")
    print(f"  Vijkl[0] value: {vijkl_n[0][4]:.6f} (= -U/4)")
    print(f"  Delta anomalous max: {np.max(np.abs(delta_n[:, :nband, nband:])):.2e}")

    write_solver_input(output_file, tij_n, omega, delta_n, vijkl_n)
    print(f"\nWritten: {output_file}")
