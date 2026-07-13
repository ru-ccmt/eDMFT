#!/usr/bin/env python3
"""Expand Nambu sig.inp from K=3 (25 cols) to K=4 (33 cols).

Splits NN sigind into NN_x and NN_y for d-wave support.

Input (K=3, 12 Nambu components, 25 cols):
  sigind order: Σ_pp(loc,NN,NNN) | Φ(loc,NN,NNN) | Φ†(loc,NN,NNN) | Σ_hh(loc,NN,NNN)

Output (K=4, 16 Nambu components, 33 cols):
  sigind order: Σ_pp(loc,NN_x,NN_y,NNN) | Φ(...) | Φ†(...) | Σ_hh(...)

For each block: NN_x = NN, NN_y = NN (same value, since normal state has C4).
The key difference: they now have INDEPENDENT sigind numbers,
so the solver/sgather can evolve them separately during Nambu iterations.

Usage:
    python3 expand_siginp_K3_to_K4.py sig.inp [output_sig.inp]
"""
import numpy as np
import sys


def expand_K3_to_K4(infile, outfile):
    # Read header
    with open(infile) as f:
        hdr1 = f.readline()  # # s_oo= [...]
        hdr2 = f.readline()  # # Edc= [...]

    s_oo_old = eval(hdr1.split('=', 1)[1].strip())  # 12 elements (4 blocks × 3)
    Edc_old = eval(hdr2.split('=', 1)[1].strip())    # 12 elements

    assert len(s_oo_old) == 12, f"Expected 12 s_oo elements (K=3), got {len(s_oo_old)}"
    assert len(Edc_old) == 12, f"Expected 12 Edc elements (K=3), got {len(Edc_old)}"

    # Read data
    data = np.loadtxt(infile, comments='#')
    nfreq, ncol = data.shape
    assert ncol == 25, f"Expected 25 columns (8*3+1), got {ncol}"

    omega = data[:, 0]
    # 12 component pairs (Re, Im): 4 blocks × 3 components
    old_pairs = data[:, 1:]  # (nfreq, 24)

    print(f"Input: {nfreq} freqs, {ncol} cols, K=3 (12 Nambu components)")

    # Expand each block from 3 components to 4 components
    # Old order per block: [loc, NN, NNN] (Re,Im pairs)
    # New order per block: [loc, NN_x, NN_y, NNN]
    # NN_x = NN_y = NN (initial value, they diverge during iteration)

    new_pairs = np.zeros((nfreq, 32))  # 16 components × 2 (Re, Im)

    for block in range(4):  # pp, Phi, Phi†, hh
        old_base = block * 6   # 3 components × 2 (Re, Im) = 6 columns per block
        new_base = block * 8   # 4 components × 2 (Re, Im) = 8 columns per block

        # loc: old[0,1] → new[0,1]
        new_pairs[:, new_base + 0] = old_pairs[:, old_base + 0]  # Re(loc)
        new_pairs[:, new_base + 1] = old_pairs[:, old_base + 1]  # Im(loc)

        # NN → NN_x: old[2,3] → new[2,3]
        new_pairs[:, new_base + 2] = old_pairs[:, old_base + 2]  # Re(NN_x) = Re(NN)
        new_pairs[:, new_base + 3] = old_pairs[:, old_base + 3]  # Im(NN_x) = Im(NN)

        # NN → NN_y: old[2,3] → new[4,5] (same value, independent component)
        new_pairs[:, new_base + 4] = old_pairs[:, old_base + 2]  # Re(NN_y) = Re(NN)
        new_pairs[:, new_base + 5] = old_pairs[:, old_base + 3]  # Im(NN_y) = Im(NN)

        # NNN: old[4,5] → new[6,7]
        new_pairs[:, new_base + 6] = old_pairs[:, old_base + 4]  # Re(NNN)
        new_pairs[:, new_base + 7] = old_pairs[:, old_base + 5]  # Im(NNN)

    # Expand headers: 12 elements → 16 elements
    def expand_header(old_12):
        new_16 = []
        for block in range(4):
            base = block * 3
            loc = old_12[base + 0]
            nn  = old_12[base + 1]
            nnn = old_12[base + 2]
            new_16.extend([loc, nn, nn, nnn])  # loc, NN_x, NN_y, NNN
        return new_16

    s_oo_new = expand_header(s_oo_old)
    Edc_new = expand_header(Edc_old)

    # Write output
    with open(outfile, 'w') as f:
        f.write(f'# s_oo= {s_oo_new}\n')
        f.write(f'# Edc= {Edc_new}\n')
        for i in range(nfreq):
            parts = [f'{omega[i]:20.15f}']
            for j in range(32):
                parts.append(f'{new_pairs[i, j]:20.15f}')
            f.write(' '.join(parts) + '\n')

    print(f"Output: {nfreq} freqs, 33 cols, K=4 (16 Nambu components)")
    print(f"s_oo: {len(s_oo_new)} elements")
    print(f"Edc:  {len(Edc_new)} elements")
    print(f"Written to: {outfile}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python3 expand_siginp_K3_to_K4.py sig.inp [output]")
        sys.exit(1)
    infile = sys.argv[1]
    outfile = sys.argv[2] if len(sys.argv) > 2 else infile
    expand_K3_to_K4(infile, outfile)
