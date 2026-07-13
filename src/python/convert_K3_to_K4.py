#!/usr/bin/env python3
"""Convert normal-state files from K=3 to K=4 for d-wave Nambu CDMFT.

Splits NN sigind into NN_x (x-direction) and NN_y (y-direction) so that
d-wave pairing (cos kx - cos ky) can be expressed in the Nambu formalism.

Converts three files in-place:
  1. sig.inp:    7 cols (2K+1, K=3) → 9 cols (2K+1, K=4)
  2. lda.indmfl: sigind max=3 → max=4
  3. lda.indmfi: sigind max=3 → max=4

2×2 plaquette convention:
  sites: 0=(0,0), 1=(1,0), 2=(0,1), 3=(1,1)
  K=3: loc=1, NN=2, NNN=3
  K=4: loc=1, NN_x=2, NN_y=3, NNN=4

  NN_x pairs: (0,1),(1,0),(2,3),(3,2) — x direction
  NN_y pairs: (0,2),(2,0),(1,3),(3,1) — y direction

Usage:
    python3 convert_K3_to_K4.py [workdir]
    # workdir defaults to current directory
"""
import numpy as np
import sys
import os
import re


# NN_y positions in the 4-site cluster block (modulo cluster_size)
NN_Y_POSITIONS = {(0, 2), (2, 0), (1, 3), (3, 1)}
CLUSTER_SIZE = 4


def convert_siginp(filepath):
    """Convert sig.inp from K=3 (7 cols) to K=4 (9 cols)."""
    with open(filepath) as f:
        hdr1 = f.readline()
        hdr2 = f.readline()

    m1 = re.search(r's_oo\s*=\s*(\[.*?\])', hdr1)
    m2 = re.search(r'Edc\s*=\s*(\[.*?\])', hdr2)
    s_oo = eval(m1.group(1))
    Edc = eval(m2.group(1))

    if len(s_oo) != 3:
        print(f"  sig.inp: already K={len(s_oo)}, skipping")
        return False

    data = np.loadtxt(filepath, comments='#')
    omega = data[:, 0]
    # K=3: [Re(loc), Im(loc), Re(NN), Im(NN), Re(NNN), Im(NNN)]
    loc_cols = data[:, 1:3]
    nn_cols  = data[:, 3:5]
    nnn_cols = data[:, 5:7]

    # K=4: [loc, NN_x=NN, NN_y=NN, NNN]
    new_data = np.column_stack([omega.reshape(-1,1), loc_cols, nn_cols, nn_cols, nnn_cols])
    new_s_oo = [s_oo[0], s_oo[1], s_oo[1], s_oo[2]]
    new_Edc = [Edc[0], Edc[1], Edc[1], Edc[2]]

    with open(filepath, 'w') as f:
        f.write(f'# s_oo= {new_s_oo}\n')
        f.write(f'# Edc= {new_Edc}\n')
        for i in range(len(omega)):
            f.write(' '.join(f'{v:20.15f}' for v in new_data[i]) + '\n')

    print(f"  sig.inp: 7→9 cols (K=3→4)")
    return True


def convert_sigind_matrix(sigind, ndim):
    """Convert a sigind matrix from K=3 to K=4 (split NN→NN_x/NN_y).

    Two input conventions are supported:
      - K-native DIAGONAL Sigind [1,2,2,3] in K-order [Gamma,X,Y,M]:
        relabel each nonzero diagonal slot by its position 1,2,3,4 ->
        [1,2,3,4], lifting the X=Y degeneracy (matches
        init_dmft_kspace.build_plaquette_sigind(lift_xy=True)).
      - Site-basis CIRCULANT Sigind (off-diagonal NN bonds): split NN by
        NN_Y_POSITIONS (legacy iNORG path), NNN 3->4.
    """
    new_sigind = sigind.copy()

    # K-native diagonal case: lift X=Y by sequential position relabel.
    # Row order is K-order [Gamma,X,Y,M]; sequential counter over nonzero
    # diagonal slots gives [1,2,2,3] -> [1,2,3,4], byte-identical to
    # build_plaquette_sigind(lift_xy=True).
    if np.array_equal(sigind, np.diag(np.diag(sigind))):
        counter = 0
        for i in range(ndim):
            if sigind[i][i] != 0:
                counter += 1
                new_sigind[i][i] = counter
        return new_sigind

    # Site-basis circulant case (legacy iNORG): split NN -> NN_x / NN_y, NNN -> 4.
    nn_value = 2
    nnn_value = 3
    for i in range(ndim):
        for j in range(ndim):
            if sigind[i][j] == nn_value:
                ci, cj = i % CLUSTER_SIZE, j % CLUSTER_SIZE
                if (ci, cj) in NN_Y_POSITIONS:
                    new_sigind[i][j] = 3   # NN_y
                # else: NN_x stays as 2
            elif sigind[i][j] == nnn_value:
                new_sigind[i][j] = 4       # NNN shifts to 4

    return new_sigind


def convert_indmfl(filepath):
    """Convert lda.indmfl sigind from K=3 to K=4."""
    with open(filepath) as f:
        lines = f.readlines()

    # Find sigind matrix blocks and header lines, convert K=3→K=4
    modified = False
    new_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]

        # Update header lines: "N  dim  3" → "N  dim  4" (num-independent-components)
        # Also update the top-level header "1  20  3  # ..." → "1  20  4  # ..."
        line_nocomment = line.split('#')[0].strip()
        parts_nc = line_nocomment.split()
        comment = '#' + line.split('#', 1)[1] if '#' in line else ''
        if len(parts_nc) == 3 and all(p.isdigit() for p in parts_nc) and parts_nc[2] == '3':
            # Could be a header line like "1  20  3" (cix-num, dim, K)
            dim = int(parts_nc[1])
            if dim >= 4:  # reasonable dimension for sigind
                new_lines.append(f'{parts_nc[0]:>5s}{parts_nc[1]:>6s}   4       {comment}\n')
                i += 1
                modified = True
                continue

        # Update component labels: add NN_y
        if "'x2y2-NN'" in line and "'x2y2-NNN'" in line and "'x2y2-NN_y'" not in line:
            new_lines.append(line.replace("'x2y2-NNN'", "'x2y2-NN_y' 'x2y2-NNN'"))
            i += 1
            modified = True
            continue

        # Check if this line starts a sigind matrix block (all-integer, len>=4).
        # NOTE: the test must assemble the FULL block first, then test the whole
        # block for (max==3 and a 2 present).  A per-row test never fires on a
        # K-native DIAGONAL block ([1,2,2,3]: no single row has both max==3 and a 2),
        # so convert_sigind_matrix would never be called -> the K-native lift would
        # silently no-op.  This mirrors convert_indmfi's whole-block gate.
        parts = line_nocomment.split()
        if len(parts) >= 4 and all(p.lstrip('-').isdigit() for p in parts):
            vals = [int(p) for p in parts]
            ndim = len(vals)
            # Tentatively assemble an ndim x ndim block of all-integer rows.
            block_lines = [line]
            ok = True
            for j in range(1, ndim):
                if i + j >= len(lines):
                    ok = False
                    break
                bl_nc = lines[i + j].split('#')[0].strip().split()
                if len(bl_nc) < ndim or not all(p.lstrip('-').isdigit() for p in bl_nc[:ndim]):
                    ok = False
                    break
                block_lines.append(lines[i + j])

            if ok:
                block = np.array([[int(x) for x in bl.split()[:ndim]] for bl in block_lines])
                max_val = int(np.max(np.abs(block)))
                # Trigger on a K=3 block (max==3 and a 2 present) -- works for both
                # diagonal K-native [1,2,2,3] and site-basis circulant matrices.
                if max_val == 3 and (2 in block):
                    new_block = convert_sigind_matrix(block, ndim)
                    for r in range(ndim):
                        new_lines.append(' '.join(f'{int(new_block[r][c]):2d}' for c in range(ndim)) + '\n')
                    i += ndim
                    modified = True
                    continue

        new_lines.append(line)
        i += 1

    if modified:
        with open(filepath, 'w') as f:
            f.writelines(new_lines)
        print(f"  {os.path.basename(filepath)}: sigind K=3→4")
    else:
        print(f"  {os.path.basename(filepath)}: no K=3 sigind found, skipping")

    return modified


def convert_indmfi(filepath):
    """Convert lda.indmfi sigind from K=3 to K=4."""
    with open(filepath) as f:
        lines = f.readlines()

    # indmfi format:
    # line 0: number of sigind blocks
    # for each block: dimension line, then dim×dim matrix
    new_lines = [lines[0]]
    i = 1
    modified = False

    while i < len(lines):
        dim_line = lines[i]
        parts = dim_line.split()
        if not parts or not parts[0].isdigit():
            new_lines.append(dim_line)
            i += 1
            continue

        ndim = int(parts[0])
        new_lines.append(dim_line)
        i += 1

        # Read ndim rows of sigind
        block = []
        for r in range(ndim):
            if i + r < len(lines):
                row = [int(x) for x in lines[i + r].split()[:ndim]]
                block.append(row)

        if len(block) == ndim:
            block = np.array(block)
            max_val = int(np.max(np.abs(block)))
            if max_val == 3 and 2 in block:
                new_block = convert_sigind_matrix(block, ndim)
                for r in range(ndim):
                    new_lines.append(' '.join(f'{int(new_block[r][c]):2d}' for c in range(ndim)) + '\n')
                modified = True
            else:
                for r in range(ndim):
                    new_lines.append(lines[i + r])
        i += ndim

    if modified:
        with open(filepath, 'w') as f:
            f.writelines(new_lines)
        print(f"  {os.path.basename(filepath)}: sigind K=3→4")
    else:
        print(f"  {os.path.basename(filepath)}: no K=3 sigind found, skipping")

    return modified


def main():
    workdir = sys.argv[1] if len(sys.argv) > 1 else '.'
    print(f"Converting K=3 → K=4 in: {workdir}")

    sig_path = os.path.join(workdir, 'sig.inp')
    indmfl_path = os.path.join(workdir, 'lda.indmfl')
    indmfi_path = os.path.join(workdir, 'lda.indmfi')

    converted = 0
    for path, converter in [
        (sig_path, convert_siginp),
        (indmfl_path, convert_indmfl),
        (indmfi_path, convert_indmfi),
    ]:
        if os.path.exists(path):
            if converter(path):
                converted += 1
        else:
            print(f"  {os.path.basename(path)}: not found, skipping")

    print(f"\nDone. {converted} file(s) converted.")
    if converted > 0:
        print("Verify: head -2 sig.inp  →  s_oo should have 4 elements")
        print("        grep max indmfl  →  sigind max should be 4")


if __name__ == '__main__':
    main()
