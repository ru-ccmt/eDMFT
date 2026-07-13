#!/usr/bin/env python
"""
nambu_indmfl_gen.py — Generate Nambu-extended case.indmfl from normal-state indmfl.

For superconducting DMFT, the orbital space doubles to include particle-hole channels.
The Nambu sigind has 2N x 2N structure with 4 blocks:
    | pp (original)    | Phi_upper (+K)   |
    | Phi†_lower (+2K) | hh (+3K)         |
where K = max(original sigind) and zeros remain zeros.
"""
import sys
import os
import numpy as np


def _setup_import_path():
    """Add directory containing indmffile.py to sys.path.

    Search order: script's own directory, then <parent>/bin-cmt.rutgers/,
    then <parent>/bin/. The script's own directory is the canonical location
    in this project; the bin* fallbacks are kept for upstream Rutgers eDMFT
    layouts.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    candidates = [script_dir] + [os.path.join(base_dir, s) for s in ('bin-cmt.rutgers', 'bin')]
    for candidate in candidates:
        if os.path.isfile(os.path.join(candidate, 'indmffile.py')):
            if candidate not in sys.path:
                sys.path.insert(0, candidate)
            return candidate
    raise ImportError(
        "Cannot find indmffile.py in script directory, bin-cmt.rutgers/, or bin/"
    )


def expand_sigind_nambu(sigind):
    """Expand a normal-state sigind matrix into Nambu 2N x 2N form.

    Args:
        sigind: (N, N) int array with values 0..K

    Returns:
        (2N, 2N) int array with values 0..4K in four blocks:
          pp (1..K), Phi (K+1..2K), Phi_dag (2K+1..3K), hh (3K+1..4K)
    """
    K = int(np.max(np.abs(sigind.ravel())))
    ndim = sigind.shape[0]
    nambu = np.zeros((2 * ndim, 2 * ndim), dtype=int)
    mask = sigind != 0

    nambu[:ndim, :ndim] = sigind                      # pp block
    nambu[:ndim, ndim:][mask] = sigind[mask] + K      # Phi upper
    nambu[ndim:, :ndim][mask] = sigind[mask] + 2 * K  # Phi_dag lower
    nambu[ndim:, ndim:][mask] = sigind[mask] + 3 * K  # hh block

    return nambu


def expand_cftrans_nambu(cftrans):
    """Block-diagonalize crystal-field transformation for Nambu doubling."""
    ndim = cftrans.shape[0]
    nambu = np.zeros((2 * ndim, 2 * ndim), dtype=complex)
    nambu[:ndim, :ndim] = cftrans
    nambu[ndim:, ndim:] = cftrans
    return nambu


def expand_legends_nambu(legends):
    """Expand legend strings with Nambu block prefixes."""
    return (
        [f'Spp:{l}' for l in legends] +
        [f'Phi:{l}' for l in legends] +
        [f'PhiD:{l}' for l in legends] +
        [f'Shh:{l}' for l in legends]
    )


def generate_nambu_indmfl(normal_indmfl_path, nambu_indmfl_path, case):
    """
    Read a normal-state case.indmfl and generate a Nambu-extended version.

    The Nambu indmfl doubles the sigind dimension and creates 4 blocks:
      - pp (particle-particle): original sigind values 1..K
      - Phi_upper (anomalous):  sigind values K+1..2K
      - Phi_dag_lower:          sigind values 2K+1..3K
      - hh (hole-hole):         sigind values 3K+1..4K

    Where K = max(original sigind values), and zeros remain zeros.

    The cftrans (crystal-field transformation) is block-diagonalized:
      | cftrans_orig    0           |
      | 0               cftrans_orig |

    Args:
        normal_indmfl_path: path to normal-state case.indmfl
        nambu_indmfl_path: path for output Nambu case.indmfl
        case: Wien2k case name (e.g., "CCOC")
    """
    _setup_import_path()
    from indmffile import Indmfl

    inl = Indmfl(case)
    inl.read(normal_indmfl_path)

    for icp in list(inl.siginds.keys()):
        inl.siginds[icp] = expand_sigind_nambu(inl.siginds[icp])
        inl.cftrans[icp] = expand_cftrans_nambu(inl.cftrans[icp])
        inl.legends[icp] = expand_legends_nambu(inl.legends[icp])

    inl.write(nambu_indmfl_path)
    return inl


def detect_nambu_indmfl(indmfl_path, case):
    """Check if an indmfl file is already in Nambu format.

    Returns:
        (is_nambu, K_normal) where K_normal is the max sigind of the
        normal-state part.
    """
    _setup_import_path()
    from indmffile import Indmfl

    inl = Indmfl(case)
    inl.read(indmfl_path)

    for icp in inl.siginds:
        legends = inl.legends[icp]
        if any(l.startswith('Phi:') for l in legends):
            K_total = int(np.max(np.abs(inl.siginds[icp].ravel())))
            K_normal = K_total // 4
            return (True, K_normal)

    K_normal = max(int(np.max(np.abs(s.ravel()))) for s in inl.siginds.values())
    return (False, K_normal)


def _self_test():
    """Quick self-test of the expansion logic (no file I/O needed)."""
    print("Running self-test...")

    # Test 1: Basic 4x4 sigind with K=3
    sigind = np.array([
        [1, 2, 2, 3],
        [2, 1, 3, 2],
        [2, 3, 1, 2],
        [3, 2, 2, 1],
    ])

    nambu = expand_sigind_nambu(sigind)
    assert nambu.shape == (8, 8), f"Expected (8,8), got {nambu.shape}"

    # pp block (upper-left): unchanged
    assert np.array_equal(nambu[:4, :4], sigind), "pp block mismatch"

    # Phi upper block: +3
    expected_phi = np.array([
        [4, 5, 5, 6],
        [5, 4, 6, 5],
        [5, 6, 4, 5],
        [6, 5, 5, 4],
    ])
    assert np.array_equal(nambu[:4, 4:], expected_phi), "Phi upper block mismatch"

    # Phi_dag lower block: +6
    expected_phid = np.array([
        [7, 8, 8, 9],
        [8, 7, 9, 8],
        [8, 9, 7, 8],
        [9, 8, 8, 7],
    ])
    assert np.array_equal(nambu[4:, :4], expected_phid), "Phi_dag lower block mismatch"

    # hh block: +9
    expected_hh = np.array([
        [10, 11, 11, 12],
        [11, 10, 12, 11],
        [11, 12, 10, 11],
        [12, 11, 11, 10],
    ])
    assert np.array_equal(nambu[4:, 4:], expected_hh), "hh block mismatch"

    # Value range
    nonzero = nambu[nambu != 0]
    assert nonzero.min() == 1 and nonzero.max() == 12, \
        f"Expected range 1..12, got {nonzero.min()}..{nonzero.max()}"

    print("  sigind expansion: PASS")

    # Test 2: Zeros preserved
    sigind_sparse = np.array([
        [1, 0, 0, 2],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [2, 0, 0, 1],
    ])
    nambu_sparse = expand_sigind_nambu(sigind_sparse)
    assert np.all(nambu_sparse[1, :] == 0), "Zero row not preserved"
    assert np.all(nambu_sparse[:, 1] == 0), "Zero col not preserved"
    assert np.all(nambu_sparse[5, :] == 0), "Doubled zero row not preserved"
    assert np.all(nambu_sparse[:, 5] == 0), "Doubled zero col not preserved"

    print("  zero preservation: PASS")

    # Test 3: cftrans block-diagonal
    cft = np.eye(4, dtype=complex) * (0.5 + 0.3j)
    nambu_cft = expand_cftrans_nambu(cft)
    assert nambu_cft.shape == (8, 8), f"cftrans shape: expected (8,8), got {nambu_cft.shape}"
    assert np.array_equal(nambu_cft[:4, :4], cft), "cftrans upper-left"
    assert np.array_equal(nambu_cft[4:, 4:], cft), "cftrans lower-right"
    assert np.all(nambu_cft[:4, 4:] == 0), "cftrans upper-right should be zero"
    assert np.all(nambu_cft[4:, :4] == 0), "cftrans lower-left should be zero"

    print("  cftrans expansion: PASS")

    # Test 4: legends
    legends = ['x2y2-loc', 'x2y2-NN', 'x2y2-NNN']
    nambu_leg = expand_legends_nambu(legends)
    assert len(nambu_leg) == 12, f"Expected 12 legends, got {len(nambu_leg)}"
    assert nambu_leg[0] == 'Spp:x2y2-loc'
    assert nambu_leg[3] == 'Phi:x2y2-loc'
    assert nambu_leg[6] == 'PhiD:x2y2-loc'
    assert nambu_leg[9] == 'Shh:x2y2-loc'

    print("  legends expansion: PASS")
    print("All self-tests PASSED.")


if __name__ == '__main__':
    if len(sys.argv) == 4:
        generate_nambu_indmfl(sys.argv[1], sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 1:
        _self_test()
    else:
        print(f"Usage: {sys.argv[0]} <normal_indmfl> <nambu_indmfl> <case>")
        print(f"       {sys.argv[0]}               (run self-test)")
        sys.exit(1)
