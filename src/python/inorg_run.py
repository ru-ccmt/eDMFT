#!/usr/bin/env python3
"""
inorg_run.py -- Bridge script to run iNORG solver from eDMFT.

coded by Jia-Ming Wang (jiaming.w@rutgers.edu, Rutgers University, U.S.) date 2025 - now

Workflow:
  1. Read mpi_prefix.dat for MPI command
  2. Read Trans.dat for sigind matrix
  3. Convert input.tsv -> solver_input.tsv
     (site-basis sigind-compressed -> K-basis nband x nband, Python complex format)
  4. Run: {mpi_prefix} {inorg_exe}
  5. Postprocess: io/seimp.txt, io/gfimp.txt -> Sig.out, Gf.out
"""

import sys
import os
import subprocess
import numpy as np

# Add this directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from inorg_edmft_io import (
    read_sigind_matrix,
    read_params_norg,
    convert_edinput_to_solver_input,
    postprocess_from_norg,
)


def main():
    print("=" * 60)
    print("inorg_run.py: iNORG solver bridge for eDMFT")
    print("=" * 60)

    # 1. Read MPI prefix
    mpi_prefix = ""
    if os.path.exists("mpi_prefix.dat"):
        with open("mpi_prefix.dat", 'r') as f:
            mpi_prefix = f.read().strip()
        print(f"MPI prefix: {mpi_prefix}")
    else:
        print("Warning: mpi_prefix.dat not found, running without MPI")

    # 2. Read Trans.dat for sigind matrix
    sigind = None
    nband = None
    if os.path.exists("Trans.dat"):
        sigind = read_sigind_matrix("Trans.dat")
        # Determine nband from sigind diagonal
        ntot = sigind.shape[0]
        half = ntot // 2
        nband = np.count_nonzero(np.diag(sigind[:half, :half]))
        print(f"Read Trans.dat: ntot={ntot}, nband={nband}")
        print(f"  sigind shape: {sigind.shape}")
    else:
        print("Warning: Trans.dat not found, using defaults")

    # 3. Check that PARAMS.norg exists and read sym parameter
    if not os.path.exists("PARAMS.norg"):
        print("ERROR: PARAMS.norg not found")
        sys.exit(1)
    params_norg = read_params_norg("PARAMS.norg")
    sym = params_norg.get('sym', 'psc')
    print(f"PARAMS.norg found, sym={sym}")

    # K-native discrimination via PP-BLOCK Sigind diagonality (works for normal
    # AND Nambu).  When the particle-particle (normal-state) Sigind block is
    # DIAGONAL (plaquette K-basis [1,2,2,3] or lifted [1,2,3,4]), dmft1/dmftgk
    # already delivered a K-diagonal Delta.imp (the cftrans does real->K inside
    # dmft1) and prepare_EDinput builds a K-diagonal input.tsv.  The redundant
    # Python site<->K rotation of the kinetic/Delta block (and the inverse K->site
    # of Sig.out/Gf.out/Phi.out/F.out) must then be SKIPPED.  The interaction
    # (Vijkl) transform is still applied (Uc is built site-basis).  A site-circulant
    # Sigind is the legacy hybrid path -> keep the transform.  (sym stays
    # 'normal-pm'/'K-dsc'; the C++ binary forms the K-token from basis=K + sym.)
    # Nambu sigind is 2N x 2N with the anomalous Phi OFF the main diagonal; the pp
    # block sig[:N,:N] carries the K-diagonal normal structure -> test that block.
    is_k_native = False
    if sigind is not None and sigind.size > 0:
        _half = sigind.shape[0] // 2
        _sig_pp = sigind[:_half, :_half] if _half >= 1 else sigind
        _offdiag = _sig_pp - np.diag(np.diag(_sig_pp))
        is_k_native = (np.count_nonzero(_offdiag) == 0)
    if is_k_native:
        print("  K-native mode (diagonal pp-block Sigind): skipping Python site<->K "
              "rotation of Tij/Delta and the inverse K->site of Sig.out/Gf.out")

    # 4. Convert input.tsv -> solver_input.tsv
    if os.path.exists("input.tsv"):
        print("\nConverting input.tsv -> solver_input.tsv ...")
        convert_edinput_to_solver_input(
            "input.tsv", "solver_input.tsv",
            sigind=sigind, transform_to_K=True, sym=sym,
            transform_kinetic=(not is_k_native))
    else:
        print("ERROR: input.tsv not found")
        sys.exit(1)

    # 5. Find and run C++ solver
    # Pick the variant from params.dat:
    #   nambu = True  -> inorg-nambu (superconducting)
    #   otherwise     -> inorg-normal (normal-state)
    # Fall back to plain "inorg" if the specific variant is missing.
    wien_root = os.environ.get("WIEN_DMFT_ROOT",
                                os.path.dirname(os.path.abspath(__file__)))
    # The launcher (inorg_nambu.py / inorg.py) runs us with cwd=imp.k/,
    # but params.dat lives in the parent (DMFT run) directory. Check both.
    is_nambu = False
    params_path = None
    for candidate in ("params.dat", "../params.dat"):
        if os.path.exists(candidate):
            params_path = candidate
            break
    if params_path is not None:
        ns = {}
        try:
            exec(compile(open(params_path).read(), params_path, "exec"), ns)
            is_nambu = bool(ns.get("nambu", False))
            print(f"Read {params_path}: nambu={is_nambu}")
        except Exception as e:
            print(f"Warning: failed to parse {params_path} ({e}); assuming normal state")
    else:
        print("Warning: params.dat not found in . or ..; assuming normal state")
    preferred = "inorg-nambu" if is_nambu else "inorg-normal"
    print(f"selecting solver variant: {preferred}")

    inorg_exe = None
    for candidate in (
        os.path.join(wien_root, preferred),
        os.path.join(".", preferred),
        os.path.join(wien_root, "inorg"),
        "./inorg",
    ):
        if os.path.exists(candidate):
            inorg_exe = candidate
            break
    if inorg_exe is None:
        print(f"ERROR: no inorg binary found. Tried: {preferred}, inorg "
              f"(in {wien_root} and .)")
        sys.exit(1)
    print(f"Using solver binary: {inorg_exe}")

    cmd = f"{mpi_prefix} {inorg_exe}".strip()
    print(f"\nRunning C++ solver: {cmd}")
    print("-" * 60)
    ret = subprocess.call(cmd, shell=True)
    print("-" * 60)

    if ret != 0:
        print(f"ERROR: iNORG solver returned exit code {ret}")
        print("Check nohup_imp.out for details.")
        sys.exit(ret)

    print("\nC++ solver completed successfully")

    # 6. Postprocess: io/seimp.txt, io/gfimp.txt -> Sig.out, Gf.out
    print("\nPostprocessing output ...")
    postprocess_from_norg(".", sigind=sigind, nband=nband, sym=sym,
                          transform_to_site=(not is_k_native))

    print("\n inorg_run.py completed successfully")


if __name__ == "__main__":
    main()
