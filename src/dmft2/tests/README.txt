Independent regression tests for the 7 November 2026 stress audit fixes
=====================================================================
Requirements: Python with NumPy and SciPy (sph_harm_y), a Fortran compiler,
and BLAS/LAPACK. macOS uses Accelerate; Linux uses -llapack -lblas.
FC may select a compiler. No WIEN material data are required.

Run from any directory, using the paths below relative to this file:

  export EDMFT_STRESS_TEST_DIR=$(mktemp -d)
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  python build_stress_harness.py
  python check_stress_formulas.py
  python check_rank_deficient_gradient.py
  python check_projector_symmetry.py
  python check_complete_g2.py

The build scripts extract the actual production routines, without changing
source, and compile temporary copies with bounds checking. The full G2 test
uses the complete production module with a small radial/atomic fixture; only
private-cache visibility is changed in its test copy. All generated files and
results go to EDMFT_STRESS_TEST_DIR, never the source tree.

Tests and independent references:
* Fixed-index window response: exact small Hamiltonian endpoint diagonalization,
  generalized overlap, arbitrary complex KS phase changes, and functional FD.
* Global normalization: endpoint inverse square roots, Frechet/adjoint identity,
  and selected generalized-eigenvector span with its within-window metric.
* Normalization symmetry: exact adjoint of the full forward symmetry map,
  including noncommuting equivalent-block and orbital averages; endpoint
  inverse-square-root functional, default APPLY and diagnostic OFF/DIAG modes.
* Projector strain frame: global reciprocal vectors with fixed physical local
  axes in a nonorthogonal cell; reproduces the old equivalent-atom xz/xy error.
* Direct C: finite log-determinant derivatives on Matsubara and real-frequency
  contours; a rank-deficient projector whose physical gradient contains a
  component invisible to the embedded B=U*C matrix.
* Old and new angular kernels: independent SciPy spherical harmonics and
  numerical angular derivatives, including the largest-l intermediate and
  the strain-allowed s/s/Y20 channel absent from an ordinary .nsh record list.
* Complete G2: exact gamma-function radial integrals, u/udot and two genuine
  LO radial functions per l, independent angular quadrature, arbitrary complex
  pseudo-bands, a nonidentity local-axis rotation, and all six contractions.

These tests establish the implemented matrix algebra and local integrals.
They do not certify a complete material stress against charge-relaxed frozen-
Sigma endpoints, the frozen-radial approximation, all RELA contributions,
changing R_MT, or a fixed-index window cut through a degeneracy.
