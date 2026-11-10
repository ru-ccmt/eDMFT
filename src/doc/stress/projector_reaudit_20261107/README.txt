Genuinely eDMFT projector stress audit, 2026-11-07

Read edmft_projector_reaudit.pdf (editable source: edmft_projector_reaudit.tex
plus material_results.tex). No production source or original material input
was changed by this audit. Material runs used /private/tmp/projector_reaudit_20261107.

Main reproducible records:
  joint_projector_results.json          Combined functional derivative, static/dynamic
  symmetry_preserving_sigma.json       Same failure with symmetry-compatible Sigma
  potential_stationarity_results.json  Finite-window versus full-space stationarity
  material_summary.json                Material diagnostics and source/binary hashes
  harness/results.json                 Individual production-Fortran checks
  harness/rank_deficient_C.json         Rank-deficient direct gradient check
  harness/source_hashes.json            Sources used by the extracted harness

New tests call production Fortran routines extracted and compiled with bounds
checking. Their finite-difference endpoint oracle independently diagonalizes
both the generalized Hamiltonian and the inverse-square-root normalization.
The dynamic oracle uses a finite auxiliary bath, so its free energy is exact
without a Matsubara-tail approximation.

Reproduction of small algebraic tests (use a fresh scratch directory):
  export EDMFT_STRESS_TEST_DIR=/private/tmp/projector_reaudit_reproduce/harness
  python /Users/haule/projects/eDMFT_new/src/dmft2/tests/build_stress_harness.py
  python /Users/haule/projects/eDMFT_new/src/dmft2/tests/check_stress_formulas.py
  python /Users/haule/projects/eDMFT_new/src/dmft2/tests/check_rank_deficient_gradient.py
Copy joint_projector_test.py and potential_stationarity_test.py into
/private/tmp/projector_reaudit_reproduce, then execute them there. They use
its harness/audit_routines.dylib and write results into that scratch directory.
The production test builder resolves sources relative to its repository path;
do not move it and expect its source-root detection to remain correct.

Material protocols:
  material_diagnostics.py copies the two central inputs, generates vectors
  with x_dmft.py -p lapw1, and runs the unchanged installed dmft2 with
  EDMFT_STRESS_NORM_SYMMETRY=DIAG and APPLY. Diagnostic marker files enable
  first-k raw-projector finite differences and normalization dumps.
  fe_atom_shear_probe.py restricts the scratch input to one correlated Fe
  atom and one k point to isolate geometry only. These are NOT valid
  pressure calculations. With MPI size one dftKS writes serial vector names,
  so these probes prepare serial dmft2 vector definitions.
  summarize.py assembles small results and regenerates material_results.tex.

Only small diagnostic outputs are archived here, not vector files or the
entire material directories. Input hashes and the original calculation paths
are in the provenance records. Installed binary hashes are recorded. The
harness directory also retains its extracted Fortran as a snapshot; its
compiled library is intentionally omitted from the permanent report.

Limitations: no proof of the material size of the induced-potential chain;
no completed radial-response repair; no MnO/MnTe or new force benchmark in
this pass. The symmetry-chain correction is not a solution to the pressure
mismatch. FeSe shear failures are localized to the symmetry-related Fe atom,
but physical derivative versus restricted finite-difference-helper conventions
still need a consistent broken-symmetry test.
