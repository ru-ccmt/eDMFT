MPI strain response: implementation and usage (9 November 2026)
=============================================================

The analytic response is distributed over k points. Each MPI rank reads and
prepares only its owned eigenvectors, H/O matrices, complement, projector and
density maps. Global k weights are never renormalized within a pool. Shared
radial/potential/core data and the real-space Krylov vectors are replicated.
The existing nonmagnetic, fixed-N, projector-5, primitive orthogonal scope and
symmetry-enforced diagonal tensor validation remain in force.

Execution
---------
The existing run_dmft.py final-measurement policy and x_dmft.py dmft2 immediate
policy are preserved. STRT/STRFT still request the complete response. Add these
solver controls to case.stress_response.json (merge with any existing settings):

{
  "solver": {
    "mpi": ["/Users/haule/miniconda3-clean/envs/work311/bin/mpiexec", "-n", "4"],
    "native_mpi": ["/opt/homebrew/bin/mpirun", "-np", "4"],
    "cache_mb": 256,
    "basis_cache_mb": 128
  }
}

These are the verified launchers on the laptop: Python mpi4py uses MPICH,
whereas the native dmft2 worker uses Open MPI. On rupc select launchers matching
its libraries after source ~/wmodules. The lists contain executable arguments,
not shell expressions. solver.mpi launches both the Python response and the
Python wrapper managing independent serial WIEN LAPW1 k pools. native_mpi is
only for the two native dmft2 passes. Ordinary production dmft2 retains its
existing MPI launch. Empty lists (the defaults) retain serial execution.

No material input file was changed during this installation. Activating MPI
for a material is an explicit configuration choice. At least one k point per
response rank is required; several k points per rank are preferable for scaling.
The two validation meshes contain 75 and 140 irreducible points, respectively.
They cannot demonstrate scaling to 300 ranks. A denser cluster benchmark is
still required before quoting a 100-300-core efficiency.

The two cache budgets are per rank, in MiB, for electronic frequency/contraction
caches and central spherical products. They are ADDITIONAL to the distributed
basis storage, not total process memory limits. Setting them to zero disables
these extra caches. Default caches do not retain all band wavefunctions on FFT
grids. FFT dimensions are chosen with a global maximum so every rank uses the
same alias-free Fourier convention. One BLAS/OpenMP thread per MPI rank is used;
newer macOS Accelerate additionally uses BLASSetThreading(single), since the
environment variable alone did not reliably constrain its newer interface.

Preparation
-----------
The native central pass remains isolated in case.response_runs/ID/central:
  lapw0_response -> LAPW1 k pools -> dmft2_number -> dmft2_export -> lcore_response
There is no strained SCF loop. The two dmft2 passes first establish mu and then
export that same state at fixed mu; the analytic derivative subsequently imposes
fixed N globally. No physics was changed to obtain the parallelism.

LAPW1 pools use the unchanged WIEN worker. Each has private files under
central/lapw1_pools/rankNNNNN. Original k labels, integration weights, and first-
line energy limits are preserved. The merger verifies equal radial headers and
identical NSH files, keeps one radial header, and concatenates the k records in
the original order. Matrix values are unchanged; their integer global k/rank
identifiers are remapped. Root's SCF1 supplies the common radial searches once.
Individual worker output logs remain in their own directories.

This removes the serial native eigensolve bottleneck without changing the
Hamiltonian quadrature or introducing a distributed dense eigensolver. Ordinary
production LAPW1 exports are not yet reused: the response continues to build a
fresh, consistent isolated snapshot. Such reuse is an additional future speed
optimization, not required for distributing this native preparation pass.
LAPW0, LCORE, small shared radial setup, stream merging and the common real-space
potential operations remain outside the k decomposition.

Response and exact adjoints
---------------------------
Write J=N[REB+D+C]K, and (I-J)^T lambda=c. Here K is the density-to-potential map,
B the native potential-to-H/O/raw-projector map, E the frozen-Sigma electronic
response, R the band-to-spatial-density reconstruction, D its direct radial
basis contribution, C core response, and N physical charge normalization.

MPI collectives have these meanings:
  B: shared potential -> local k matrices; B^T SUMS local potential covectors.
  R: local band densities -> SUMMED real-space density; R^T distributes the
     same real-space covector to each rank.
  D: both input and output shared; SUM local contributions in BOTH directions.
  C, K, N: shared operations counted ONCE, not summed over MPI ranks.

The electronic map E contains its own global reductions:
  S = Sym sum_k w_k P_k^dagger P_k, X=S^(-1/2).
  dS is globally summed before the Frechet derivative dX. In reverse, the
  covector of dX is globally summed before normalization/symmetry pullbacks.
  Fixed-N dmu uses global dN and global dN/dmu; its transpose uses the global
  mu-response contraction. The actual exported real-linear symmetry transpose
  is retained. The minimum window gap is also a global reduction.

Krylov vectors live in shared real-space coordinates. Their dot products are
ordinary Euclidean products, NOT an MPI sum of identical copies. Distributed
band-coordinate dot tests do sum local pairings. response_parallel.py exposes
both types explicitly; ChargeResponseChain.check_transposes uses the correct
input/output distribution for every leaf.

Shared source contractions
--------------------------
With one converged adjoint, a full electronic response per strain is unnecessary.
For fixed radial energies,
  s_nu = N[RE b_nu + (REB+D+C)v_nu + g_nu] + h_nu,
where b_nu is explicit H/O/raw-projector geometry, v_nu the explicit potential,
g_nu direct spatial geometry, and h_nu normalization geometry plus inverse
material-coordinate transport. Prepare ONCE
  eta=N^T lambda,
  b_bar=E^T R^T eta,
  v_bar=B^T b_bar + D^T eta + C^T eta.
Then
  lambda.s_nu = <b_bar,b_nu> + <v_bar,v_nu> + <eta,g_nu> + <lambda,h_nu>.
K-local contractions are summed; shared terms are counted once. The already
validated searched-energy term lambda.P A v_nu is added with its matching direct
functional contribution. This is an exact transposition of the original source,
not a new stress approximation. No second central-cell normalization is applied
after material-coordinate transport.

The full-source diagnostic remains available:
  lapw_total_stress_driver.py ... --full-sources
Default source*.npz archives contain the induced scalar contraction; their
source array is empty because the full vector was not formed. The report says
sources_contracted_by_adjoint=true. --full-sources restores the explicit vectors.
Use --save-geometry to retain large band-geometry arrays; under MPI these are
rank-local geometryN.rankRRRRR.npz files with original global k IDs. Default
production output avoids gathering or storing all those arrays on root.

Central k coordinates, site frames, local-orbital indexing, VNS data, and the
central H/O reconstruction check are reused across strain directions. Native
matrix files are not reread for each direction. Fingerprints are streamed and
distributed, then checked again before a result is accepted.

Validation and records
----------------------
See RESULTS.json / REPORT.txt for measurements. Native LAPW1 pool vectors,
energies and NSH are byte-identical for SrVO3 NREL (complex vectors) and FeSe RELA
(real vectors); H/O matrix values are exactly equal. Full SrVO3 native MPI
preparation also gives byte-identical response metadata and k exports compared
with serial LAPW1 plus MPI dmft2. Earlier serial versus MPI dmft2 differs only at
roundoff in mu and normalization.

The complex synthetic MPI test deliberately uses a non-self-adjoint real-linear
orbital symmetry map, unequal k weights, moving metric/projectors, fixed-N and
fixed-mu ensembles, and both cached and uncached frequency paths. Run it with
1, 2 or 4 ranks using tests/mpi_response_check.py.

Installation copies only the changed helpers/tests into the repository and
helpers into ~/dbin. Existing native executables and live material inputs are
unchanged. INSTALL_MANIFEST.json and before_installation preserve identities and
previous files; mpi_response.patch records source differences.
