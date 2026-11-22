# Experimental dftKS eigensolvers

Development started 20 November 2026. Both real and complex executables support
the following environment switch; unset means the unchanged legacy solver.

| `DFTKS_EIGENSOLVER` | Behavior |
| --- | --- |
| `GVX` | Original `DSYGVX` / `ZHEGVX`; default and reference. |
| `GVD` | `DSYGVD` / `ZHEGVD`, divide and conquer. |
| `DAVIDSON` | Cached generalized block Davidson, with dense seeding/fallback. |
| `AUTO` | Davidson when a compatible cached subspace occupies at most 25% of the basis; otherwise GVD. This initial threshold is a heuristic, not a measured universal optimum. |

For example, in an **isolated test calculation**, using a newly built executable:

```sh
DFTKS_EIGENSOLVER=GVD /path/to/new/lapw1c lapw1.def
DFTKS_EIGENSOLVER=DAVIDSON /path/to/new/lapw1c lapw1.def
```

The first Davidson call seeds its cache using GVD. Later invocations reuse the
orbitals at each k point. No changes to `case.in1`, `case.in1c`, `params.dat`,
or the Python driver are required. The switch takes effect only in a rebuilt
dftKS executable; it does not select a native WIEN2k solver.

## Band and coefficient conventions

All modes preserve the existing selections:

* With `FIXED_NBANDS`, return the lowest requested number of eigenpairs.
* Otherwise return eigenpairs in `(Elow, Eup]`, like the original GVX call.
  The numerical fourth field in the K-VECTORS record alone does not enable
  fixed-band selection.

GVD computes the full eigensystem and then applies that selection. Its larger
workspace can make it unsuitable on memory-constrained large cells, even if
it is faster. LAPACK workspace queries and allocation results are checked.

The backend explicitly reconstructs the physical Hermitian H and O from the
existing two-triangle LAPW storage **before** the legacy conjugation step.
The vector/energy/NSH formats and TAPEWF conjugation convention are unchanged.
The matrix-only response path continues to bypass diagonalization entirely.

## Generalized Davidson algorithm

1. Validate the cached basis dimension, local-orbital count, full integer basis
   labels, and fractional k point. Discard an absent, incompatible or truncated
   cache and use a dense seed.
2. Construct and factor the current `H - sigma*O` with pivoted Hermitian LDL
   (`SYTRF` / `HETRF`). Its inertia counts the eigenvalues below sigma, assuming
   the physical LAPW overlap is positive definite.
3. Use the same factors to precondition correction blocks (`SYTRS` / `HETRS`).
   No explicit matrix inverse is formed. An estimated reciprocal condition
   number below `1e-12`, a singular pivot or ambiguous 2x2 inertia block causes
   dense fallback.
4. Orthonormalize cached vectors in the **current** overlap metric. Project
   corrections against the trial space twice. Diagonalize their small Gram
   matrix and discard dependent combinations using relative threshold `1e-10`.
   Do not regularize, clip or remove directions from the physical LAPW overlap.
5. Expand a limited block of unconverged residuals. Diagonalize the projected
   Hamiltonian, retain guard bands and restart at the configured subspace size.
6. Require converged generalized residuals for **all retained bands including
   guards**, metric orthonormality, and a band-count match with the inertia.
   In the energy-window case sigma is Eup. In fixed-band mode sigma is chosen
   inside the cached guard spectrum. The inertia count must cover the requested
   bands and be strictly smaller than the retained subspace.
7. If any check fails or iteration stagnates, solve the original physical H/O
   problem with GVD and refresh the cache. Diagnostics state the reason.

For a returned vector c, the residual criterion is

```
||Hc - e Oc||_2 / max(1, ||Hc||_2 + |e| ||Oc||_2) <= tolerance.
```

The inertia check matters: after a band crossing, an old invariant subspace
can have zero residual while omitting a new low eigenstate. The regression
suite constructs precisely this case and verifies dense fallback.

The current factorization is recomputed at each k point on each invocation.
It is still cubic, though cheaper than a complete dense eigensolution in many
regimes. This first version therefore does **not** remove all cubic cost or
avoid dense Hamiltonian/overlap construction. Factor reuse and cheaper LAPW
preconditioners are future optimizations after full SCF validation.

## Options and diagnostics

| Variable | Default | Meaning |
| --- | --- | --- |
| `DFTKS_DAVIDSON_TOL` | `1e-10` | Relative residual threshold above. |
| `DFTKS_DAVIDSON_BLOCK` | `16` | Maximum correction block size. |
| `DFTKS_DAVIDSON_MAXITER` | `40` | Maximum Rayleigh-Ritz iterations before dense fallback. |
| `DFTKS_DAVIDSON_MAXSPACE` | automatic | Maximum trial-space dimension; automatic is `min(N, max(2P, P+2*block))`, where P includes guards. A value not exceeding P triggers dense fallback. |
| `DFTKS_EIG_CACHE` | vector-file name | Cache filename prefix. Use a distinct prefix for each calculation/spin channel. Parent directories must exist. |

Each cache is named `PREFIX.dftks-r.k00000001` or
`PREFIX.dftks-c.k00000001`. It contains vectors and eigenvalues, not a dense
inverse. Writes use a PID-specific temporary file and atomic C `rename` on
POSIX systems. Each k point has one owner under existing MPI distribution.
The default prefix can change when vector-file rank suffixes change; this
causes cache misses, not reuse of mismatched data. An explicit stable prefix
allows reuse across MPI layouts. Do not run two calculations writing the same
cache prefix concurrently. Native binary caches are disposable and should be
regenerated after an incompatible compiler/ABI or format change.

`:DFTKS_SOLVER`, `:DFTKS_SOLVER_REASON`, `:DFTKS_DAVIDSON` and
`:DFTKS_SOLVER_WALL` identify the actual route, fallback reason, residual and
metric errors, and wall time including cache I/O. Invalid option values fail
explicitly. Production default remains GVX pending broader validation.

## Build and verification

Only `diagonalize.FPP` and its call in `lapw1.FPP` change existing production
source. The solver module is in the same compilation unit as the wrapper, so
no Makefile or `configure.py.gnu` edits are needed. Use serial `make` in this
directory because the existing build has incomplete module dependencies.

Standalone tests (NumPy, SciPy, gfortran and LAPACK required):

```sh
python tests/test_eigensolvers.py
```

On Linux the default link options are `-llapack -lblas`; on macOS they are
`-framework Accelerate`. Override `FC` and `LAPACK_LIBS` as needed. Tests build
in a temporary directory and enable Fortran bounds/runtime checks.

Full LAPW1 checks on copies of material inputs:

```sh
python tests/check_material_solvers.py CASE_DIRECTORY --executable ./lapw1c
python tests/check_material_solvers.py CASE_DIRECTORY --executable ./lapw1 --ranks 2
```

These copy structure/potential inputs, keep the first two k points, and compare
GVX, GVD, dense seeding, cached Davidson and AUTO. They independently read the
native vector files and dumped H/O to check residuals and normalization in
the actual output convention. The original calculation is never modified.

The initial 60 standalone cases cover synthetic real/complex generalized problems,
changing overlap, degeneracies, hidden band crossings, cache incompatibility
and truncation, ill-conditioning, forced nonconvergence, invalid settings,
and dense fallback. Real MnTe and complex displaced Ce samples agree with
GVX to about `1e-12 Ry` or better, including a two-process MnTe MPI check.
The material warm-cache calls hold H/O fixed; the standalone suite also
perturbs H and O between calls. These checks do not yet establish speedup in
large cells, convergence of full DFT+eDMFT cycles, or unchanged force/stress
accuracy. Those are the next acceptance tests before production use.
