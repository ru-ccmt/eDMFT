# dftKS portability and MPI threading

`lapw1` and `lapw1c` distribute k-points over MPI ranks. Each rank keeps its
own dense Hamiltonian, overlap matrix, eigenvectors, and LAPACK workspace.
Threading BLAS independently on every rank therefore causes CPU
oversubscription and, more importantly, multiplies temporary memory use.

For MPI runs, `mpi.FPP` now calls `limit_blas_threads` before the first BLAS
operation. The implementation in `cputim.c` sets the common OpenMP,
OpenBLAS, MKL, BLIS, and Apple Accelerate thread controls to one. Serial runs
are unchanged. To deliberately use more than one BLAS thread per MPI rank,
set, for example:

```sh
export EDMFT_BLAS_THREADS=2
```

Use this only when the product of MPI ranks and BLAS threads fits the
available cores and memory.

The largest matrix allocation now reports its dimensions, approximate
per-rank storage, and the allocator error instead of failing without useful
context. LAPACK workspace allocation and query failures are also checked.

GNU optimization settings are maintained in the repository-level
`configure.py.gnu`. Do not hand-edit generated `Makefile.in` compiler flags;
run the normal `setup.py` configuration process after changing the
configuration file.

The current `dftKS/Makefile` does not describe all Fortran module
dependencies. Consequently, a clean parallel build such as `make -j4` can
race while producing `mpi.mod`. Use a serial `make` for this directory until
the generated build system gains explicit module dependencies.
