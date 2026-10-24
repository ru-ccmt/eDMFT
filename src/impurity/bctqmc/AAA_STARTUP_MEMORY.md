AAA startup memory revision (2026-10-18)
========================================

No PARAMS changes are needed. OCA/TCA/FCA pole budgets, AAA tolerances,
balanced-channel selection, and residue/pole-refinement objectives are unchanged.

Implementation
--------------
* Ignore matrix-element rows that are EXACTLY zero at every conjugate-extended
  Matsubara sample. For the impurity adapter's seven-channel diagonal auxiliary
  matrix this means seven fitted elements instead of 49. There is no magnitude
  cutoff: arbitrarily small nonzero anomalous elements remain in the fit.
* Fill the real conjugate-weight Loewner matrix directly, without the two large
  complex intermediate matrices. Keep the public square-matrix representation
  and all symmetry/conjugation conventions unchanged.
* Allocate one maximum-sized real buffer and use block views at each support
  count. Perform column-pivoted QR in place; compute the smallest right singular
  vector from its small R factor, restoring the column permutation. This is a
  QR/SVD solve, NOT normal equations. A^T*A would square the condition number.
* Report per-process RSS and peak RSS; on macOS also report physical footprint.
  Diagnostics surround each QR/SVD and the nonlinear refinement. Only the MPI
  master fits and prints. Linux uses /proc/self/statm plus getrusage; unavailable
  process-memory fields are omitted rather than treated as zero.

The memory diagnostic helper is in aaa_memory.h. Copy it along with aaa.h and
imag_axis.h when transferring the updated source to another machine.

Why reuse matters
-----------------
The UTe2 beta=100 input is interpolated by the existing solver to beta=300,
giving 14,321 positive frequencies. Before this revision the final AAA iteration
had roughly 4.18 GiB of large matrices live. Isolated measurements also showed
memory growing across successive support counts; the pristine old fit exceeded
an 8 GiB watchdog limit and was stopped before finishing. Merely dropping the
zero rows and two intermediates still peaked around 4.0 GiB. Reusing the buffer
and doing QR in place removes this accumulation of growing tall workspaces.

Validation
----------
Workspace: /Users/haule/Documents/ChatGPT/HBCTQMC/aaa_startup_memory
Selected implementation: aaa.h matches aaa.inplace_v2.h byte-for-byte.
Selected executable: bctqmc.inplace_v2
SHA256: 2eb984cda51efea483b84983c03b4a9760eaed62e800ce98cc3f7eb6fe0a6bd1

Full UTe2 hybridization-only test (not a full impurity calculation):
* 14,321 positive samples; seven independent channels; all 25 support pairs.
* Complete joint + balanced fitting and 200-step refinements.
* Measured peak RSS 315.64 MiB; sampled peak physical footprint 271.78 MiB.
* Selected balanced, 23 real poles, relative RMS fit error 2.55648e-5.
* About 56 seconds in this test; concurrent development work was running, so
  this is not a controlled startup speed benchmark.
* See ute2_final.log. The full old fit was stopped safely; do not interpret
  these numbers as a completed old-versus-new full UTe2 solver benchmark.

A safe direct old/new comparison used the SAME 1,500-point subset of the
interpolated UTe2 data, retaining low frequencies and spanning the entire tail:
* Pristine old fit peak RSS 5,119 MiB; new fit peak RSS 55.27 MiB.
* Both retained 23 real poles.
* Matsubara fit difference: global relative 8.29e-6, worst channel 1.02e-5.
* Imaginary-time hybridization difference: global relative 1.12e-5,
  worst channel 1.74e-5 (1,001 time samples including both endpoints).
This subset is a validation fixture ONLY; production still uses every sample.

Two-rank fixed-input OCA solver comparisons used isolated copies of saved
single-band, normal cluster and SC cluster restarts (DLR tolerance 1e-7).
Global relative differences in analytical order-two output:
                  aStval.2        ag00.2
single band       5.86e-6          3.29e-6
normal cluster    9.59e-5          1.05e-5
SC cluster        1.31e-7          7.52e-8
SC worst pseudo-output column: 1.08e-5; physical G column: 5.18e-7.
Normal worst pseudo-output column: 5.21e-4; physical G column: 1.76e-5.
NCA order-one analytical outputs were identical.

These are recomputed AAA fits, not fixed-pole diagram comparisons. Removing
zero rows and factoring QR separately changes floating-point roundoff in the
near-null singular vector, which can change the subsequent finite-iteration
L-BFGS trajectory. The normal-cluster auto selector chose balanced rather than
joint; its relative fit RMS changed from 5.68818e-4 to 5.71611e-4. SC stayed at
23 poles and joint, with relative RMS 1.43344e-5 -> 1.43293e-5. No tolerance was
relaxed to obtain the memory saving. Outputs are close, not bitwise identical.

The *_checked.log reports use norms of EACH complex physical channel for
Gf.out/Sig.out, not separate relative errors of symmetry-zero real/imaginary
parts. This still protects small anomalous channels. Earlier *_validation.log
files retain the overly strict separate-component checks for auditability.

The existing AAA/DLR/OCA algebra tests also passed, including complex off-diagonal
matrix hybridizations, real-pole recovery, pole refinement, convolutions, and
physical/Bulla contractions. test_aaa_memory covers exactly-zero baths, all-zero
input, off-diagonal conjugation, and a nonzero channel of amplitude 1e-16.

Limits
------
This revision fixes the AAA startup allocation pattern. It does not assert a
316 MiB bound for the full UTe2 solver: pseudo-particle initialization, Dyson,
enumeration and final physical evaluation have their own memory requirements.
No full UTe2 solver job or self-consistency loop was launched for this change.
User PARAMS, hybridizations, restart files and existing results were untouched.
