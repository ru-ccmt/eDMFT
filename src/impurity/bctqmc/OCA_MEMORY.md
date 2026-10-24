OCA bounded self-energy scratch and streamed pair families
========================================================
2026-10-18

Defaults (PARAMS entries are optional)
  "oca_dlr_self_scratch_mb": 64.0,
  "oca_dlr_stream_pairs": true

For the previous full-table evaluation/lifetime policy:
  "oca_dlr_self_scratch_mb": 0.0,
  "oca_dlr_stream_pairs": false

No approximation, pole truncation, or topology cutting is added. NCA, TCA,
FCA, Dyson, and the existing OCA MPI ownership/reduction policy are unchanged.
The compiled OCA plans still stay cached when N_min_order == Norder.

Outer contraction batching
--------------------------
The default whole-pole/serial evaluator processes bounded contiguous rows of
each cut's outer-pair table. It performs the bath contraction and DLR solve
for that batch, then the final convolutions and accumulation. Contiguous
equal-outer-propagator subranges share their left convolution factor. A batch
may split an exceptionally large group without changing its formula.

Sparse self-energy contractions are stored row-major instead of column-major
so a row batch touches only its own nonzeros. No second sparse copy is kept.
This storage choice also applies to the fallback. Thus the fallback restores
the previous full-table algorithm, not necessarily bitwise identical rounding
to the older executable; comparison with that executable is tested separately.

The scratch target allows twelve rank-by-batch double matrices for temporary
contractions, DLR transforms, and the current convolution implementation.
At DLR rank 49, the default gives 14266 rows per batch. At least one row is
always allowed. Targets must be finite and nonnegative.

This is NOT a total-memory limit. Persistent plans, sparse construction work,
pair tables, MPI buffers, input/output arrays, and allocator/linear-algebra
overhead are outside the target. The diagnostic cooperating-team path
(oca_dlr_whole_poles=false with team size >1) keeps its existing per-outer
block contractions; it does not use this new batch-row target.

Pair streaming
--------------
Self-only calls build catalog 1, use it for cuts 0 and 1, and release it
before building catalog 0 for cuts 2 and 3. Each catalog is still constructed
once per pole, with the same factor dressing and convolution. MPI team ranks
follow the same collective order, including empty shards/catalogs.

Physical calls need both families simultaneously, so streaming is disabled
automatically there. The self-energy portion of such calls still uses bounded
outer batches. Physical reflection/mapped-table batching is a separate future
optimization and is NOT part of this change.

The log prints the scratch target, whole-pole row limit, and effective pair
streaming policy at every OCA evaluation. Restart with the newly installed
executable to use this implementation; replacing the executable does not
change code already loaded in an existing process.

Validation and limitations
--------------------------
Staging tests/logs:
/Users/haule/Documents/ChatGPT/HBCTQMC/oca_memory
See validation.txt for completed tests and measured memory/timing results.
No full UTe2 run is launched by this installation. Its actual memory reduction
and speedup must be measured on the next production run; avoiding compression
pressure can help, but batching is not guaranteed to improve CPU time.
