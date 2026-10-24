Bounded-memory physical OCA (2026-10-18)
======================================

PARAMS defaults
---------------
  "oca_dlr_physical_mode": "auto",
  "oca_dlr_physical_batch_paths": 100000,
  "oca_dlr_batch_directory": ""

Mode is "auto", "batched", or "full". Auto selects batched physical OCA
when CIX has >=4096 symmetry-reduced pseudo propagators. This is a global,
rank-independent heuristic, not a measurement of available RAM. UTe2's
19685 pseudos select batched automatically. Smaller normal/SC cluster
cases keep the existing full-plan path; explicitly select batched to test
them. Self-only iterations, NCA/TCA/FCA, Dyson, and legacy OCA are unchanged.
Legacy-comparison diagnostic mode deliberately uses the full path.

Memory design
-------------
The previous physical transition retained the self-only plan while building
one large measurement hash table, weights, compact paths, and new plans.
The new path drops the old plans first. It walks the same atomic contributions
using a fixed-capacity hash table and weight vector. When full, it spools
partial sums to private local files grouped by reduced bath pair. There is
no combinatorial count pass and no dynamically growing global map in this
mode; the original full path still uses its two-pass algorithm.

Each bath-pair file is reduced by COMPLETE packed key. An oversized file is
partitioned into two children by successive bits of a full-key hash. This
revisits a disk file, not the atomic diagram loops. All duplicate weights meet
before a key is evaluated, including physical-only Bulla paths. Fully zero
records are removed only after summing. Floating-point addition order can
change; no poles, paths, observables, or importance rules are approximated.

Each reduced batch builds compatible self/physical plans, evaluates ALL poles,
adds into the DLR node outputs, and frees its records/plans. The final summed
nodes are interpolated just once. At most batch_paths distinct keys live in
one enumeration/reduction map or plan. This caps the growing topology part,
NOT total process memory in bytes: pseudo arrays, DLR rank, pair tables,
compiler/allocator overhead, sparse construction, and MPI buffers still count.
Lower batch_paths if needed; higher values may reduce repeated preparation
and disk traffic but consume more RAM. Very small values can be much slower.

Physical numerical scratch is improved in both modes: one reflected pair
family at a time, and row-batched sparse contractions controlled by the existing
oca_dlr_self_scratch_mb (default 64). Zero disables row batching. Physical
evaluation still needs both original pair families; bounded plan size limits
their growth. The self-only pair-streaming optimization remains intact.

MPI and disk
------------
Whole reduced bath pairs remain assigned to ranks by the existing cost model.
Enumeration, spooling, partitioning, plan building, and all-pole evaluation
are LOCAL. There are no path gathers or batch-by-batch collectives. One final
reduction phase comprises four MPI_Allreduce calls: Sigma nodes, G, F00, F01.
Ranks without bath pairs participate only at the final reduction. Thus this
mode uses at most the number of reduced bath pairs (49 for this UTe2 input)
for physical work; local hash partitions are not dynamically dispatched.
oca_dlr_mpi=false retains replicated evaluation and performs no final sum.

Set oca_dlr_batch_directory to an EXISTING node-local disk scratch directory
on Beowulf. The empty default uses std::filesystem::temp_directory_path()
(normally TMPDIR or /tmp). Do NOT choose RAM-backed tmpfs: that defeats the
memory saving. Prefer local disk over a slow shared filesystem. The setting
is a literal path, not a shell expansion; use a resolved path in PARAMS.

Each rank creates a unique oca-physical-XXXXXX child directory. Files are
temporary native-format working data, not restart files. The .self persistent
cache is unchanged; the old .physical cache is bypassed in batched mode because
its loader materializes an entire shard. This implementation does not yet
persist the bounded batches between solver invocations.

Writes are checked, and free space is checked periodically with 1 GiB headroom
(not an MPI-wide reservation). A write/quota/disk-space error aborts the MPI
job rather than leaving ranks hanging at reduction. Scratch is removed on
normal completion or C++ exception. SIGKILL or reboot cannot run cleanup;
each rank prints its private scratch path so leftovers can be identified.
Never remove a shared scratch root indiscriminately.

Progress reports include partial records and bytes written during enumeration,
then completed local batches and unique paths during evaluation. Reports occur
at flush/batch boundaries approximately every 60 seconds. Rank 0's log is local
progress, not a global percentage. End-of-work reduction time includes waiting
for slower ranks, not just network transfer.

Scope
-----
This fixes the unbounded physical-enumeration/evaluation storage. It does not
fix or explain the separate ~48 GiB initialization peak, and does not establish
that a strict 2 GiB-per-rank job limit is safe. The full UTe2 physical run must
still be measured. Keep the saved PPSigma restart; do not overwrite completed
results for a benchmark.
