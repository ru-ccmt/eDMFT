Enumeration-free block-matrix OCA (default)
==========================================

PARAMS:
    "OCA_solver": "dlr",
    "oca_dlr_backend": "matrix"

The default is "oca_dlr_backend": "matrix". To use the previous evaluator,
set "oca_dlr_backend": "paths" explicitly. This restores the established
full/bounded scalar-path evaluators and their disk topology caches. No cache
is deleted by the matrix backend. An existing explicit "paths" setting still
overrides the default; removing the setting now selects matrix OCA.

The matrix backend computes pseudo self-energy and, when requested by the
driver, physical G and Bulla F00/F01. It uses the same OCA AAA fit, DLR basis,
signed Sindex, pseudo degeneracies, bath degeneracies, Q normalization, and
Nambu conventions. It does not borrow TCA's reduced-pole fit or change Dyson.
No atomic-state truncation or hybridization approximation is introduced.

Geometry and signs
------------------
The existing Phi representative has forward lines (0,2),(1,3). Rotating its
four pseudo cuts to an open chain gives root/split direction pairs 00,01,10,11.
All four are evaluated; there is no additional factor 1/2. The root is (0,2),
the split line (1,3). One crossing supplies -1, and each Delta uses the shared
engine's -residue convention. Physical cuts place creation at zero and
annihilation at two and give the two split directions 0 and 1.

The root bath is summed inside a two-sided sparse-operator contraction. The
other bath and one AAA pole index label independent tasks. Two matrix-valued
DLR convolutions replace each scalar path integral. Physical output uses two
independently convolved arcs and sparse boundary traces, with the FULL dense
Bulla block (not TCA's fixed-transition-support projection).

Each task retains a few block-matrix functions; it never creates the enormous
four-pseudo path map or scalar-pair catalogs. Signed Sindex labels reconstruct
matrix entries from the already-shared unique pseudo coefficients. The first
implementation intentionally keeps the four cuts explicit for auditing.

MPI
---
Tasks (cut,pole,split bath) are assigned round-robin. There are 4*L*Nbath self
tasks and 2*L*Nbath physical tasks. For 23 poles and 14 baths: 1288+644 tasks.
Only final compact DLR node arrays are summed across ranks. No communication
occurs inside convolutions; no topology or pair tables are exchanged.
"oca_dlr_mpi": false computes the complete answer independently on each rank.

Fallbacks and limits
-------------------
An active G_importance/OCA_cut_sum filter generally cannot be applied after
matrix summation. If 4*max(G_importance)>OCA_cut_sum, select the established
path backend with an explicit log message. Legacy diagnostics also fall back
to the established path workflow. A conservative fallback may occur even
when no surviving path would actually be cut; it never silently ignores cuts.

oca_dlr_physical_mode (auto/full/batched), topology sharding, whole-pole teams,
and disk enumeration caching apply to the paths backend, not matrix tasks.
All existing parameters must still be valid. The matrix backend's static
scheduler is deliberately simple; its speed depends on atomic block sizes.
There is no automatic switch based on model size at this stage. Matrix OCA
is the default even for small normal-state clusters, where paths can be faster.

Validation/benchmark tools
--------------------------
check_solver.py creates isolated single/normal/SC fixtures from starting_OCA
and compares all self-energy/physical solver files to the paths reference.
bench_oca_matrix.cc reads CIX/Sindex/saved PPG and evaluates a bounded sample
of tasks. A pole table can supply a real fit; otherwise the benchmark clearly
labels its four poles synthetic. It never changes production inputs/results.
See OCA_MATRIX_VALIDATION.txt for measured comparisons and timing limitations.
The isolated test driver and its logs are retained in
/Users/haule/Documents/ChatGPT/HBCTQMC/oca_matrix_prototype/.
