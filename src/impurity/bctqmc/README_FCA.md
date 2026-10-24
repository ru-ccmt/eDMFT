# FCA: fourth-order DLR reference implementation

This implementation adds all six supplied fourth-order Phi skeletons (A-F).
It is a matrix evaluator with optional exact cut reduction and physical-chain
reuse. Four-pole normal-state self-energy comparisons against fresh MC are
encouraging; pole convergence and statistical errors still require care.
NCA/OCA/TCA remain separate; the new sum-reuse optimizations are FCA-only.

## Closing-vertex orientation sums and E/F physical presumming

Two independently switchable paths extend the previous prefix-cache work:

```json
{
  "fca_self_direction_sum": true,
  "fca_physical_ef_presum": true,
  "fca_physical_ef_boundary": false
}
```

The self-energy task now fixes two internal arrows and explicitly evaluates
BOTH orientations of one selected line. Its sparse opening/closing operators,
positive/negative-frequency decorations, residues and kernel normalization
are applied separately. Only after that line closes are contributions with
the same remaining bath labels merged. The common closing convolution and
suffix run once. This is linearity, not a symmetry assumption or removal of
an orientation. Grouped-pole tasks combine this with their existing pole sum;
the two-pole case still respects `fca_self_pole_sum_min_poles=3` and groups only
arrows, not poles. Shared A--D families likewise merge both first-line arrows
before each member's own closing convolution. The MPI task count halves, but
the expensive interval is still evaluated for both arrows: this does NOT
imply a factor-two speedup. Prefix-budget failures explicitly sum the old
single-arrow paths. Disable the new switch for the previous schedule.
For independent-pole tasks the existing adaptive bath anchor is preserved;
previously unstreamed geometries keep two independent arrow evaluations.
This prevents grouping from making the middle of a chain more expensive.

E/F physical families use the existing canonical crossing-bath-pair traversal.
The q=3 E members' closed right-arc inputs are added with their INDIVIDUAL
normalizations and cut multiplicities before the common final convolution.
F remains separate to preserve topology-resolved output. The short-arc
operation/decorations and final-edge equality are checked before reuse.
With `fca_physical_ef_boundary=true`, external-bath-summed boundary weights
are reused across the family: `(C L A)^T`
for a shared left arc, or `(A R C)^T` for the mirrored q=5 shared right arc.
G, F00 and F01 retain the same connector, channel and sign conventions.
This extra boundary preparation is opt-in: E/F families have only two or
three members, so its setup cost is not necessarily amortized. It is separate
from the existing A--D `fca_shared_boundary` optimization, which stays enabled.

The closed E accumulator reserves part of the existing physical shared-cache
allowance; it does not allocate a bath-pair table. Insufficient space uses
separate final convolutions. Boundary weights respect `fca_boundary_cache_mb`
and fall back before accumulating any output if admission fails. The log
reports saved final convolutions and boundary builds/hits. No AAA fit, DLR
tolerance, cut factor or production input is changed by either optimization.

Coefficient-space convolution chaining is NOT implemented by these changes.

## Exact self-energy pole sums and shared input transforms

Two additional optimizations are enabled by default:

```json
{
  "fca_self_pole_sum": true,
  "fca_self_pole_sum_min_poles": 3,
  "fca_shared_branch_transforms": true
}
```

`fca_self_pole_sum` schedules a rooted cut, internal arrows, and TWO fixed
poles. The third pole and its bath are streamed through that line's interval.
Its prefix is evaluated once; contributions are merged at the closing vertex
with all remaining bath labels still distinct; its suffix is evaluated once.
The pole-dependent kernel normalization is applied BEFORE this merge.
Positive-frequency edge factors and negative-frequency endpoint factors are
both included. No frequency, bath, or diagram is dropped. This changes the
sum order, not the AAA fit or the DLR tolerance.

Pole grouping activates automatically at three or more poles. A bounded SC
two-pole benchmark found that changing the streamed bath could INCREASE work,
especially in topology F. One/two-pole calculations therefore keep independent
triples by default, while still allowing shared transforms. Setting
`fca_self_pole_sum_min_poles=1` explicitly forces grouping for convergence and
performance experiments. The startup `self_pole_sum=yes/no` reports the actual
choice, not merely the requested switch. The direct benchmark group method
always exercises grouping, independently of this production guard.

With cut reduction there are `216*p*p` self tasks, compared with the previous
`216*p*p*p`. Each grouped task is larger; the expensive middle still scales
cubically in poles. Thus task-count reduction is NOT the speedup. Physical
tasks remain `216*p*p*p`: their complementary-arc pole dependencies must not
be prematurely summed. The existing dynamic MPI dispatcher is retained.

The streamed line is chosen from estimated live bath-state work. At most two
variable bath labels are stored; no three-pole/bath tensor is materialized.
If the retained prefix exceeds `fca_self_prefix_mb`, the task falls back to
independent pole triples with their original adaptive bath-reuse plan, not
directly to the slower Cartesian reference. `fca_self_reuse=false` also uses independent triples
inside the grouped task. For the completely previous schedule and evaluator,
set BOTH new switches false.

`fca_shared_branch_transforms` uses the identity `P[H(tau) F_b] = P[H] F_b`
for the time-independent sparse bath operator. The scalar endpoint decoration
is included before the common interpolation/Hilbert transforms. Sparse
operators act on both transformed arrays, which feed the subsequent convolution
without another input transform. The output remains node values: no
coefficient-space chaining or replacement of the interpolation solve is used.
The streamed bath shares this transformed prefix once per pole as well.
This describes the branch-transform optimization, not the separate physical
coefficient-chaining option described below.

`fca_physical_coeff_chain=true` now keeps intermediate FCA physical arcs in
DLR coefficients and returns node values only at their final convolution.
It covers ordinary arcs, shared A--D prefixes and pre-summed E/F tails.
The default is true, but it is enabled only for `dlr_tolerance >= 1e-7` and
a successful finite precomputed transform. Complete SC bath sums exposed
sensitivity at tighter tolerances that isolated chains did not reveal;
those cases automatically use the previous node path. Set this switch to
false for the previous path at any tolerance. Pseudo self-energy, NCA/OCA/TCA,
pole counts and MPI scheduling are unchanged. See `FCA_COEFF_CHAIN.txt` for
the algebra, representation boundaries, numerical safeguards and tests.

An experimental extension, `fca_extended_branch_transforms=true`, applies this
identity to the shared A--D self-energy family, the shared physical family,
and internal-arc bath openings of the ordinary physical evaluator (including
E/F). **The extension defaults to false:** initial normal/SC cluster samples
did not establish a worthwhile overall speedup. The older ordinary-self
optimization remains enabled as before. Closures still sum node values before
convolution; there is no independent summation of correlated crossing baths.
Setting `fca_shared_branch_transforms=false` also disables the extension and retains the previous node-based
path in all these evaluators, without disabling the independent prefix/arc
sharing switches.

In the shared physical family, coefficients and their convolution-Hilbert
transform are cached for the common input to the innermost bath branch.
The cache holds immutable input owners as well as the two transformed
functions, so pointer identity cannot be invalidated by node-cache eviction.
It is local to one pole/direction task and its accounted retained bytes are
bounded by `fca_physical_cache_mb` (64 MiB by default), ADDITIONAL to the
`fca_physical_shared_mb` node-prefix cache and the prepared-propagator cache.
Zero `fca_physical_cache_mb` disables this transform cache. Live branches,
convolution scratch, and transient packing arrays are not a total-RSS budget.
Outer bath-prefix outputs continue to use the node cache; their transformed
inputs are not retained because distant sibling visits caused cache thrashing
in initial timing tests.

Shared-family logs report `shared_branch_inputs`; shared physical logs also
report `branch_transform_hits` and `transform_cache_peak_mb`. Ordinary physical
profiles report `physical_shared_branch_inputs`. These count actual reuse-path
work, not the nominal number of baths. Bounded real-input A/B modes are
`bench_fca CIX Sindex PPG poles branch-prefix 2` and `branch-physical 2`;
both hold the existing family-sharing algorithm fixed and toggle only the
common-input transform, alternating benchmark execution order by sample.

Initial extension validation (2026-10-16): saved normal-cluster CIX/Sindex/PPG
and its actual four-pole FCA fit, plus saved superconducting blocks with a
signed synthetic two-pole fit. These are bounded single-rank samples, NOT
full solver-phase or many-rank scaling measurements. The synthetic SC fit
tests the algebra, not the physical accuracy of using two hybridization poles.

| Sample | Previous seconds | Extended seconds |
| --- | ---: | ---: |
| Normal shared self, two fixed-pole pairs | 17.503 | 18.242 |
| SC shared self, two fixed-pole pairs | 56.230 | 57.334 |
| Normal shared physical, two triples | 9.276 | 9.343 |
| SC shared physical, two triples | 91.941 | 95.917 |
| Normal E/F ordinary physical, eight cuts/sectors | 5.163 | 4.977 |
| SC E/F ordinary physical, eight cuts/sectors | 60.929 | 59.580 |

Small percentage differences are not strong timing evidence; some compiler
and validation processes were running concurrently. They certainly do not
justify enabling another per-rank cache by default. The initial unrestricted
physical transform cache was about 12% slower in both cluster samples; the
retained implementation limits it to the inner bath's adjacent siblings.
Largest normwise A/B difference in the retained implementation was 9.7e-11
(in tiny SC physical contributions); self samples differed by at most
6.5e-15. Both Bulla components are compared, not just G or Phi.
The final one-band solver regression compared 68 diagram files and the final
physical outputs with the previous executable at 1, 2, and 4 ranks. It also
checked disabled, zero-budget, tiny-budget, bounded-cache, and non-family
routes. Worst diagram-file difference was 3.3e-13; the default two-rank route
and enabled-family two-rank route matched the baseline numeric tables exactly.
All analytic, mixed-pole, reflected-convolution, and quadrature tests passed.

Sharing is automatic only when compatible branch matrix entries exceed 1.5
times incoming entries, and an estimate including target block sizes, cached
prefix coefficients, and temporary arrays fits `fca_self_prefix_mb`.
Otherwise the previous node-based batched transform is used. This budget
controls additional reuse buffers, not total solver memory. Caches are local
to a task and cannot leak across pseudo-self-energy iterations.

The shared convolution API additionally handles an all-zero right propagator
by clearing its output (rather than returning unmodified input). Existing
nonzero TCA calculations keep their original numerical path.

The startup log prints both switches. Self profiles add `pole_groups`,
`shared_branch_inputs`, and `shared_transform_seconds`. The latter is a subset
of vertex time; it is not included in the downstream convolution-transform
time. All profile seconds are worker sums, not elapsed wall time.

For bounded validation, `bench_fca ... sums` compares 48 pole groups across
all six topologies; `sumsquick` checks one group per topology; `sumsall` checks
all 384 rooted-cut/arrow sectors (intended for one band). Each compares the
previous sum of independent triples, pole grouping alone, and both changes,
with additional memory-fallback and map/no-scratch checks. These are sampled
groups, not complete production evaluations.

Validation of the sum-reuse change (2026-10-15):

- All 384 one-band cut/arrow sectors, mixed-sign two-pole fit: maximum
  relative group difference `1.37e-9` (includes small cancelling groups).
- Complete one-band sum versus the Cartesian reference, including physical
  G/F: maximum relative difference `1.84e-14`.
- Analytic simplex and mixed-sign/zero-frequency quadrature tests retained
  their previous accuracy. Prepared-left rectangular convolutions agree
  exactly in the reflection/scratch unit test; the all-zero right-block
  regression also passes.
- Isolated full solver A/B runs: 68 order/topology/cut output files agree
  within `1.57e-12`, including forced fallback and 1/2/4 MPI ranks. Final
  physical outputs are unchanged at printed precision for two ranks.
  Final four-pole solver runs (automatic grouping enabled) agree within
  `3.36e-13` across the same switches/rank counts. One/two-pole default runs
  explicitly report `self_pole_sum=no`; forcing the minimum to one reports
  `self_pole_sum=yes` and also passes the comparison.
- Normal cluster: 48 four-pole groups on saved `FCA_partial/PPG.0` give maximum
  relative difference `1.12e-10`. Worker sums were 170.66 s before, 151.46 s
  with pole grouping, and 153.66 s with both optimizations. This is about
  10% less time overall in that sample, not a full production wall-time
  prediction. The extra shared-transform benefit is not resolved by that
  timing sample; its effect depends on block shapes and memory fallback.
- SC cluster: six cut/arrow groups on the saved real matrix-valued PPG and a
  synthetic signed two-pole fit agree within `6.67e-14`. With the optimized
  fallback, forced pole grouping still took 129.17 versus 124.56 worker seconds
  in this sample. This is the reason for the default minimum of three poles;
  it is NOT claimed as an SC speedup. It is a bounded topology check, not a
  full SC evaluation or a test of physical accuracy of a two-pole fit.
- A complete three-pole one-band comparison (negative, zero, positive
  frequencies) agrees with the Cartesian reference within `3.02e-15`.

## Running one deterministic evaluation

Start from a COPY of a converged restart and its matching Delta/CIX/Sindex:

```json
{
  "N_min_order": 4,
  "Norder": 4,
  "Nitt": 1,
  "OCA_solver": "dlr",
  "TCA_solver": "dlr",
  "tca_accuracy": "veryhigh",
  "FCA_solver": "dlr",
  "fca_accuracy": "veryhigh",
  "fca_dlr_mpi": true,
  "fca_progress_seconds": 60.0,
  "fca_store_self_diagnostics": true
}
```

These are additions/overrides to the model's existing PARAMS, not a complete
model input. All lower orders are included once. Physical G/F are computed
only on the existing final-measurement iteration. With Nitt > 1 the usual
dense-mesh Dyson iteration now includes fourth order.

`fca_accuracy` controls **only FCA**: its hybridization pole cap and the DLR
basis used to evaluate fourth-order diagrams. OCA, TCA, and the dense-mesh
Dyson solver keep their own settings. The presets are:

| `fca_accuracy` | FCA AAA poles | FCA DLR tolerance | Purpose |
| --- | ---: | ---: | --- |
| `low` | 3 | `1e-4` | Fast exploratory SC calculations |
| `modest` | 4 | `1e-5` | First production candidate; validate for each model |
| `high` | 6 | `1e-6` | More accurate fourth-order check |
| `veryhigh` | all available | `1e-7` | Reference calculation, potentially very expensive |
| `custom` | `fca_aaa_poles` | `fca_dlr_tolerance` | Explicit independent controls |

The default remains `standard` for compatibility with existing PARAMS files.
It preserves the old behavior: `fca_aaa_poles=0` uses all available poles
and `fca_dlr_tolerance=0` reuses the shared OCA/TCA DLR basis. `standard`
also accepts manually set values, but new inputs should use `custom` when
overriding either control. For example:

```json
{
  "fca_accuracy": "custom",
  "fca_aaa_poles": 5,
  "fca_dlr_tolerance": 1e-6
}
```

In `custom`, `fca_aaa_poles=0` means the full accurate AAA pole table, while
`fca_dlr_tolerance=0` means the shared DLR basis. A positive pole count caps
an independent channel-balanced FCA fit; a positive FCA tolerance creates an
FCA-only basis if it differs from the shared tolerance. Named presets replace
both low-level controls, even if they also appear in PARAMS. The selected
settings and actual pole count, fit errors, DLR tolerance and rank are printed
in `ctqmc.log`. Pole count and tolerance remain fixed across self-consistency
iterations; there is no implicit last-iteration promotion. Do not change
them in PARAMS while a job is running.

These levels are **cost/accuracy targets, not universal error guarantees**.
The SC three-pole hybridization fit had a large relative fit error, and its
good rank-22 versus rank-37 comparison tested the DLR reduction at fixed
three poles, not the pole approximation against a full-pole answer. Because
three hybridization lines are split, pole work grows roughly as the cube of
the pole count; `veryhigh` can be impractical for a large SC problem.

Compression minimizes the worst relative channel fit error using the existing
greedy residue-refit algorithm. When the source has both frequency signs and
the requested count is at least two, both signs are retained. Per-channel
fit errors and the actual frequencies are printed. This does NOT guarantee
that a two-pole fit is physically adequate, particularly for anomalous channels.
The unexpanded root uses the SAME fit as the three expanded lines.

For an independent MC reference use `N_min_order=2, Norder=4` (OCA anchor)
or `N_min_order=3, Norder=4, TCA_solver="dlr"` (TCA anchor), with the
appropriate MC sampling parameters. For fifth-order MC use
`N_min_order=4, Norder=5, TCA_solver="dlr", FCA_solver="dlr"`.

`N_min_order` identifies the last analytic order. MC visits that order as
its normalization sector and all higher orders through `Norder`; only the
higher sectors are added to the analytic result. The sampled anchor is still
written to `Stval.<order>` and `mg00.<order>` for comparison with its analytic
counterpart. Analytic TCA/FCA and MC are therefore not mutually exclusive.

In mixed analytic/MC runs, TCA and FCA use the full fitted AAA representation,
all pseudo-particles and no frozen tail. TCA accuracy presets and both FCA
accuracy controls are available for purely analytic runs, but are bypassed
for the fourth-order normalization anchor when fifth-order MC is enabled:
an approximate anchor would rescale every higher-order MC contribution.
FCA then uses the shared DLR basis. Invalid/zero analytic or
sampled normalizations are rejected rather than producing zero/NaN output.
An unavailable analytic anchor (e.g. `N_min_order=5`, or `TCA_solver="mc"`
with `N_min_order=3`) is rejected before sampling.

## Output

`aStval.4` is now saved immediately after all FCA pseudo-self-energy tasks
and MPI reductions finish, BEFORE physical G/F evaluation. Wait for
`FCA self checkpoint saved: file=aStval.4 ...` in the log, not merely the
self-task 100% progress line. The file is closed and atomically renamed
before that marker. Phi and optional topology/cut diagnostics are available
at the same phase boundary. This also applies on pseudo-only iterations.

Stopping the job after this marker preserves the complete FCA self-energy
for comparison with `Stval.4`. It does NOT finish the Dyson iteration or save
a resumable physical-phase checkpoint. No new `ag00.4` is produced until the
physical phase completes; any existing `ag00.4` may belong to an older run.
Use a separate test directory or check timestamps to avoid mixing runs.

- `aStval.4`: fourth-order pseudo self-energy only, on the usual dense mesh.
- `ag00.4`: fourth-order physical G, F00 and F01, in Print_Physical format.
- `aStval.4.topologyA` ... `topologyF`: topology-resolved pseudo self-energy.
- `aStval.4.topologyA.cut0` ... `.cut7`: all eight cuts of each topology.
- `ag00.4.topologyA` ... `topologyF`: topology-resolved physical estimators.
- Logs include all six Phi contributions, cut-sum checks, D/F C4-orbit
  differences, phase progress and per-iteration `t(FCA)`.

Topology/cut files require `fca_store_self_diagnostics=true` (reference default).
Physical files are written only on measurement iterations. As with the older
orders, these filenames describe the latest evaluation, not iteration archives.
The cut files ALREADY include the topology symmetry factor. Their sum is the
topology file. Reduced D uses 8*weighted_cut0; F uses
4*(weighted_cut0 + weighted_cut1). The saved cut files already contain their
original Phi symmetry factors, which must not be applied a second time.

## Algebra and implementation

`fca_topology.h` contains the edge matchings and cut rotations. `fca_dlr.h`
contains the arbitrary-pole evaluator and integer-task MPI scheduler. Matrix
convolution, signed atomic operators, root contraction, Sindex projection,
and physical/Bulla closure reuse the tested TCA primitives. The only TCA
header change exposes its order-independent physical boundary helpers.

The all-arrow Phi factors are [1/2,1,1/2,1/8,1,1/4] for A-F. Every pseudo
cut is retained. The crossing parities [-,+,-,+,-,+] are checked against the
generic MC SignSort routine. These are separate from atomic and bath signs.
`fca_reduce_cuts=true` now groups identical rooted chord geometries under
rotation. The self cut representative counts for A-F are 4,8,4,1,8,2.
Physical cuts include the external orientation in their geometry and have
the same representative counts. Internal bath, direction, and pole sums
remain complete; the original symmetry factors are multiplied by orbit
multiplicity. This never identifies actual atomic states via Sindex.

For pseudo self-energy, rotate the cut to the closing edge, leave the line
incident at vertex zero unexpanded, and split the other three. For effective
frequency kappa > 0, distribute K over every edge inside the line; for
kappa <= 0, use the two-endpoint identity. Denominators always use a
nonnegative frequency. The root bath and its two orientations are contracted
inside the matrix chain. Split bath labels are summed immediately at closure;
the reference streams the latest-closing label so at most two other labels
remain live. The default optimized evaluator instead shares the prefix before
the streamed line opens, streams only its bath-dependent interval, then sums
that bath BEFORE its closing-edge convolution and common suffix. Other open
bath labels remain separate until their own closing endpoints. This is exact
linearity; no pseudo states, diagrams, orientations, or poles are discarded.

The self-energy controls (all default to the optimized path) are:

```json
{
  "fca_self_reuse": true,
  "fca_self_adaptive_anchor": true,
  "fca_reuse_prepared": true,
  "fca_self_prefix_mb": 64.0
}
```

`fca_self_adaptive_anchor` compares the three possible streamed labels using
a geometry-based estimate of bath-state operations and root contractions.
Equal-cost plans prefer fewer live states and then a smaller retained prefix;
this prevents an arbitrary label ordering from retaining a large SC prefix.
It also considers an unstreamed chain when no more than two labels are ever
open. The selection is a heuristic, not a measured timing guarantee; it never
permits an Nbath^3 state table. Set it false to retain the latest-closing
anchor while still sharing prefix/suffix work.

The retained prefix is checked against `fca_self_prefix_mb` before entering
the streamed interval. If too large, it is released and the reference evaluator
is used. This is a bound on additional retained prefix storage, NOT on total
process memory: live states, accumulated suffix inputs, coefficient blocks,
and convolution buffers also consume memory. No prefix persists across tasks
or self-consistency iterations. A zero budget forces fallback for streamed
plans; unstreamed plans do not require a retained prefix.

`fca_reuse_prepared` shares dressed pseudo coefficients and their lazy block
caches between edges whose sampled decoration vectors are exactly identical.
It applies to both optimized self and physical evaluators. It uses no
tolerance-based merging and does not cache across pole triples or iterations.

For physical cuts, remove one Delta and split the remaining three. Positive
effective frequencies decorate the direct arc; negative frequencies decorate
the complementary arc. Evaluate the tau and beta-tau chains, reflect the
second and close ordinary/Bulla estimators with a shared external bath sum.
With `fca_physical_reuse=true`, a bath entirely inside one arc is summed
at its closing endpoint before subsequent convolutions. Baths joining the
two arcs are kept correlated until the ordinary/Bulla boundary contraction.
Prefixes depending on fewer than all crossing baths are shared. No complete
three-bath chain table is cached. Internal states are convolved in batches;
the final convolution of the second arc also performs its reflection.

`fca_physical_cache_mb=64` bounds the retained prefix cache per MPI worker
and per task (not total process RSS). Two independently testable options
improve its use without increasing this budget:

- `fca_physical_cache_resume=true`: find the furthest available prefix before
  evaluating an arc. This avoids rebuilding a large missing intermediate
  immediately before hitting a later, already bath-summed cached result.
- `fca_physical_cache_retain=true`: a later prefix replaces earlier prefixes
  with exactly the same crossing-bath dependencies. At a dependency boundary,
  a new prefix may evict least-recently-used bath-dependent entries; temporary
  branches cannot evict those entries. Bath-independent seeds are protected.
  If the later prefix cannot fit, earlier fallbacks remain available.

Both default to true. Retention requires resumption; setting resume to false
restores the old forward traversal and fill-once policy. Zero cache budget
disables storage, not exact early bath sums or batching. Cache lifetimes are
bounded by a task: nothing persists into another iteration or later MC work.

`fca_physical_ef_shared=true` groups matching E/F physical cuts with q=3 or q=5
and the same canonical crossing-line mask. It is independent of the older
q=4 A--D flag `fca_physical_shared`. Canonical line labels follow increasing
opening vertex; both directions and pole labels are mapped back for the
ordinary fallback. All orientations and orbit multiplicities are retained.
Within one crossing-bath pair, matching two-stage prefixes are reused across
members. The first crossing label encountered by the longer arc runs in the
outer bath loop; this order reverses between the q=3 and q=5 families.
Cache keys compare the complete operator history, internal pole
indices, sampled edge decorations, and reflection flags, so this includes
negative-pole decorations without assuming they equal positive-pole ones.
Global normalization and crossing-line residues remain outside these cached
values. Full atomic blocks are preserved; Sindex is not an equivalence test
for diagram paths.

The E/F cache budget combines `fca_physical_shared_mb` (128 MiB by default)
and `fca_physical_cache_mb` (64 MiB). Half the shared allowance is reserved
for the current bath-pair's cross-diagram prefixes; the other half supplements
the ordinary allowance, which is split between members. This lets each member
retain both a branched seed and a later bath-summed prefix between pairs,
without multiplying the total cap by the member count. With a zero ordinary
allowance, it remains disabled and all shared storage serves the current pair.
Cross-diagram entries are released after each pair; member caches after the
family task. The combined retained cap is 192 MiB at the default settings.
Zero shared budget or `fca_physical_ef_shared=false` uses ordinary evaluation.
No complete bath-cubed table is retained. These caps exclude live chains,
convolution workspaces and the separately bounded prepared-propagator cache.

The log reports resumed vertices and evictions for ordinary work, and a
separate `FCA E/F shared physical` summary with worker time, convolution
count, cross-diagram hits, and peak retained shared storage. Grouped tasks
include multiple representatives, so task counts cannot be compared directly
with the old one-representative-per-task scheduler.

For isolated numerical/timing comparisons using a saved CIX, Sindex and PPG,
`bench_fca ... cache-resume 5` compares the old ordinary evaluator with the
new cache policy; `bench_fca ... ef-prefix 2` compares E/F sharing with the
already improved ordinary evaluator. These two benchmark modes use DLR
tolerance 1e-10. A `synthetic2` pole table tests computation and timing only,
not whether two poles accurately approximate the physical hybridization.

To recover the unreduced Cartesian reference for A/B comparisons, set
`fca_reduce_cuts=false`, `fca_physical_reuse=false`, and `fca_self_reuse=false`.
All switches default
to the optimized path. `fca_profile=true` prints per-topology sums of worker
seconds, executed tasks, physical convolution counts, prefix hits, and peak
retained cache. Worker seconds are NOT wall times. Phase progress lines
remain wall times. Equivalent per-cut diagnostics in reduced mode are
reconstructed, so their equality is not an independent symmetry test;
use unreduced mode to test those identities.

Self-energy profiling is printed immediately after its checkpoint, before
physical evaluation. It reports convolved functions, convolution batches,
root contractions, unique/reused prepared edges, prefix-budget fallbacks,
and unstreamed task counts separately for A-F.

Further exact kernel optimizations default to:

```json
{
  "fca_root_plans": true,
  "fca_indexed_states": true,
  "fca_convolution_scratch": true
}
```

Root plans lazily cache compatible atomic-block/operator connections and
precompute the root hybridization weights. This is NOT a dense four-index
superoperator and never uses Sindex equality to merge atomic paths. Retained
connection metadata is capped conservatively at 16 MiB per evaluator/rank;
when full, uncached connections are constructed transiently. The evaluator
is private to an MPI worker, reused sequentially across its tasks, and dies
before a subsequent iteration/MC phase. Do not share it between threads.

Indexed bath-state storage uses a direct-address table for at most two
variable bath labels, with a vector of active functions. A streamed label
is fixed and checked separately. No Nbath^3 table is allocated. Sparse atomic
block maps remain unchanged; only the outer bath-state map is replaced.
Iteration/accumulation order can differ, so comparisons allow roundoff.

The convolution scratch workspace retains batch descriptors and matrices
within a task, including the small matrix-product buffers formerly allocated
at every DLR node. Packing/unpacking traverses contiguous columns of the
rank-by-M buffers, and fully assigned output matrices are not zero-filled
first. The existing 32768-scalar-entry chunk limit remains (one oversized
atomic block is still processed alone). This is not a total-RSS bound.
Buffers are caller-owned and never shared through the DLR basis. Existing
TCA calls pass no scratch and keep the original DLR numerical kernel.

Set all three flags false to reproduce the previous optimized FCA path.
The self profile additionally prints `seconds_root_vertex_convolution` and,
with scratch enabled, `convolution_seconds_pack_transform_products_unpack`.
These are summed worker times, not wall times; allocation cost is included
in the relevant section, not measured as a separate allocator-only metric.
The section breakdown instruments the optimized self path. If the prefix
budget forces reference fallback, its full cost is still included in worker
seconds, but not in that section breakdown; check `prefix_fallbacks` first.

No scalar pseudoparticle-path catalog is built. Pole dressing and DLR solves
are shared over baths within each task. Each MPI task is a topology/cut/
direction/pole-triple sector. Rank zero dispatches and workers compute, with
O(1) scheduler storage and only compact output reductions. This deliberately
simple reference scheduler uses size-1 compute workers when size > 1; later
performance work can reuse a compute-capable hybrid dispatcher. The nominal
task counts for each of the self and physical phases are 216*p^3 with cut
reduction, versus 384*p^3 in reference mode; cheap Pauli
zero sectors are still visited. There is no hard-coded two-pole assumption.

## Tests

Root-plan/indexed-state/scratch validation (2026-10-15): 384 one-band sectors
agreed with the previous optimized path within 1.6e-12; full one-band evaluator
sums agreed with the Cartesian reference within 3.8e-14. Independent simplex
and mixed-sign quadrature checks retained errors below 5.3e-9. End-to-end
tests compared 68 analytic files with the previous executable, including
individual switches, the prefix-budget fallback, and 1/2/4 ranks; the worst
relative difference was 1.6e-12. On two ranks PPSigma.0, Sig.out and Gf.out
were identical at printed precision; MPI summation-order differences in
Sig.out were at most 1.1e-11 absolute in the other rank-count checks.

The 48-task four-pole normal-state sample agreed within 1.5e-9 for individual
sectors. The six-task SC self and physical checks used saved SC PPG and a
synthetic signed two-pole hybridization (not a production fit); differences
were below 6.2e-14 and 1.1e-14, respectively, including anomalous F01.
Initial incremental timings were approximately 5-6% for self-energy and
about 2% for the SC physical sample; these are small gains, not a claimed
factor-of-two or complete-solver speedup. Root/state metadata work is not
the dominant cost. The new section timings expose the remaining DLR
transforms and matrix-product cost directly.
With compilation and the other validation jobs finished, the clean 48-task
normal-state repeat took 42.80 versus 40.17 summed worker seconds (1.065x,
about 6% less time). Reference/root/indexed/scratch totals were
42.80/42.78/42.49/40.17 seconds. Thus scratch reuse provides most of this
incremental gain; root plans and outer indexing alone are close to timing
noise in this case. Remaining convolution time was 29.18 seconds: 16.14 in
DLR transforms, 8.01 in matrix products, and 3.61 in packing/unpacking, with
the remainder in batch assembly and related work. This is a sampled-worker
measurement, not a complete-run wall-time prediction.

Optimization validation (2026-10-15): the old executable, unreduced reference,
cut-only, reuse-only and combined paths were compared in isolated one-band
two-pole runs. All 68 analytic result files agreed within 1.1e-12 relative;
PPSigma.0 was identical and physical Sig/Gf differences were below 1.2e-13.
Independent constant-matrix/simplex and mixed-sign quadrature checks passed
on two and four MPI ranks (largest quadrature error 5.3e-9).

For the archived normal-state cluster (rank 39, 8 baths, 236 pseudos, actual
four-pole fit), a 48-task physical A/B sample took 153.98 versus 34.60 summed
worker seconds: 4.45x for chain reuse alone, before cut reduction. A repeat
including explicit physical rotation/pole-label mapping checks gave 4.53x.
Largest G/F relative difference was 5.7e-9; peak retained cache was 60.33 MiB.
A six-task SC block check (rank 37, 12 baths, 474 pseudos, saved PPG but a
synthetic signed two-pole hybridization) gave 2.83x and a maximum difference
of 2.2e-9, including nonzero F01. This is NOT a production SC timing estimate.
All tasks, pole fits, and inputs are fixed between each old/new pair. No full
optimized normal-state/SC solve was used to claim an end-to-end speedup.

The subsequent self-energy optimization was compared on 48 normal-cluster
tasks covering every A-F cut with varied arrows and four-pole triples.
The final sample gave 96.24 versus 60.63 summed worker seconds (1.59x), with
per-topology speedups 1.56, 1.55, 1.60, 1.38, 1.80, 1.76. The maximum relative
self-energy difference was 7.7e-11. This measures the additional self-chain
optimization, NOT the earlier cut reduction or a complete solver wall time.
The prefix/suffix rearrangement provides most of this gain; anchor selection
and prepared-edge sharing did not give a separate measurable speedup in this
sample, but reduce retained data and repeated coefficient preparation.
An exhaustive one-band 384-sector check (all six topologies, eight cuts,
eight arrow assignments, varied signed poles) passed with maximum relative
difference 2.3e-9 for individual sectors. Complete evaluator sums agreed
within 3.7e-14. End-to-end comparisons with the previous executable checked
68 analytic output files for reference, fixed/adaptive anchor, coefficient
sharing, forced prefix fallback, and 1/2/4-rank execution; the worst relative
difference was 1.6e-12. The final independent analytic and mixed-sign
quadrature errors remained 6.3e-11 (self), 1.3e-9 (physical), and 5.3e-9
(quadrature). None of these tests overwrote production run inputs or results.
The updated prepared-edge cache was also checked on six SC physical tasks
with saved SC PPG and synthetic signed two-pole hybridization: nonzero G,
F00 and anomalous F01 agreed with the Cartesian reference within 2.2e-9.
Six SC self tasks (one cut per topology, saved SC PPG, synthetic signed two
poles) took 126.02 versus 95.39 summed worker seconds, 1.32x overall, with
maximum relative difference 2.5e-14. All default plans avoided prefix-budget
fallback; the retained prefix was 1.03 MiB or less except topology E (10.71
MiB). Individual speedups varied from about 1.0x to 6.3x. These bounded
samples do not establish a full-production SC wall-time speedup.

`make bench_fca` builds the bounded A/B profiler. It takes CIX, Sindex,
saved PPG, a pole/residue table (omega followed by independent-channel
residues), and a mode. `sample` checks 48 physical tasks, `quick` checks one
task per topology, and `full` compares complete reference/cut-only/reuse-only/
combined evaluations (use the full mode only for small one-band cases).
Self modes `self` (48 tasks), `selfquick` (six tasks), and `selfall` (384
topology/cut/arrow sectors, intended for one band) now compare the previous
optimized path, root plans, indexed states, and convolution scratch cumulatively.
They also explicitly exercise the zero-prefix-budget fallback. Physical sample
modes compare the previous optimized physical evaluator to the new one;
the Cartesian numerical reference remains available via `full` mode.
`synthetic2` in place of the pole-table path supplies a mixed-sign two-pole
algebra test; it is not a fit to the input hybridization. Sample modes also
test mapped rotational representatives and selected zero-cache fallbacks.
The driver uses beta from the PPG grid, omega_max=30, tolerance=1e-7 and
panel_order=24; it is not an automatic reproduction of arbitrary PARAMS.

Example (no solver iteration, no production results overwritten):

```sh
mpiexec -np 4 ./bench_fca tests/FCA_partial/bctqmc_U_8.0.cix \
  tests/FCA_partial/Sindex.dat tests/FCA_partial/PPG.0 \
  tests/FCA_partial/Delta_FCA_4p.poles.dat sample
```

Build `make test_fca`, then run:

```sh
./test_fca test1/single_band.cix test1/Sindex.dat
mpiexec -np 4 ./test_fca test1/single_band.cix test1/Sindex.dat
mpiexec -np 4 ./test_fca tests/starting_OCA/bctqmc_U_8.0.cix tests/starting_OCA/Sindex.dat partial
mpiexec -np 4 ./test_fca test_SC/starting_OCA/bctqmc_U_8.0_SC1.cix test_SC/starting_OCA/Sindex.dat partial
```

The exhaustive one-band test compares all cuts and G/F estimators with an
independent full-Hilbert-matrix/simplex-volume reference, then tests negative,
zero and positive poles against six-dimensional Gauss-Legendre quadrature.
It does not use the FCA convolution algorithm to generate expected integrals.
The partial cluster tests use synthetic constant pseudo matrices and a common
zero pole to compare the full D/F cut0 and cut2 arrays, including every bath
and orientation. They are algebra checks, NOT convergence tests on the actual
cluster hybridization or propagators.

Initial single-band checks (2026-10-15): analytic self error 6.3e-11,
physical G/F error 1.3e-9, mixed-sign quadrature error 5.3e-9, serial/MPI node
agreement around roundoff. A full 17-pole, eight-rank one-iteration calculation
completed in about 234 seconds (O2 reference build), with Phi4=0.0049398563.
D/F cut-orbit differences were at most 2.2e-11. Lower-order outputs agreed
with the pre-FCA executable to about 2.3e-9 relative or better.

The normal-state and SC synthetic-propagator partial tests also passed:
cut0/cut2 relative differences for D/F were below 1e-15. The final O3 serial
build reproduced the O2 four-rank two-pole run within 4.7e-11 relative in
the fourth-order output values.

### Cross-task dressed-propagator cache

`fca_prepared_cache_mb` defaults to 128 MiB **per MPI rank**, in addition to
the existing prefix budgets and normal solver working storage. Set it to 0
to restore task-local preparation. `fca_reuse_prepared=false` also bypasses
this cache. Neither switch changes pole counts, Sindex grouping, diagram
weights, DLR accuracy, or the set of paths being evaluated.

The cache reuses dressed pseudo-propagator DLR coefficients, their Hilbert
transforms, and lazy atomic-block expansions across pole tasks, cuts, and
the self/physical phases. Keys are exact sampled scalar decorations, with
byte equality checked after hashing: no rounded-frequency identification is
used. It is local to the evaluator and discarded at the end of an FCA call,
so the next self-consistency iteration cannot reuse transforms of stale G.

Transforms still operate on **unique signed-Sindex scalar functions**, not
all atomic matrix entries. Atomic blocks are expanded afterwards using the
signed labels; distinct operator paths are never identified merely because
their scalar propagators agree. The cache retains those expansions too.

Entries reserve conservatively for their complete eventual block expansion
before admission. Least-recently-used unpinned entries are evicted; objects
held by a current chain are not removed from accounting. If an entry cannot
fit, the task uses ordinary uncached preparation. The budget is an accounting
bound on retained cache storage, not a bound on process RSS or task scratch.

`FCA prepared cache` log lines after self and physical evaluation report
cumulative all-rank hits, misses, evictions and bypasses, plus the largest
per-rank reservation. Misses count lookups that failed; bypasses also include
disabled/oversized entries, and thus are not a disjoint request category.
`bench_fca ... cache` checks bounded normal/SC samples against cache-off
evaluation, including self, G, F00 and F01. An optional final sample count
(1..12) limits the work, e.g. `bench_fca ... synthetic2 cache 1` for one SC
self/physical comparison. `test_fca` additionally checks
tiny budgets, eviction, cache disabling, and a fresh pseudo-propagator update.

Cache regression (2026-10-16): one-band mixed-sign/zero-pole cache-on/off
samples agree exactly, including eviction and tiny-budget fallback. A saved
normal-cluster four-pole sample across six topologies also agrees exactly;
it used 14 entries and 20.45 MiB of conservative reservation. The isolated
one-band four-pole full solver matches the previous executable exactly at
the same MPI rank count in 68 diagram files; changing among 1/2/4 ranks gives
at most 3.4e-13 relative differences in those files. These are correctness
checks, not full normal/SC runtime benchmarks.
One bounded SC topology-A self/physical sample with saved PPG and a synthetic
signed two-pole fit also agrees exactly (13.86 MiB reservation). This checks
SC block/sign handling, not convergence in the physical hybridization fit.

The TWO-pole version of that same input gave Phi4=0.00000519487 and roughly
99.7% relative changes in aStval.4/ag00.4. Its hybridization fit error was
about 56%. Thus two poles are a useful plumbing test but NOT an acceptable
approximation for that input. No general accuracy claim is made for SC or
normal-state clusters; pole convergence and MC comparisons remain necessary.

### Experimental streamed shared self-energy prefixes

`fca_shared_prefix` defaults to `false`, preserving the established evaluator
and its task partition. Setting it to `true` enables a compatible subset of
rooted self-energy cuts across A--D. Physical G/F evaluation is unchanged.
The new path requires cut reduction and self reuse; otherwise it uses the
established path. `fca_shared_prefix_mb` defaults to 256 MiB per MPI rank.
Larger family tasks are dispatched before ordinary tasks, allowing the
remaining ordinary work to fill MPI workers as shared tasks finish. The
number of family tasks is `8 * Npole^2`; this reduces available concurrency
as well as total arithmetic, so large-rank scaling must be measured rather
than inferred from a single-worker speedup.

The eligible cuts have three hybridization openings at vertices 1, 2, 3,
and the summed root at vertex 4. With split lines relabeled by their opening
vertex, their decorated chain through vertex 4 is identical. Its matrix
products, convolutions and root contraction are evaluated once and then
branched into the different closing sequences at vertices 5, 6, 7. Diagram
signs, kernel normalization factors and cut multiplicities remain separate.
The implementation checks matching prefix decorations explicitly.

A task fixes two canonical pole labels and the three directions, streams
the first opening's pole and bath, and immediately sums them at that line's
closing vertex in each continuation. Thus subsequent convolutions are not
repeated unnecessarily. No bath-cubed or pole-triple tensor is retained.
`Sindex` still compresses propagator transforms; complete operator histories
are not identified merely because individual propagators have equal Sindex.

Admission uses a block-size estimate. If necessary it excludes the two
members whose first-opening line closes at vertex 5, because their merged
continuations retain two bath labels. Those cuts use the old evaluator;
four later-closing members can still share a prefix. If even that family
does not fit, no shared tasks are scheduled. During evaluation a retained
state byte guard also checks the prefix plus all continuation accumulators.
Exceeding the budget releases these states and recomputes that task using
the existing grouped evaluator. This is an accounting limit, not a hard RSS
limit: matrix-operation temporaries, convolution scratch, and the separate
prepared-propagator cache are additional. The first rejected allocation can
exceed the accounting target before fallback.

`FCA shared prefix` logs the number of admitted members, worker seconds,
tasks, runtime budget fallbacks, prefix/tail convolution function counts,
root contractions, and maximum accounted state storage. Shared worker time
is reported separately, not assigned arbitrarily to individual topologies.
The largest accounted size includes rejected attempts.

`bench_fca CIX Sindex PPG pole-table|synthetic2 prefix [sectors]` compares
one to eight complete canonical direction sectors against the old evaluator
using its original adaptive pole anchor and full pole sum. It compares each
cut separately, not just the summed Phi. This benchmarks only the eligible
family, not the complete solver or physical evaluation. Keep the feature
opt-in until its cost on the intended input and MPI size is established.

Validation (2026-10-16): mixed negative/zero/positive one-band poles and all
eight direction sectors agree to 2.7e-15 in the unit test, including forced
runtime budget fallback. The isolated four-pole one-band solver agrees with
the previous executable in all 68 diagram files to 2.7e-13 over 1/2/4 MPI
ranks, including both six- and four-member admission. Switch-off and tiny
admission-budget runs are identical to the baseline at the same rank count.
The regression also checks final physical outputs and unchanged restart input.

A saved normal-cluster PPG with its actual four-pole hybridization, one full
canonical direction sector, and all six eligible cuts gives a maximum
per-cut relative difference of 4.1e-14. Accounted shared states peak at
85.66 MiB with zero fallbacks. Observed old/new subset times were 514.35/88.30
seconds. This is not a full-solver benchmark: other cuts and physical work
are unchanged, and the reference overlapped a production run that finished
before the shared portion. Treat the timing ratio as indicative, not a clean
hardware-controlled measurement or a prediction of total speedup.

The saved SC PPG with a synthetic signed two-pole hybridization admits four
members. One complete canonical direction sector agrees per cut to 1.9e-14,
with 154.12 MiB accounted states and no fallbacks. Observed subset times were
601.32/95.06 seconds, with the same contention caveat. This tests matrix and
anomalous-channel sign handling, NOT the physical adequacy of a two-pole fit.
No full SC solver run or large-rank scaling benchmark was done for this path.

### Experimental shared physical convolution arcs

`fca_physical_shared` defaults to `false`, independently of the PP
`fca_shared_prefix` switch. When true, six compatible reduced physical cuts
across A--D share a complete left arc and matching short right prefixes.
They have external endpoints 0,4 and three internal lines opening at 1,2,3
and closing at 5,6,7. Relabeling internal lines by opening vertex makes the
left arc identical for fixed canonical poles/directions/baths. External-line
orientation still determines the creation boundary; Bulla insertions are
not moved, and each orbit's sign, normalization and multiplicity are retained.

Each task streams the correlated bath triples and immediately joins the
common left arc with each right arc for G, F00 and F01. Cross-bath labels
are never summed independently between arcs. There is no retained bath-cubed
tensor. Exact prefix plan IDs compare operator-label histories and sampled
decorations; equality of Sindex labels alone never merges matrix histories.
The existing unique-Sindex dressed-propagator cache remains in use.

`fca_physical_shared_mb` defaults to 128 MiB per worker. This bounds cached
short-arc prefixes, excluding the live chains, prepared propagators and
convolution scratch. A full cache evicts least-recently-used entries; oversized
entries are computed without caching. No numerical terms are dropped. Zero selects the previous
evaluator. Disabling physical reuse or cut reduction also keeps the previous
route. Unrelated cuts always use the existing physical evaluator.

Shared families are scheduled first. Their task count is `8*Npole^3`, with
remaining ordinary cuts filling workers as families finish. Early task-count
ETAs are biased high by these larger tasks. `FCA physical sharing` logs the
selected member count; `FCA shared physical` reports family worker seconds,
convolution counts, cache hits/bypasses, complete-left reuses, and retained
cache peak. This family time is separate from ordinary topology profiles.

### Shared boundary contractions (item 2)

`fca_shared_boundary=true` (default) enables a second level of reuse inside
the opt-in `fca_physical_shared` family. Item 1 remains disabled by default:
`fca_extended_branch_transforms=false`. No TCA, PP self-energy, independent
FCA cut, or physical-family scheduling formula changes in this optimization.

For each fixed correlated bath triple, all six compatible members have the
same left arc L but different right arcs R. Instead of recomputing
`Tr[C_b L A_b R]` for every member and external bath b, prepare the transpose
of `C_b L A_b`, and sum it over external baths contributing to the same
physical channel. Joining R then needs only elementwise dot products at DLR
nodes. The ordinary creation matrix and its Bulla replacement produce
separate weights for G, F00, and F01. Signed bath-channel factors are included
in these weights. Each member's pole residues, diagram sign, inverse kernel
normalization, and symmetry factor remain outside. Orbit multiplicities
remain in the caller, as before. Individual family outputs are preserved;
different topologies or correlated internal baths are never merged.

Weights are built lazily only for right atomic blocks actually encountered.
They use exact sparse operator entries for both normal and general Bulla
connectors, avoiding a new dense cubic product. Immutable geometry routing
has a separate 16 MiB cap per evaluator. Numeric weights survive only for the
current left arc and are capped by `fca_boundary_cache_mb` (default 64 MiB per
rank), additional to node-prefix and prepared-propagator caches. The cap is
accounted retained storage, not a total RSS limit. Admission for a whole
member finishes before its output is accumulated: if the cap is insufficient,
that member is evaluated once with the unchanged original boundary routine.
Zero budget or `fca_shared_boundary=false` selects that reference directly.

The shared-physical log reports boundary block builds/hits/fallbacks, peak
numeric storage, and separate prepare/contract/reference worker seconds.
`bench_fca CIX Sindex PPG pole-table|synthetic2 boundary 2` compares the new
and old boundary paths with identical shared arcs, item 1 OFF, alternating
execution order. It compares every family's G/F00/F01 independently. These
bounded samples are not whole-solver or many-rank speedup measurements.

Boundary validation (2026-10-16), two canonical samples per cluster:

| Case | Previous family time | Shared-boundary family time | Family speedup | Boundary-only speedup |
| --- | ---: | ---: | ---: | ---: |
| Normal, saved actual four-pole fit | 8.627 s | 7.204 s | 1.20x | 3.09x |
| SC, signed synthetic two-pole fixture | 90.442 s | 80.320 s | 1.13x | 1.99x |

The boundary-only measurement includes preparing AND consuming the new
weights. Additional numeric storage peaked at 2.805 MiB (normal) and
8.797 MiB (SC); neither sample used the fallback. Every member's G/F00/F01
agreed with the reference to at most 2.1e-14 normwise. The synthetic SC fit
checks matrix/sign algebra, not whether two hybridization poles are a good
physical approximation. Small compiler/validation jobs ran concurrently,
so these numbers should not replace a full application benchmark.

One-band analytic, mixed-pole, quadrature and memory-fallback tests pass.
The isolated four-pole solver regression compares 68 diagram files plus the
final physical outputs at 1, 2 and 4 MPI ranks. At two ranks the default,
disabled, zero-budget, tiny-budget and bounded boundary paths match the
previous executable's numeric diagram tables exactly; the worst difference
across all tested routes/ranks is 3.3e-13. Item 1 is off in these solver and
cluster benchmarks. Production PARAMS and restart/results were not modified.

`bench_fca CIX Sindex PPG pole-table|synthetic2 physical-family [samples]`
checks one to eight bounded canonical pole/direction samples, each summing
all baths and comparing G, F00 and F01 separately for every family member.
Old/new evaluation order alternates between samples. These are subset tests,
not full-solver speedup or hybridization-pole convergence measurements.

Initial validation (2026-10-16): the mixed-sign/zero-pole one-band test gives
maximum G/F00/F01 difference 4.3e-16 over eight direction sectors, explicit
old-route fallback, tiny-cache bypass, and complete pole sums. Two saved
normal-cluster four-pole family samples agree within 4.5e-11 using the test's
1e-12 norm floor (the largest normalized discrepancy is in a tiny contribution).
Both samples exercise LRU eviction. Combined old/shared sample times are
14.11/9.10 seconds, or 1.55x for this subset only, with at most 128 MiB of
retained prefix cache. Complete physical-phase speedup has not been measured.
The isolated four-pole one-band full solver passes comparisons in 68 diagram
files: at most 5.9e-15 at the same rank count and 3.4e-13 across 1/2/4 ranks.
Final physical outputs also pass. Both self and physical sharing are enabled in
these comparisons. Physical switch-off and zero-budget runs match the old
executable exactly at the same rank count. Tiny and eviction-limited caches
also pass. Production inputs/results are not modified by these tests.
Two saved-SC samples with synthetic signed two-pole hybridization also pass,
including nonzero F01. Maximum relative difference is 1.9e-9 in a contribution
of order 1e-11 (absolute difference about 4e-20); larger contributions agree
near 1e-14 relative. Old/shared sample times total 194.78/86.51 seconds, or
2.25x for this family subset, with cache usage below 128 MiB. The synthetic
fit checks contraction/sign correctness, not two-pole physical adequacy.
Neither a complete cluster physical phase nor large-rank scaling has been
benchmarked for this prototype.

## Combining the final physical convolution (2026-10-16)

`fca_physical_final_sum=true` enables an exact linearity optimization within
the existing shared physical family. It requires `fca_physical_shared=true`;
that broader family switch remains opt-in. Set `fca_physical_final_sum=false`
to retain separate final convolutions, with all earlier sharing unchanged.
Neither the pseudo self-energy nor TCA is altered.

The six compatible A-D physical representatives contain one A, two B, two C
and one D. For fixed canonical baths, poles and directions, the last dressed
propagator is independent of the internal closing order. After applying the
last operator separately, combine the two B inputs and the two C inputs:
`a*(H1*G) + b*(H2*G) = (a*H1 + b*H2)*G`. This reduces six final right-arc
convolutions to four, not six complete diagrams to four. Normalization and
orbit multiplicity are included in the weighted input. The actual final-edge
decorations are compared before grouping, and incompatible members use the
separate route. Negative-pole complementary intervals, reflection, signed
baths and the G/F00/F01 boundary operators retain their existing conventions.

`physical_family()` still returns independent, unweighted member results for
diagnostics. The explicit `physical_family_topologies()` entry point returns
six A-F totals with multiplicities already included. Production uses the
latter, avoiding a second application of orbit weights. The accumulation is
streamed per bath tuple and topology: no bath-cubed cache or collection of six
complete right arcs is retained. The combined sparse function can contain
the union of its two members' blocks; its live storage is not part of the
retained prefix-cache cap.

The physical-stage log reports `final_sum=yes/no`; its final shared-family
profile reports `final_sum_groups` and `final_convolutions_saved`. The latter
counts nonempty member convolutions eliminated, not a predicted speedup.
`bench_fca CIX Sindex PPG pole-table|synthetic2 final-sum [samples]` compares
weighted topology totals with the switch off/on, holding previous caches
and boundary sharing fixed and alternating evaluation order.

Initial isolated validation: all eight one-band direction sectors pass with
negative/zero/positive poles, weighted-member reference sums, deliberately
unequal orbit counts, scratch on/off, extended transforms on/off, and
disabled/tiny caches. The combined shared-physical test has maximum normalized
difference 3.2e-14. Two normal-cluster samples using the archived actual
four-pole fit take 7.10 s with separate tails versus 6.00 s with combined tails
(1.18x); maximum normalized G/F00/F01 difference is 1.3e-14. Two SC samples
with synthetic signed two-pole hybridization take 82.36/68.29 s (1.21x).
The largest SC normalized difference is 6.9e-10 in a topology-B contribution
of order 1e-11 (absolute norm difference below 3e-20); larger contributions
agree near 1e-14 relative. These are paired subset measurements, not full-phase
timings or a test of the physical adequacy of a two-pole SC hybridization.

The full four-pole one-band solver regression also passes at 1, 2 and 4 MPI
ranks: 68 diagram tables, plus Gf/Gimpw/Sbulla/Sdyson/Sig/PPG outputs, are
compared against the previous executable. The largest normalized diagram
difference is 3.4e-13. Explicit final-sum disable, zero/tiny/bounded boundary
caches, and the ordinary non-family route are included; the saved PPSigma.0
restart agrees byte-for-byte. Validation files and paired benchmark logs are
kept in `/Users/haule/Documents/ChatGPT/HBCTQMC/fca_final_sum/`. Production
PARAMS, restarts and numerical results are not modified by these tests.

