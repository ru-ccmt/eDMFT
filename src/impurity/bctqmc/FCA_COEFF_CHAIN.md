FCA physical coefficient-space chaining
======================================

PARAMS switch:
    "fca_physical_coeff_chain": true

Default: true. Set false to restore the previous node-based physical
convolutions without disabling prefix reuse, E/F presumming or MPI scheduling.
This option changes only FCA physical G/F evaluation, not FCA pseudo
self-energy, NCA, OCA or TCA. No pole count or DLR tolerance is changed.

Activation and fallback
-----------------------
Requires fca_physical_reuse=true, dlr_tolerance>=1e-7, and a successful
finite B*M=T solve with relative residual <=1e-12. Otherwise the evaluator
uses the previous node path. The FCA evaluator log prints requested and
effective physical_coeff_chain settings, and explains a basis-guard fallback.
The intended production setting dlr_tolerance=1e-7 is eligible.
The cutoff is conservative: complete SC bath sums at 1e-10 failed the
equivalence threshold even though isolated chains and one-band tests passed.
The discrepancies were concentrated in tiny cancelling task contributions
(~1e-14 to 1e-13 norms, absolute changes ~1e-17), not a demonstrated large
error in the complete Green function. The restriction is deliberately cautious.
This is not a guarantee for arbitrary parameters. The residual checks
the precomputed transform, not the accuracy of each diagram. If a new case
shows disagreement, use false and compare fixed-input outputs. There is no
automatic dual evaluation of every diagram or general roundoff detector.
DLR tolerance controls basis compression, not a uniform relative-error
bound on every integrated diagram. In particular, agreement between these
two implementations does not establish the absolute accuracy of the DLR
approximation, the few-pole hybridization fit, or the diagram-order cutoff.

Algebra and representations
---------------------------
Let B map DLR coefficients to interpolation-node values and T denote the
equal-frequency convolution term. Existing convolution evaluates
    values = beta*(T*D+B*O),
where D_l=f_l*g_l and O_l=(H*f)_l*g_l+f_l*(H*g)_l. Atomic matrix order is
preserved. Solve B*M=T once, without forming an inverse. Intermediate
convolutions can then return
    coefficients = beta*(M*D+O).
Constant F operators and closed-bath sums commute with the DLR transform.
The final convolution returns ordinary or reflected node values through
the existing formula. G/F boundary contractions always receive node values.
Time-dependent pseudo-self-energy root contractions are deliberately not
changed by this implementation.

Both ordinary physical arcs and shared A--D and E/F families use this path.
Their pre-final weighted sums remain in coefficients when enabled. Cached
intermediate prefixes contain coefficients; final prefixes contain node
values. Cross-diagram cache keys distinguish these representations.
Initial dressed propagators are fitted in unique signed-Sindex space and
then expanded into atomic blocks. Existing cache budgets remain in force.
Additional basis storage is only one rank-by-rank real matrix; no new
bath-cubed table is introduced. Older experimental branch transforms are
bypassed on coefficient arcs, since their extra node-to-coefficient solve
would misinterpret the data and is unnecessary here.

Testing
-------
test_fca defaults to dlr_tolerance=1e-10; FCA_TEST_DLR_TOLERANCE overrides it
for isolated tests. It compares both paths across all six topologies,
arrow sectors, negative/zero/positive poles, scratch modes, map/indexed
states, tiny caches and presummed physical families. A separate analytic
matrix-exponential/simplex oracle checks six successive convolutions when
the basis guard enables chaining (run with FCA_TEST_DLR_TOLERANCE=1e-7).

bench_fca accepts mode "coeff-chain" and an optional sample count. It
compares complete bath sums for selected ordinary and shared-family tasks
using saved CIX/Sindex/PPG inputs. FCA_BENCH_DLR_TOLERANCE selects the test
tolerance (default for this mode: 1e-7). These benchmarks do not run an
impurity iteration or write production results.
