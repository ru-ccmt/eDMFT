Enumeration-free OCA validation, 2026-10-18
=========================================

All solver tests use isolated copies of starting_OCA inputs, Norder=2,
N_min_order=2, no importance pruning, dlr_tolerance=1e-7. No live run was
stopped and no production PARAMS or scientific outputs were changed.

Initial -O2 candidate, two MPI ranks, one fixed-input iteration:
                         aStval.2 relative L2      ag00.2 relative L2
single band                     1.02e-13                 5.78e-14
normal cluster                  1.28e-13                 4.81e-15
SC cluster                      6.97e-13                 5.74e-14

References: ../oca_persistent_batches/{single,normal,sc}_before.
These are the established path evaluator, using the identical restart data.
Also checked aStval.1, ag00.1, PPG.0, PPSigma.0, Sig.out, Gf.out, Gimpw.out,
Sdyson.out, Sbulla.out. SC final Sbulla relative L2 = 1.01e-13. The small
anomalous real column agrees to 1.72e-13 relative. Its near-zero imaginary
column (norm 1.43e-7) differs by at most 3.53e-15 absolute.
SC worst nonzero-column errors: aStval.2 2.41e-12, ag00.2 7.72e-13.

Additional checks:
- MPI disabled with three ranks: full replicated answer agrees.
- Active importance cuts: explicit fallback to paths, no silent loss of cuts.
- Nitt=2: first iteration has physical=0 and no physical output evaluation;
  second has physical=1. Compared matrix/3 ranks with paths/2 ranks:
  aStval.2 7.63e-10, ag00.2 7.66e-11, all checked files <1e-8 relative.
  This includes propagation of small numerical differences through Dyson.

UTe2 bounded benchmark (saved PPG.0, actual CIX/Sindex, DLR rank 49):
- 14 original baths, 19685 unique pseudo functions.
- Recovered the previously validated 23-pole fit's real residues from its
  saved fitted values and pole positions (not a new fit to raw Delta).
- Eight tasks cover different baths, both directions, and positive/negative
  poles from about -39 to +10.4. All outputs finite.
- Self task 0.266--0.401 s; physical task 0.415--0.528 s.
- Peak benchmark process RSS 740950016 bytes (about 707 MiB), INCLUDING
  loading the saved dense PPG and its temporary fitting buffers.
- Logs: ute2_realfit.log; pole table: ute2_poles.dat.

Final -O3 bounded benchmark also covers ALL FOUR root/split direction pairs:
self task 0.230--0.517 s, physical task 0.405--0.513 s; peak RSS 810680320
bytes (about 773 MiB), still including the input loader/fitting workspace.
See ute2_realfit_final.log. No full-size path enumeration is involved.

This is NOT a completed UTe2 solver iteration or measured MPI scaling.
The simple extrapolation is about 12 worker-minutes for 1288 self + 644
physical tasks. Perfect 12-rank scaling would be about one minute; actual
concurrency, memory bandwidth, task balance, and startup costs must be timed
in the user's next full run. Measurements above ran alongside the existing
12-rank job and therefore are not clean machine-exclusive timing comparisons.

The paths backend remains the default and all its caches remain untouched.

Final production -O3 executable, three MPI ranks, all three fixtures passed:
                         aStval.2 relative L2      ag00.2 relative L2
single band                     8.20e-14                 1.21e-16
normal cluster                  2.55e-13                 5.84e-15
SC cluster                      8.00e-13                 5.00e-15
All eleven output-file checks passed; SC Sbulla.out relative L2 1.29e-14.
Normal/SC final comparisons use fresh path-backend reference runs.
SHA256: aebbfe50d4b7ad1b05bf00153f4c0c6d751327d8910005a45f12f6b9d4f4b2d0

Timings are not clean scaling tests (concurrent UTe2 and compilation):
O2/two ranks, normal paths 1.23 s including enumeration vs matrix 3.05 s;
SC paths 43.30 s including enumeration vs matrix 14.28 s.
Final O3/three-rank SC matrix 3.98 s, self work 2.49--2.50 s/rank,
physical work 1.49--1.50 s/rank. Do not interpret the difference between O2
and O3 runs as a measured optimization or MPI speedup: load differed.
