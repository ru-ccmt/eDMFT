# AAA refinement and balanced fitting

## What the modes mean

The default is now `"balanced"`. An explicit setting in PARAMS takes
precedence; existing `"auto"` settings are not rewritten.

`aaa_fit_mode="joint"` fits all independent hybridization channels with one
shared pole set, minimizing an absolute matrix residual. Large channels
therefore tend to dominate the optimization.

`aaa_fit_mode="balanced"` first computes each channel's RMS amplitude

    s_b = sqrt(sum_n |Delta_b(iw_n)|^2 / N).

AAA and the real-pole refinement then operate on Delta_b/s_b. The scale is
bounded below by 1e-12 times the largest channel scale, to avoid unbounded
amplification of essentially absent channels. This approximately minimizes
the sum of squared relative channel errors. After selecting the shared real
poles, residues are fitted to the original, unscaled hybridization.

This changes the relative weighting, not the requested pole count or the
physical amplitudes. It can protect a small anomalous channel, but does not
guarantee that every channel improves. Fitting is nonlinear, and the number
of retained real poles can also differ. The new log reports each selected
channel's relative RMS error in every mode.

`aaa_fit_mode="auto"` first fits joint. If its relative RMS error exceeds
`aaa_channelwise_trigger`, it also fits balanced. It selects balanced if
the global error improves OR the worst significant channel error improves
by at least 5%. Otherwise it keeps joint. This comparison costs two full
fits when triggered; the explicit balanced mode only performs one.

## Exact computational optimizations

- Pack only exactly nonzero matrix elements during real-pole refinement and
  residue fitting. Inspect the complete positive/negative-frequency sample
  set so off-diagonal Hermitian partners are retained. Tiny nonzero elements
  are never pruned. UTe2 uses seven active entries instead of 49 RHS columns.
- Compute each pole derivative by contracting the frequency index before
  multiplying by residues. Do not materialize a full frequency-by-channel
  derivative matrix for every pole.
- Keep column-pivoted QR for the least-squares solve. Do not form normal
  equations or relax fit tolerances.

The objective is unchanged. Summation order changes at roundoff level, so
L-BFGS-B can follow a different path; pole-by-pole identity is not an
appropriate regression requirement. Test fitted functions and observables.

## Diagnostics

`AAA candidate start` identifies joint, balanced, or channelwise fitting.
`AAA support construction done` reports that phase's wall time. Refinement
reports active/full element counts and progress approximately every 30 s
at objective-call boundaries. `evaluations` includes line-search trials,
whereas the final `iterations` counts accepted optimizer iterations.
The final selected-channel report includes total fitting wall time.

Only the master MPI rank fits. Other ranks receive the compact result by
broadcast; no additional collective communication is introduced.

## Reproducible validation

`bench_refinement.cc` checks the optimized objective and gradient against
the original full-matrix formulas, finite differences, zero inputs, tiny
channels and complex off-diagonal inputs. Build with
`EDMFT_AAA_REFERENCE_HEADER="original/aaa.h"` for the untouched reference.
The benchmark fits both joint and balanced modes on each supplied input.
`check_solver.py` runs old and new executables in isolated restart copies;
the `_balanced` suffix tests explicit balanced against original auto.

No user PARAMS, restart files, or remote running jobs are altered by these
tests. Test results and timings are recorded separately in VALIDATION.md.
