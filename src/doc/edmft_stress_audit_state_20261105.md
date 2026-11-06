# Analytic stress in DFT+eDMFT: audit state and correction algorithm

Date: 2026-11-05

This note records the current state of the stress investigation.  It separates
facts established by controlled tests from working hypotheses and proposes a
reference algorithm that can be made exact before it is optimized.

## 1. Executive conclusion

The explicit derivative of the normalized DMFT projector is not the dominant
unresolved error.  In its controlled, fixed-coefficient test it agrees with a
joint finite difference.  The omitted radial projector response is real, but
it changes the hydrostatic mean by only -0.006 GPa in MnO and about +0.50 GPa
in SrVO3.

The serious problem is the one-particle LAPW response, especially the way the
generalized-overlap and moving-window pieces are reconstructed from Blaha's
rearranged stress kernels.  Direct finite differences of the complete LAPW
Hamiltonian and overlap matrices are well conditioned and satisfy all
algebraic checks.  The present isolated analytic `M_Sigma` overlap channel does
not agree with this direct metric derivative by a material-independent sign or
factor.  In MnO, the separately evaluated window and normalization channels
are +62.25 and -61.19 GPa and cancel to +1.05 GPa.  That cancellation is too
fragile to use as the foundation of a production stress tensor.

The recommended correction is therefore:

1. retain a DFT-like LAPW stress contracted with the interacting density
   matrix;
2. evaluate the genuinely new energy-moment term with the **complete**
   generalized-overlap derivative;
3. replace the separate projector, window, and normalization additions by one
   derivative of the complete normalized projector;
4. first implement these two objects by a small internal central finite
   difference at the converged central potential; only after closure is
   achieved should they be replaced by analytic matrix derivatives.

No two-particle vertex and no self-energy derivative is required.  Their
absence follows from stationarity of the Klein functional, not from an
approximation.

**Update after the joint endpoint test.**  The complete-metric
(M_\Sigma O_q) replacement is stable.  A joint endpoint projector derivative
is not generally stable unless the cross-strain LAPW metric is known.  Direct
Procrustes alignment in projector space works smoothly for SrVO3, is only
marginal for MnO, and fails decisively for FeSe.  It is therefore retained as
an opt-in reference diagnostic, not as the production implementation of
Eq. (10).  The production form should use the central Sternheimer derivative
in Eq. (16), with degenerate spectral-projector blocks.

## 2. Stationary equations

Let the LAPW coefficient basis at k be nonorthogonal,

\[
H_k A_k=O_k A_k\epsilon_k,\qquad A_k^\dagger O_k A_k=1.
\]

In the stored-band basis define the frequency-summed density and energy-moment
matrices

\[
 \rho_k=T\sum_n G_k(i\omega_n)e^{i\omega_n0^+},
 \qquad
 \rho_{E,k}=T\sum_n(i\omega_n+\mu)G_k(i\omega_n)e^{i\omega_n0^+}.
\]

For a strain coordinate \(q\), use the energy-derivative convention verified
below against finite differences of LAPW1 eigenvalues:

\[
 D_q^{1p}=\sum_k\operatorname{Re}\operatorname{Tr}
 \left[\rho_k H_{k,q}-\rho_{E,k}O_{k,q}\right].                 \tag{1}
\]

Equation (1) is the nonorthogonal-basis form of the stationary one-particle
derivative.  In the noninteracting diagonal limit it reduces to

\[
 D_q^{1p}=\sum_{kn} f_{kn}\epsilon_{kn,q}.                     \tag{2}
\]

The final sign used in a WIEN2k `:STRESS` record must be calibrated once,
because pressure, Cauchy stress, and energy-derivative conventions differ by a
minus sign in common implementations.  This convention issue cannot alter any
of the operator identities tested below.

The energy moment has the exact decomposition

\[
 \rho_{E,k}=\frac12\{\epsilon_k,\rho_k\}+M_{\Sigma,k}.          \tag{3}
\]

Consequently

\[
 D_q^{1p}=D_q^{\mathrm{DFT-like}}
 -\sum_k\operatorname{Re}\operatorname{Tr}
       [M_{\Sigma,k}O_{k,q}],                                  \tag{4}
\]

with

\[
 D_q^{\mathrm{DFT-like}}=
 \sum_k\operatorname{Re}\operatorname{Tr}\left[
 \rho_kH_{k,q}-\frac12\{\epsilon_k,\rho_k\}O_{k,q}\right].   \tag{5}
\]

Equation (4) is the cleanest minimal port to the present implementation: the
first term has the algebraic structure of DFT stress but uses the interacting
density matrix, while the second term is the unique new energy-moment overlap
term.  It must use the **full LAPW metric derivative**, not only the overlap
piece that can be extracted from one muffin-tin correction kernel.

Let \(P_k\) be the raw local-orbital projection in the fixed band-index window,

\[
 S=\sum_k w_k P_k^\dagger P_k,\qquad X=S^{-1/2},\qquad U_k=P_kX. \tag{6}
\]

The complete normalized derivative is

\[
 U_{k,q}=P_{k,q}X+P_kX_q,                                      \tag{7}
\]

\[
 S_q=\sum_k w_k(P_{k,q}^\dagger P_k+P_k^\dagger P_{k,q}),      \tag{8}
\]

and, with \(Y=S^{1/2}\), \(X_q\) is the Hermitian solution of

\[
 YX_q+X_qY=-XS_qX.                                             \tag{9}
\]

If \(C_k\) denotes the frequency-summed quantity conjugate to \(U_k\), the
single explicit projector derivative is

\[
 D_q^{\mathrm{proj}}=2\operatorname{Re}\sum_k
       \operatorname{Tr}[C_k^\dagger U_{k,q}].                 \tag{10}
\]

The raw derivative in (7) must contain all of

\[
 P_{k,q}=P_{k,q}^{\mathrm{geometry}}
        +P_{k,q}^{\mathrm{radial}}
        +P_{k,q}^{\mathrm{subspace}}
        +P_{k,q}^{\mathrm{axes}}.                              \tag{11}
\]

The last term is zero only for fixed local axes.  A later implementation for a
self-energy truncated to diagonal local components should use co-moving local
axes and include their derivative.

At a stationary point, variations proportional to \(G_{,q}\),
\(\Sigma_{,q}\), and the impurity two-particle vertex cancel.  Equations
(1)-(11) contain explicit basis, metric, and projector derivatives only.

## 3. Controlled total-stress comparisons

All entries below are hydrostatic mean stresses in GPa.  The numerical value
uses a symmetric \(\epsilon=10^{-3}\) derivative.  The DMFT endpoints use the
frozen-Sigma Klein functional and a relaxed DFT charge at fixed Sigma.

| material | level | analytic | numerical | analytic - numerical |
|---|---:|---:|---:|---:|
| FeSe | DFT | 12.002640 | 12.138883 | -0.136243 |
| SrVO3 | DFT | 8.514917 | 10.115045 | -1.600128 |
| MnO | DFT | 35.081380 | 34.011832 | +1.069548 |
| FeSe | frozen-Sigma DMFT | -2.642237 | -1.049104 | -1.593133 |
| SrVO3 | frozen-Sigma DMFT | 8.176837 | 4.617396 | +3.559440 |
| MnO | frozen-Sigma DMFT | -92.752420 | -72.347433 | -20.404987 |

After subtracting the DFT baseline, the residual specifically associated with
the current DMFT correction is -1.456890 GPa for FeSe, +5.159568 GPa for
SrVO3, and -21.474536 GPa for MnO.  MnO is therefore the decisive failure.
The DFT baseline is also not yet sufficiently accurate in SrVO3 and MnO and
must remain a separate work item.

The controlled FeSe total comparison is RELA.  The operator-level FeSe test in
the next section is an NREL diagnostic copy; its numbers must not be combined
with the RELA total-stress row as if they came from one run.

## 4. Direct matrix diagnostic

Default-off LAPW1 diagnostics write complete \(H_k\) and \(O_k\) matrices for
the same ordered APW+LO coefficient basis at \(-h,0,+h\).  dmft2 writes the
central \(A_k,\rho_k,\rho_{E,k}\).  The checker then evaluates (1) without any
analytic stress kernel.

### 4.1 Algebraic checks

| material | k points | max \(|A^\dagger OA-1|\) | max \(|A^\dagger HA-\epsilon|\) Ry |
|---|---:|---:|---:|
| MnO | 72 | 6.72e-15 | 2.23e-14 |
| SrVO3 | 75 | 4.17e-14 | 3.28e-14 |
| FeSe NREL | 140 | 1.89e-14 | 1.15e-13 |

Thus the matrix convention, complex conjugation, local-orbital ordering, and
central eigenvectors are internally consistent.

### 4.2 Complete stationary pair

| material | Tr(rho Hq) Ry | Tr(rhoE Oq) Ry | pair, Eq. (1), Ry | one-sided pairs Ry |
|---|---:|---:|---:|---:|
| MnO | 102.110534 | 99.504610 | +2.605924 | 2.582594, 2.629253 |
| SrVO3 | 160.857881 | 162.059509 | -1.201628 | -1.811486, -0.591769 |
| FeSe NREL | 387.119364 | 382.822155 | +4.297209 | 4.219982, 4.374436 |

The large individual terms show why an isolated overlap contribution is not a
safe quantity to modify without a stationary paired formulation.  SrVO3's
pair is especially sensitive to cancellation, so the reference calculation
must use two strain steps and Richardson/step convergence.

### 4.3 Independent eigenvalue check

With ordinary DFT occupations, Eq. (1) must equal the weighted central
difference of the endpoint LAPW1 eigenvalues.  It does:

| material | matrix Eq. (1), Ry | endpoint eigenvalues, Ry | difference Ry |
|---|---:|---:|---:|
| MnO | 0.846619306 | 0.846699241 | -7.99e-5 |
| FeSe NREL | 3.841594824 | 3.841881359 | -2.87e-4 |

This is decisive evidence that the direct \(H/O\) derivative is physical and
is not an artifact of reciprocal-vector labeling or eigenvector gauge.

### 4.4 Energy-moment overlap

The checker reports the contribution in the Eq. (1) convention,

\[
 D_q^{M}=-\operatorname{Re}\operatorname{Tr}[M_\Sigma O_q].    \tag{12}
\]

| material | direct Eq. (12), hydrostatic mean GPa | current analytic diagnostic GPa |
|---|---:|---:|
| MnO | +112.8499 | +81.1073 |
| SrVO3 | +2.3726 | +0.2424 |
| FeSe NREL | +19.0611 | +9.2112 |

The final WIEN sign must still be applied consistently, but the ratios between
the two columns are material dependent.  Therefore the disagreement cannot be
repaired by one global sign or factor of two.  The likely reason is that the
current `K(0)-K(1)` extraction sees only the overlap content of the selected
rearranged LAPW correction kernel, not the full generalized APW+LO metric
derivative.

## 5. Projector diagnostics

### 5.1 Normalized fixed-coefficient derivative

For MnO, the central analytic contraction of the complete normalized
projector derivative is -0.161414 Ry for hydrostatic strain.  The two one-sided
finite-difference contractions are -0.161691 and -0.161137 Ry.  The maximum
relative endpoint change of the normalized projector is 0.00205.

This test validates the raw geometry derivative, the inverse-square-root
normalization derivative in Eq. (9), and the present explicit-projector sign
for this controlled fixed-coefficient problem.  It does **not** yet validate
the complementary-band derivative of a moving fixed-index band subspace.

### 5.2 Radial projector response

| material | baseline explicit mean GPa | with endpoint radial response GPa | change GPa |
|---|---:|---:|---:|
| MnO | -5.3367 | -5.3423 | -0.0056 |
| SrVO3 | -3.5925 | -3.0953 | +0.4972 |

The radial term must eventually be retained, but it cannot explain the MnO
error.

### 5.3 Window/normalization fragility

In the separated detailed MnO calculation,

| channel | hydrostatic mean GPa |
|---|---:|
| fixed-X complementary-window response | +62.2473 |
| normalization response | -61.1943 |
| sum | +1.0530 |

The fused implementation obtains the sum more safely, but the existence of
two large convention- and gauge-sensitive terms is a warning.  The corrected
formulation should construct one full \(U_q\) and contract it once as in
Eq. (10), rather than treating band-window and normalization responses as
independent stresses.

## 6. Other checks

Passed:

- fixed band indices, k lists, G-vector counts, and RMT in the controlled
  triplets;
- `Tr(rhoE)-Tr(epsilon rho)-Tr(M_Sigma)` to about 1e-14 Ry;
- projector reconstruction to about 1e-14--1e-15;
- MnO Matsubara 100/30 versus 300/100 changes the isolated overlap diagnostic
  by only about 0.0016 GPa;
- overlap and Hamiltonian Hermiticity in the matrix dumps;
- real and complex LAPW1 dump conventions;
- default-off diagnostic builds of dmft2, lapw1, and lapw1c.

Not passed or not yet closed:

- frozen-Sigma total Klein derivative for all three materials;
- the direct full-metric energy-moment term versus the current analytic
  `M_Sigma` kernel;
- a gauge-covariant moving fixed-index band-window derivative;
- DFT residuals below about 0.1 GPa in SrVO3 and MnO;
- RELA operator-level \(H/O\) closure;
- nonhydrostatic normal and shear components;
- magnetic two-spin assembly (MnTe), which should follow nonmagnetic closure.

## 7. Algorithm for a correct reference stress

### Stage A: internal matrix finite difference

This is the recommended next implementation because it is conceptually exact
and provides a target for every later analytic optimization.

For each of the six independent strain coordinates q:

1. At the converged central charge and frozen self-energy, build
   \(H_k(+\eta),O_k(+\eta),B_k(+\eta)\) and
   \(H_k(-\eta),O_k(-\eta),B_k(-\eta)\).
   Here \(B\) is the raw local-orbital-to-LAPW-basis projector before the
   band-window normalization.
2. Hold the central spherical and nonspherical potentials fixed.  Initially
   hold every RMT fixed.  Use identical integer G labels, LO labels,
   linearization energies, and band count.
3. Form

   \[
   H_{k,q}=\frac{H_k(+\eta)-H_k(-\eta)}{2\eta},\qquad
   O_{k,q}=\frac{O_k(+\eta)-O_k(-\eta)}{2\eta}.                 \tag{13}
   \]

4. Contract the central \(A_k,\rho_k,\rho_{E,k}\) using Eq. (1).  No
   perturbed eigenvectors and no Sternheimer solve are needed for this step.
5. For the projector only, solve the endpoint generalized eigenproblems for
   all stored bands, construct the same fixed-index subspace, parallel
   transport degenerate subspaces, construct \(P(\pm\eta),S(\pm\eta),
   X(\pm\eta),U(\pm\eta)\), and form

   \[
   U_{k,q}=\frac{U_k(+\eta)-U_k(-\eta)}{2\eta}.                \tag{14}
   \]

6. Add Eq. (10) once.  Do not add separate window or normalization tensors.
7. Repeat with \(\eta/2\).  Require agreement and optionally use

   \[
   D_q(0)\simeq\frac{4D_q(\eta/2)-D_q(\eta)}{3}.               \tag{15}
   \]

This reference path is more expensive than the desired final analytic code,
but it is run only after convergence.  Hydrostatic debugging requires two
matrix builds, not twelve.  The full tensor can be optimized after closure.

### Stage B: analytic production path

Once Stage A matches the frozen-Sigma Klein finite difference:

1. Differentiate the LAPW1 Hamiltonian and overlap construction directly in
   the common APW+LO coefficient basis.
2. Contract matrices before diagonalization; never diagonalize rho or rhoE
   merely to reuse scalar-weight stress kernels.
3. Compute the complementary-band part of the fixed-index subspace derivative
   from the stored external bands.  In a parallel-transport gauge, for
   nondegenerate bands,

   \[
   A_m^\dagger O A_{n,q}=
   \frac{A_m^\dagger(H_q-\epsilon_nO_q)A_n}
        {\epsilon_n-\epsilon_m},\quad m\ne n.                 \tag{16}
   \]

   Use degenerate spectral-projector blocks rather than individual-band
   denominators.  Convergence is controlled by the number of bands already
   written by LAPW1; request more bands when necessary.
4. Add the raw geometry, radial, local-axis, and subspace pieces in Eq. (11),
   solve Eq. (9), and make the one contraction in Eq. (10).
5. Keep the impurity functional derivative out of dmft2.  In a frozen-Sigma
   numerical validation it belongs to the small Klein-functional
   postprocessor already used for endpoint comparisons.  At self-consistency
   stationarity removes its implicit Sigma derivative.

## 8. Integration strategy

The safest minimal integration is to use Eq. (3):

- evaluate the DFT-like part, Eq. (5), with the interacting rho;
- remove the present isolated `M_Sigma` overlap reconstruction;
- add the complete Eq. (12), with the final WIEN sign calibrated in the
  sigma=0 limit;
- replace `explicit projector + window + normalization` by Eq. (10) using one
  complete \(U_q\).

A more invasive alternative is to replace the entire valence stress by
Eq. (1).  That requires a careful re-audit of the Hartree, XC, external
potential, core, nuclear, and mixer double-count terms.  It should not be the
first production change because the present DFT stress assembly contains
large cancellations across lapw0, lapw2/dmft2, lcore, and mixer.

## 9. Required acceptance gates

No production default should change until all gates pass:

1. `A^dag O A` and `A^dag H A` residuals below 1e-10.
2. DFT Eq. (2) agrees with endpoint eigenvalues for every tested component.
3. Sigma=0 dmft2 reproduces the patched WIEN2k DFT tensor component by
   component.
4. The direct \(M_\Sigma O_q\) contraction agrees with the Stage-A internal
   matrix finite difference.
5. The complete normalized projector contraction agrees with a joint
   endpoint finite difference, including moving bands and radial response.
6. Results converge with strain step and number of external stored bands.
7. Frozen-Sigma Klein closure is better than 0.1 GPa for FeSe, SrVO3, and MnO.
8. NREL and RELA pass separately.
9. Normal strains and shear strains pass before magnetic MnTe is attempted.

## 10. Diagnostic interfaces added during this audit

All are default-off:

- `DFTKS_OVERLAP_DUMP=PREFIX`: write LAPW1 O and H matrices;
- `DFTKS_CANONICAL_G_ORDER=1`: strain-stable reciprocal basis ordering;
- `EDMFT_ONEPARTICLE_DUMP=PREFIX`: write central A, rho, and rhoE;
- `EDMFT_MSIGMA_DUMP=PREFIX`: write the isolated energy moment;
- `EDMFT_PROJECTOR_RIX_MODE=WRITE|READ`: freeze projector radial integrals;
- `EDMFT_PROJECTOR_RIX_DERIVATIVE_FILE=FILE`: supply a radial derivative;
- `EDMFT_PROJECTOR_SNAPSHOT_DEGENERACY_TOL=value`: diagnostic gauge sweep;
- `EDMFT_STRESS_PROJECTOR_SIGN=+1|-1`: explicit sign audit;
- `EDMFT_STRESS_DETAILED_DIAGNOSTICS=1`: separated slow stress channels.

Persistent checkers:

- `src/python/check_oneparticle_stress_fd.py`
- `src/python/check_sigma_overlap_fd.py`
- `src/python/check_projector_joint_fd.py`
- `src/python/assemble_stress_reference_correction.py`

The diagnostic executables used for these tests were kept in
`/private/tmp/edmft_stress_diag_bin`; they were not installed as production
executables.

## 11. Joint replacement test and guarded implementation

The two replacements proposed in Section 8 were tested together.  Define

\[
 \Delta T =
   (T_{M_\Sigma}^{\mathrm{full\ metric}}-T_{M_\Sigma}^{\mathrm{current}})
  +(T_U^{\mathrm{joint}}-T_U^{\mathrm{current}}).              \tag{17}
\]

Here

\[
 T_U^{\mathrm{current}}=
 T_U^{\mathrm{explicit}}+T_U^{\mathrm{window}}+
 T_U^{\mathrm{normalization}}.                                \tag{18}
\]

The complete-metric term is numerically smooth:

| material | full-metric (M_\Sigma) mean | current analytic mean | one-sided asymmetry |
|---|---:|---:|---:|
| MnO | +112.849912 GPa | +81.107330 GPa | 0.00668 |
| SrVO3 | +2.372604 GPa | +0.242435 GPa | 0.00506 |
| FeSe NREL diagnostic | +19.061130 GPa | +9.211247 GPa | 0.00540 |

The endpoint normalized-projector derivative behaves very differently:

| material | joint projector mean | one-sided derivatives (Ry) | relative one-sided asymmetry |
|---|---:|---:|---:|
| MnO | -20.266155 GPa | -0.668703, -0.555631 | 0.1847 |
| SrVO3 | -8.592695 GPa | -0.671112, -0.670174 | 0.00140 |
| FeSe NREL diagnostic | -3.109898 GPa | +61.0217, -61.6927 | 365.8 |

The FeSe result is not a large physical response: it is a failed endpoint
gauge comparison.  Projector-space alignment makes the symmetric difference
look moderate while the center is nowhere near the midpoint of the two
one-sided contractions.  Such a result must be rejected.

For the two tests that pass the deliberately loose 25% reference gate, the
combined replacement gives:

| material | old analytic - numerical | corrected analytic - numerical |
|---|---:|---:|
| MnO | -20.404987 GPa | -4.638258 GPa |
| SrVO3 | +3.559440 GPa | +2.404998 GPa |

Thus Eq. (17) identifies important missing terms, especially in MnO, but does
not close the frozen-Sigma derivative.  It must not become the default stress.

### 11.1 Guarded code path

`assemble_stress_reference_correction.py` reads the central diagnostic
decomposition, the full-metric report, and optionally the joint-projector
report.  It checks the two one-sided derivatives and writes an extensive
3-by-3 tensor.  The present diagnostic is hydrostatic, so its trace is placed
isotropically on the diagonal; it is not a general shear correction.

`dmft2` reads this tensor only when explicitly requested:

```bash
export EDMFT_STRESS_REFERENCE_CORRECTION=/path/to/reference_correction.dat
dmft2 dmft2.def
```

The file is applied after MPI reduction and after the response-mode audit,
but before crystal symmetrization and division by volume.  The output records
the exact file, raw tensor, and mean correction with
`DMFT_STRESS_REFERENCE_CORRECTION_*` tags.  With the environment variable
unset, the code path returns immediately and the ordinary result is unchanged.

An end-to-end SrVO3 test supplied -1.154442753 GPa.  The resulting change in
the three diagonal `STR_VAL` entries was -7.847736e-5 Ry/Bohr3, corresponding
to -1.154442448 GPa after printed precision.  This verifies the units, sign,
volume conversion, and propagation to the mixer input.

### 11.2 Required production replacement

The remaining implementation should avoid endpoint band gauges.  At the
central structure, obtain (H_q) and (O_q), then form the complementary
band response directly from Eq. (16).  Within a degenerate block, transport
the spectral projector as a block rather than assigning derivatives to
individual eigenvectors.  Combine that response with the explicit geometry
and radial derivatives of the raw projector, construct one (S_q), solve one
Sylvester equation for (X_q), and contract the resulting (U_q) once.

That central construction is gauge covariant and is the next candidate for a
production default.  The guarded endpoint correction implemented here is a
regression oracle and an assembly test, not the final algorithm.
