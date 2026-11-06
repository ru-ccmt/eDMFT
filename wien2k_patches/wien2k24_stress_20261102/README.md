# WIEN2k-24 analytic-stress corrections

This directory contains portable, reviewable patches for the WIEN2k-24
analytic stress implementation used by DFT+eDMFT.  The patch base is the
WIEN2k-24 source found in three byte-identical local installations:

- `/Users/haule/w24_mini`
- `/Users/haule/Docker/w24_mini`
- `/Users/haule/software/wien_unpack_src`

The production installation used for development was
`/Users/haule/wien2k_24`.  The active cluster installation found on
`rupc-01` is `/opt/apps/ompi/wien/wien2k_24`.

`BASE_SOURCES.sha256` records the 25 pristine WIEN2k-24 source files used to
generate the cumulative patch.  `TARGET_SOURCES.sha256` records the exact 26
patched files, including the new scalar-relativistic LAPW2 source.  Site
Makefiles are excluded from those checksums so compiler and library settings
remain local; the patch changes only LAPW2's object list.

The canonical installer uses `patches/current/0001-wien2k24-complete-stress.patch`.
It also recognizes and upgrades trees produced by the former required-only or
required+debug patch sets through `0002` and `0003`.  `0004` upgrades the
previous complete target to the RELA operator cache; its exact starting
checksums are retained in `PRE_OPERATOR_TARGET_SOURCES.sha256`.  The older granular
patches remain under `patches/required` and `patches/optional` for historical
review, but are no longer the deployment path.

## Required corrections

### 0001: lapw0 built-in PBE/PBEsol VSIGMA

File: `SRC_lapw0/vxclm2.F`

The built-in PBE and PBEsol branches return before the generic LIBXC block
that supplies the derivatives of the exchange-correlation energy density
with respect to the GGA gradient invariants.  Consequently,
`sig_str_u`, `sig_str_d`, and `sig_str_ud` were zero or undefined and
the GGA contribution to the analytic stress was missing.

The patch evaluates the required PBE-family VSIGMA derivatives using the
same PBE/PBEsol constants and spin/cutoff conventions as the existing
energy and potential routines.

### 0002: lcore integration range and normalization

File: `SRC_lcore/core_corr_stress.f`

There are three independent defects:

1. The radial integrals corresponding to Blaha et al. Eqs. (E4)-(E5)
   must end at the muffin-tin radius, but the code used `nrp+141`.
2. The nonspherical potential is initialized only through `nrp`, so the
   old integration range differentiated and integrated undefined values.
3. Equation (E3) gives both nonspherical terms a factor `1/(4*pi)`.
   The value passed to mixer instead used `2/(4*pi)` for one term and
   `1/sqrt(4*pi)` for the other.

The corrected accumulation also avoids converting double precision values
through default-kind `REAL(...)`.

### 0003: mixer stress sign

File: `SRC_mixer/mixer.F`

`scfana` assembles the intermediate `Tstress`.  The physical convention
used in the paper and by finite differences is

`sigma_ab = (1/Omega) dE/d epsilon_ab = -Tstress_ab`.

The distributed source contained `+Tstress`, although the WIEN2k
`SRC_mixer/difflist` records the intended minus sign.

### 0004: lapw2 Eq. (G1) F3 factor

File: `SRC_lapw2/val_kinstrmt.F`

The first-radial-derivative part of the muffin-tin kinetic contribution
`F3` explicitly contains two identical copies of its angular tensor.  Blaha
et al. Eq. (G1), however, has one overall factor `1/2` multiplying this
term.  The historical code retained the factor on the accompanying
`-u/r^2` part, but effectively canceled it on the derivative part by
duplicating the angular expression without a second `1/2`.

The patch restores the common factor in all four realizations: APW-APW,
LO-APW, APW-LO, and LO-LO.  This is not a fitted correction.  The individual
`F1`, `F2`, `F3`, and `F6` diagnostic tensors close algebraically, and the
change removes essentially the complete traceless error in the exact
fixed-G FeSe finite-difference test.

### 0005: cubic l=8/l=10 coefficient assignment

Files: `SRC_lapw0/finl_elect_str.F`, `SRC_lcore/core_corr_stress.f`

The cubic-harmonic conversion for l=8 and l=10 needs three coefficients
`c1`, `c2`, and `c3`.  Both routines assigned the third coefficient to `c1`,
overwriting the first coefficient and leaving `c3` undefined.  The patch stores
the `m+8` coefficient in `c3`.  This is a definite source error, although it is
numerically tiny in the present SrVO3 and FeSe tests because their relevant
stress expansions do not depend appreciably on l=8/l=10 components.

### 0006: spin-down XC channel and parallel atom rotation indices

Files: `SRC_lapw0/xcpot3.F`, `SRC_lapw0/finl_elect_str.F`,
`SRC_lapw0/stress_gga.F`

The spin-down interstitial LDA stress used `VKVXC(J,1)` instead of its own
`VKVXC(J,2)`.  This contaminates spin-polarized analytic stresses but has no
effect on non-spin-polarized FeSe.  In addition, two parallel atom loops reset
the full-atom `ROTIJ` index to zero even when the MPI rank starts at a later
inequivalent atom type.  The corrected loops pre-advance the index by the
multiplicities assigned to preceding ranks.  Serial results are unchanged.

## Optional diagnostics

Five historical patches under `patches/optional` add machine-readable output
enabled only when `WIEN_STRESS_DEBUG` is nonempty and not `0`:

- `:DBG_STR_DENS:` from lapw0;
- `:DBG_STR_CORE_ATOM:` and `:DBG_STR_CORE_RADIAL:` from lcore;
- `:DBG_STR_VAL:` from lapw2;
- `:DBG_STR_MIX:` from mixer.

Normal production output is unchanged when the variable is unset.

The optional patch `0104-lapw1-fixed-g-debug.patch` adds a
serial-only finite-difference diagnostic to LAPW1.  A normal strained
calculation rebuilds the spherical `|k+G|` cutoff, allowing plane waves to
enter or leave the basis.  This produces a discontinuous energy contribution
which is not part of the analytic derivative at a fixed variational space.
The patch is completely dormant unless one of these marker files exists in
the case directory:

- `lapw1.fixgvec.write` writes the central integer G-vector lists to
  `lapw1.fixgvec`;
- `lapw1.fixgvec.read` reads that reference and reuses its exact ordered lists,
  while recomputing Cartesian `k+G` vectors in the strained reciprocal lattice.

The reference and strained calculations must use the same k-point ordering,
lattice type, symmetry ordering, and local-orbital definitions.  The diagnostic
deliberately aborts for parallel LAPW1 so that a reference cannot silently be
split among processes.

`tests/fixed_g_stress_test.py` prepares these files and compares a converged
finite difference with the analytic `STR_DENS - STR_CORE - STR_VAL` tensor:

```bash
python3 tests/fixed_g_stress_test.py prepare \
    --center /path/to/center --minus /path/to/minus --plus /path/to/plus \
    --case FeSe_orig

# Converge the minus and plus cases serially.  Keep lapw1.fixgvec.read present.

python3 tests/fixed_g_stress_test.py compare \
    --center /path/to/center --minus /path/to/minus --plus /path/to/plus \
    --case FeSe_orig --strain 0.001 --mode ab --tolerance 0.1
```

Allowed modes are `iso`, `ab`, `x`, `y`, and `z`.  The comparison refuses to
call the basis fixed unless both strained reference files have the central
SHA-256 hash and all LAPW1 matrix-size records agree.

`tests/check_core_l20.py` independently checks the angular coefficients of
the tetragonal nonspherical core correction.  Run lcore with
`WIEN_STRESS_DEBUG=1`, then pass its output and the matching structure:

```bash
python3 tests/check_core_l20.py \
    --output /path/to/case.outputc --struct /path/to/case.struct
```

For the FeSe test, the raw E4/E5 radial integrals reproduce every stored
`I1` and `I2` Cartesian contribution within `3e-17` in stress-times-volume
units.  This rules out the l=2,m=0 angular algebra and its `1/(4*pi)`
normalization as the source of the remaining c-axis error.

## Deploying on rupc-01

First load the rupc WIEN2k environment.  The available module setup found
during testing was:

```bash
source ~/wmodules
```

The `/opt/apps/ompi/wien/wien2k_24` installation is not a complete build
tree.  `SRC_lapw0` contains its sources, but the installed `SRC_lcore` and
`SRC_mixer` directories contain only `W2kutils.h`.  Therefore the complete
multi-program bundle cannot be applied directly to that directory.

Prepare a private WIEN2k-24 build tree from the same source release, copy
this bundle to rupc, and run the read-only compatibility check there:

```bash
./apply_wien2k24_stress_patches.sh --check /path/to/private/wien2k_24
```

The check reports `pristine`, `legacy-required`, `legacy-full`,
`complete-before-operator-cache`, or `applied`.
Apply or upgrade after arranging write permission:

```bash
./apply_wien2k24_stress_patches.sh --apply --install \
    /path/to/private/wien2k_24
```

The script creates a timestamped source backup before changing anything.
It never invokes `sudo` and refuses incompatible source files rather than
forcing a patch.  After testing the privately built executables, an
administrator can install `lapw0`, `lcore`, `mixer`, `lapw1`, `lapw1c`,
`lapw2`, and `lapw2c` into
`/opt/apps/ompi/wien/wien2k_24`.  `lapw2` and `lapw2c` are now required
outputs.  LAPW1's fixed-G support and all detailed diagnostics are included
but remain dormant unless their environment variables or marker files are
used.  Do not overwrite the shared executables
before preserving their current versions and completing a finite-difference
stress check.

## Current validation status

- All five affected source components compile with local gfortran; the
  resulting programs are `lapw0`, `lapw1/lapw1c`, `lcore`, `mixer`, and
  `lapw2/lapw2c`.
- Installed executables match the rebuilt source-tree executables.
- The lcore normalization diagnostic closes at approximately `1e-13`.
- With no marker file, patched LAPW1 energies and vectors are bit-for-bit
  identical to the original executable.  Central write and read modes are
  likewise bit-for-bit identical.
- An exact fixed-G SrVO3 isotropic finite difference gives 10.114602 GPa.  The
  earlier scalar-RKMAX compensation gives 10.114538 GPa, a difference of only
  0.000064 GPa.  The corrected analytic value remains 8.207874 GPa, proving
  that its approximately -1.907 GPa common stress error is genuine rather than
  a plane-wave cutoff crossing.
- For the FeSe LDA, RKMAX=7.3 test at `z_Se=0.24`, the exact fixed-G
  finite difference gives:

| component | numerical finite difference | corrected analytic |
|---|---:|---:|
| sigma_xx | 12.110895 GPa | 7.700360 GPa |
| sigma_zz | 7.982898 GPa | 3.391520 GPa |

The old moving-cutoff result for `sigma_xx` was about 10.27 GPa, so keeping
the integer G list fixed remains essential for a meaningful derivative.  The
residual analytic-minus-numerical errors after the Eq. (G1) correction are
-4.410535 GPa in plane and -4.591378 GPa along z.  Their mean is
-4.470816 GPa; the traceless remainder is only
(+0.060281, +0.060281, -0.120562) GPa.  Thus the former 13 GPa c-axis failure
was the missing F3 factor, while the unresolved error is now overwhelmingly
isotropic.  The remaining xx/zz error difference is only 0.180843 GPa.

A subsequent source audit found another incompleteness: `val_corr_stress.F`
implements only the first three, spherical-potential lines of Eq. (G2), not
its last three nonspherical-potential lines.  The hydrostatic invariant of
those three lines is now independently validated (see below).  A literal use
of the complete tensor still produces an unphysical deviatoric part because
the corresponding fixed-sphere shape response is not yet available.  The
validated trace is therefore the production default; setting
`WIEN_STRESS_G2_MODE=NONE` restores the historical result.  The full
experimental tensor remains excluded from the portable production patch set.

Additional localization tests give the following results:

- Replacing tetrahedron occupations by the full `TEMPS 0.005` Mermin free
  energy changes the FeSe finite differences only modestly.  The exact fixed-G
  values are 12.198775 GPa (`xx`) and 7.740950 GPa (`zz`), versus analytic
  values 7.85424 and 3.05308 GPa.  The occupation/free-energy issue discussed
  around paper Eq. (24) is real but is not the dominant residual here.
- Increasing RKMAX from 7.3 to 8.0 changes the converged center analytic stress
  only to about 8.81 GPa (`xx`) and 3.21 GPa (`zz`); it does not drive the
  analytic result toward the numerical derivative.
- A fresh NREL, PBE, RKMAX=7.3 hydrostatic FeSe check with the complete
  cumulative patch gives `13.497083 GPa` from the exact fixed-integer-G
  finite difference and `10.74188 GPa` analytically
  (`10.58299, 10.58299, 11.05967 GPa`).  The remaining mean-pressure error is
  `-2.75520 GPa` (`-20.4%` relative to the numerical derivative).  Thus the
  previously observed hundreds-of-GPa output was stale stress data, but the
  several-GPa frozen-augmentation/basis residual is still present.
- Increasing the Weinert pseudocharge parameter from NCON=9 to NCON=12 changes
  the center density-stress contribution by only -0.073 GPa (`xx`) and
  +0.149 GPa (`zz`).  Appendix F's frozen-pseudocharge assumption is therefore
  not large enough to explain the 4--5 GPa component residual in this case.
- The raw finite-difference groups cancel at the scale of thousands of GPa:
  for the in-plane Mermin test, approximately +2312.36 GPa (density),
  -1076.67 GPa (valence eigenvalues), and -1223.49 GPa (core eigenvalues).
  These raw derivatives are not one-to-one counterparts of the analytic
  groups, but they show why sub-per-mille inconsistencies remain visible.
- The trace of the analytic electrostatic contribution satisfies the Coulomb
  homogeneity identity `Omega Tr(sigma_es) = -E_es` to about 0.4 GPa in trace.
  Thus a wholesale Madelung or interstitial Coulomb sign/factor error is ruled
  out.  The unresolved trace error is much smaller and is now concentrated in
  the delicate kinetic/Pulay and nonspherical Eq. (G2) balance.

## Uniform-RMT sensitivity diagnostic

The same-geometry radius derivative is **not** a stress correction.  Appendix
D defines the published stress at fixed spherical muffin-tin boundaries and
Eq. (D6) explicitly discards the radial linear response in the
frozen-augmentation approximation.  Changing all radii instead probes a
different basis partition:

```text
dE/deta_R = sum_a R_a dE/dR_a,
S_R = -(1/(3 Omega)) sum_a R_a dE/dR_a.
```

`S_R` is printed in GPa only as a convenient measure of sensitivity.  It must
not be added to any stress component.  `src/python/rmt_stress_correction.py`
(historical filename) interpolates the converged density to
`RMT*(1-delta)` and `RMT*(1+delta)`, converges both branches, and writes
`case.rmtstress` and `case.rmtstress.json`.  Its output is now explicitly
tagged `RMT_SENSITIVITY`; it no longer writes a corrected tensor or pressure.

The original two-cycle shortcut is not generally converged.  In the PBE FeSe
control, the apparent RELA sensitivity moved from `7.7744 GPa` after two cycles
to `9.9251 GPa` after full charge relaxation; NREL similarly moved from
`7.2641 GPa` to `9.8235 GPa`.  The helper therefore uses 20 cycles as the
default *maximum* and requires normal `run_lapw` convergence.  An earlier
successful stop is accepted even though it produces fewer than 20 `:ENE`
records; reaching the maximum without convergence is an error and reports the
last charge distance.  Histories with unequal side lengths are aligned at
their final cycles, and different record counts no longer produce a false
matrix-size mismatch.  Increase `--cycles` when 20 is insufficient and always
inspect the final-aligned sensitivity history.  The helper explicitly requests
both energy and charge convergence, defaulting to `--ec 1e-6 --cc 1e-5`;
these may be changed on the command line.  This avoids `run_lapw`'s usual
energy-only default.  The private plus/minus copies also change inherited
`STR`, `STRF`, or `FOR` selectors to `TOT`: only total energies enter the RMT
derivative, and recomputing the analytic stress in every relaxation cycle is
unnecessary.  Reference-directory inputs are never changed.

The utility also scales RKMAX by the same factor as every RMT.  This preserves
`Kmax=RKMAX/RMT_min` and the integer plane-wave list at the fixed reference
lattice.  It rejects the result if any LAPW matrix-size record differs between
the two sides.  This is essential: leaving RKMAX fixed mixed the radius
derivative with a discontinuous basis-cutoff change in the SrVO3 validation.

For the FeSe Mermin test, the fully converged diagnostic is
`sum R dE/dR = -0.449110 Ry`, or `S_R=+4.16327 GPa`.  The numerical value is
well converged, but its earlier apparent agreement after addition to the
analytic pressure was accidental and must not be used as validation.

The actual fixed-boundary stress error converges strongly with basis size.
For SrVO3 the analytic-minus-numerical hydrostatic residual changes from
`-3.5259 GPa` at RKMAX `6.31` to `+0.2224 GPa` at RKMAX `8.0`.  For MnTe it
changes from `-4.6373 GPa` at RKMAX `7.0` to `-0.8664 GPa` at RKMAX `8.0`.
This identifies frozen augmentation/basis incompleteness, not an omitted
`R dE/dR` chain-rule term, as the dominant several-GPa error.

Typical use is:

```bash
source ~/activate
~/dbin/rmt_stress_correction.py /path/to/reference -f case
```

The script is installed automatically by the existing `src/Makefile` rule
that copies `src/python/*.py` into the eDMFT executable directory.

### LAPW2 stress performance

Three mathematically exact rearrangements remove the main stress-only
bottlenecks: the G1 radial factors are evaluated once per band and multiplied
by their Cartesian tensor; the dense matching-coefficient and interstitial
step-function contractions use BLAS; and the independent RELA band responses
are distributed with OpenMP.  Separate CPU/wall timers for the NREL and RELA
muffin-tin parts are printed at the end of LAPW2.

The OpenMP RELA acceleration is enabled through the normal WIEN2k `.machines`
file, for example:

```text
omp_lapw2: 12
```

On the local RKMAX=8 regression calculations, SrVO3 decreased from about
102 s to 60 s, while RELA MnTe decreased from about 375 s to 78 s.  The BLAS
and interstitial changes preserve the immediately preceding printed `STR_VAL`
tensor to 12 decimal places.  Relative to the original scalar contractions,
the largest internal-component change was `2.6e-9 Ry/Bohr^3` (less than
`4e-5 GPa`), solely from floating-point summation order.

### Cached scalar-relativistic operators (2026-11-05)

`val_rela_kinstrmt.F` now precomputes radial/angular bilinear operator
matrices once for the current atom and radial basis. Band and k-point
work becomes a sparse contraction of LAPW/LO coefficients with those
operators. No stress formula or strain derivation was changed or fitted.
Surface terms, ordinary LOs and the existing relativistic/SDLO treatment
are included; explicit band-energy factors remain outside the cache.

The cache compares the actual radial functions, derivatives and fluxes,
potential derivative, mesh, active channels, dimensions, diagnostic scales,
and quadrature flags before reuse. Perturbing/restoring the radial basis
and changing/restoring the quadrature mode each force a rebuild; the restored
results reproduce the initial cached buffers bit for bit. Set
`WIEN_STRESS_RELA_CACHE_DIAG=1` for build/reuse timing records in
`operator-cache.log`; the default creates no additional diagnostic file.

Measured full LAPW2 wall times on macOS ARM64 with gfortran 15.2, at matched
thread counts (`VECLIB_MAXIMUM_THREADS=1`):

| case | OpenMP threads | before | cached operators | speedup |
|---|---:|---:|---:|---:|
| MnO | 1 | 40.19 s | 4.61 s | 8.7x |
| FeSe | 1 | 456.32 s | 54.98 s | 8.3x |
| MnO | 4 | 14.72 s | 4.06 s | 3.6x |
| FeSe | 4 | 148.68 s | 40.07 s | 3.7x |

These timings include operator construction and all other LAPW2 work.
Two operator builds replace repeated work across 144 MnO or 560 FeSe calls.
Real/complex MnO, diagnostic branches, and OpenMP 1/4 consistency were checked.
MPI and other compilers have not been independently validated.

Moving integration ahead of coefficient contraction changes floating-point
summation order. Every returned complex buffer was compared with the
preceding routine: maximum absolute differences were `2.3e-13` in per-band
buffers and `1.1e-15` in high-precision `:REL_KIN_*` records. Energy,
occupations, densities, printed valence stress and mixer total-stress/energy
records remained byte-identical in these tests. The user explicitly accepted
these disclosed rounding-level differences before production installation;
bitwise equivalence of the RELA buffers is not claimed.

The cumulative pristine and both legacy upgrade patches include this code.
An existing pre-cache complete installation is upgraded by
`0004-upgrade-complete-to-operator-cache.patch`; unrelated source files and
site-specific Makefile settings are preserved.

Local synchronization note: `/Users/haule/wien2k_24/SRC_lapw1/modules.F`
already contained an unrelated `fixed_nbands` declaration before this
optimization. It is deliberately not folded into this LAPW2 patch. The
strict full-tree verifier reports that local LAPW1 divergence; all LAPW2
target checksums match the newly installed operator-cache source.

### Default-off radial-response diagnostic

LAPW2 can dump the actual radial solutions generated from each converged
strained potential.  This is an IBS diagnostic precursor, not a stress term,
and has no cost or output unless explicitly requested.  Run serial LAPW2 in
the minus, center, and plus calculations with distinct filenames, for example:

```bash
WIEN_STRESS_RADIAL_DUMP=minus.radial_fd x lapw2 -f case
WIEN_STRESS_RADIAL_DUMP=center.radial_fd x lapw2 -f case
WIEN_STRESS_RADIAL_DUMP=plus.radial_fd x lapw2 -f case
~/dbin/radial_response_fd.py --minus minus.radial_fd \
  --center center.radial_fd --plus plus.radial_fd --strain 0.001
```

The analyzer reports finite-difference norms and sphere-boundary responses of
`u_l`, `udot_l`, and their radial derivatives.  It deliberately performs no
stress correction.  The next theoretical step is to contract this measured
radial response with `(H-epsilon)` and the central-state coefficients to form
the omitted frozen-augmentation IBS contribution.

## Validated nonspherical G2 trace

The last three lines of Blaha et al. Eq. (G2) contain the response of the
nonspherical muffin-tin potential: the matching-coefficient derivative, the
radial-coordinate derivative, and the angular derivative.  Their full tensor
cannot be enabled safely yet, but their trace has two strong simplifications:

```text
Tr(rhat_alpha rhat_beta) = 1,
Tr(angular generator) = 0.
```

Consequently the radial three-Gaunt recoupling must collapse to the ordinary
one-Gaunt nonspherical matrix element and the angular contribution must have
zero trace.  The implementation evaluates both radial forms independently.
For FeSe they agree to about `3e-19 Ry/Bohr^3`, while the angular trace is zero
at numerical precision.  This is a substantially stronger check than merely
reproducing the final pressure.

TRACE is now selected when the environment variable is absent.  It may also
be requested explicitly with

```bash
export WIEN_STRESS_G2_MODE=TRACE
```

when running `lapw2`.  `NONE` retains the historical result.  `MATCH`,
`RADIAL`, `ANGULAR`, and `FULL` are
diagnostic modes; in particular, `FULL` is not a production option.  Detailed
per-atom checks are printed only when `WIEN_STRESS_G2_DIAG=1` is also set.

In the FeSe Mermin regression the trace mode lowers every diagonal of
`STR_VAL` by `1.9993603e-5 Ry/Bohr^3`.  Since total physical stress is assembled
as density minus core minus valence, this raises every physical diagonal by
`0.29411605 GPa`.  This is a derived, independently checked G2 contribution.
It must be assessed directly against fixed-boundary finite differences; it
must not be combined with the uniform-RMT sensitivity to manufacture a
corrected pressure.  Remaining several-GPa low-cutoff discrepancies are now
treated as frozen-augmentation/basis-convergence errors.

## Scalar-relativistic valence stress (`CALC=RELA`)

The local production tree now contains an explicit scalar-relativistic
valence correction in `SRC_lapw2/val_rela_kinstrmt.F`.  The historical RELA
path used scalar-relativistic radial functions and normalization, but still
differentiated the nonrelativistic operator `p_alpha p_beta`.  That is not the
metric derivative of the implemented Koelling-Harmon equation.  In Rydberg
units the relevant upper-component problem is

```text
H_SR(E) = -div[A(r,E) grad] + V,
A(r,E)  = 1/M(r,E),
M(r,E)  = 1 + (E-V)/(4 c^2).
```

The new routine replaces the muffin-tin gradient tensor by the Hermitian
Koelling-Harmon tensor `Q_SR`, subtracts the historical large-plus-small
tensor `Q_old`, and adds the exact APW flux companion at the sphere.  Two
quantities caused by the energy- and potential-dependent mass are useful in
the derivation:

```text
N^E_ab = (1/(4 c^2)) Re integral (A d_a psi)^* (A d_b psi) d^3r,
rho_KH = |psi|^2 + (1/(4 c^2)) |A grad psi|^2.
```

`N^E_ab` belongs to a direct derivative of the eliminated, energy-dependent
upper-component equation.  It must *not* be added to WIEN2k's implemented
matrix problem: LAPW1 already constructs `H c = epsilon S c` with the stored
large-plus-small `RINT13` overlap, and the existing Pulay terms differentiate
that `S`.  Adding `N^E_ab` would double count an overlap that is not present in
the implemented generalized eigenproblem.  It remains available only as a
diagnostic with default scale zero.

The difference between `rho_KH` and WIEN2k's stored large-plus-small density
multiplies the full tensor `r_b d_a V`.  This mass-potential term is required
and carries an explicit factor one-half in the LAPW2 buffer because the
existing accumulation applies a common factor two appropriate to the two
kinetic metric derivatives, whereas `Delta rho * delta V` contains only one
strain derivative.  That one-half is therefore derived, not fitted.

The RELA correction is applied by default when `CALC=RELA`; it is a no-op for
`CALC=NREL`.  The exact historical RELA result can be recovered with

```bash
export WIEN_STRESS_RELA_MODE=NONE
```

`WIEN_STRESS_RELA_MODE=DIAG` evaluates and prints the new terms without adding
them, while `APPLY` explicitly selects the production default.
`WIEN_STRESS_RELA_ENNAT_SCALE` defaults to zero and is retained only to
reproduce the rejected eliminated-equation overlap.  The production
mass-potential scale `WIEN_STRESS_RELA_MASSV_SCALE` defaults to one.  The
legacy small-component Pulay sector changes the FeSe pressure by only
`0.003 GPa`; the exact Pulay-D and experimental small-surface additions are
also negligible and remain default-off.

LAPW0 and LCORE now write fresh `STR_DENS` and `STR_CORE` records for RELA.
Their density-functional and localized-core correction formulas remain valid;
the old NREL-only output guards merely allowed stale records to survive into
mixer.  Both real and complex LAPW2 variants have been rebuilt.

The fixed-integer-G FeSe RELA-minus-NREL hydrostatic comparison is the cleanest
test because it cancels most of the common, still imperfect NREL augmentation
stress.  It gives:

```text
                                      NREL          RELA        RELA-NREL
fixed-RMT numerical pressure       17.701055     18.579253       0.878198 GPa
raw analytic pressure               9.426945     10.344817       0.917872 GPa
analytic-minus-numerical differential                              0.039674 GPa
uniform-RMT sensitivity             9.823511      9.925111       0.101600 GPa
```

The uniform-RMT line is diagnostic and is not included in the differential
stress error.  The raw differential comparison isolates the new RELA sector
and reduces its former tens-of-GPa error to about `0.04 GPa`.  Absolute
pressures still require basis-convergence tests; a clean independent RELA
system and component-strain tests are also still needed.
The attempted SrVO3 check was rejected because its starting charge was a
DFT+eDMFT state incompatible with a plain-DFT RELA convergence test.

The installed local executables have backups named
`*.pre_rela_stress_20261103` and `*.pre_stress_omp_20261103`.  The cumulative
patch now contains the RELA extension, safe default TRACE, APW/G2 corrections,
LAPW1 matrix support, RELA LAPW0/LCORE output fixes, performance changes, and
the default-off radial-response dump.  Applying only historical patches
0001--0006 is intentionally no longer documented as a complete installation.

## Audit of the published formulas

The derivation in Belbase, Troester, and Blaha, Phys. Rev. B 104, 174113
(2021), is useful but should not be transcribed literally.  Confirmed or highly
probable publication problems are:

- the prose after Eq. (1) prints `r_hat_a=r/r_a`; it must use the local vector
  `r_a/r_a`;
- Eq. (A1) joins the vector derivative and the derivative of its magnitude in
  one equality.  They are separately
  `d V_vec/d eps_ab = +/- (V_a e_b+V_b e_a)/2` and
  `d|V|/d eps_ab = +/- V_a V_b/|V|`;
- Eq. (D5) is missing `R_a` when the prime on the spherical Bessel function is
  its derivative with respect to the dimensionless argument `|V|R_a`;
- Eq. (24) is dimensionally and statistically wrong as printed.  The Mermin
  functional is `F=E_tot-T*S`, not `E_tot-sum(n_vk sigma S)`, and conservation
  of electron number does not make each occupation derivative second order;
- the last term of Eq. (G5), written with
  `V_M(tau_a) Y_1t(tauhat_a)Y_1t'(tauhat_a)`, is origin dependent and undefined
  for an atom at the origin.  It cannot be a general expression for
  `d V_M(tau_a)/d eps_ab`; the WIEN2k source actually uses a different local
  sphere/boundary decomposition;
- the symbol `s` is reused for two unrelated angular-momentum roles in Eq.
  (G2), making a literal implementation ambiguous.

Equations (11), (18)--(23), Appendix C, and the overall radial/angular
structure of Eq. (G2) survive direct checks.  The assumptions in Appendix F
(strain-invariant interstitial and pseudocharge Fourier coefficients) and the
frozen-augmentation treatment are approximations, not general identities.
The paper itself reports that the earlier Klueppelberg electrostatic
derivation failed numerically; its own final Appendix-G electrostatic formula
also contains the origin-dependent misprint above.  Therefore every compact
Appendix-G term must be validated against invariants and finite differences
before it is reused in eDMFT.

These patches correct demonstrable implementation errors; they do not yet
establish that the entire WIEN2k stress implementation is correct.

When running isolated component steps, regenerate the combined SCF record
after lcore and before mixer (for example with `merge_scf.py`).  Otherwise
mixer can read stale stress components from the previously assembled SCF.

## Recommendation for eDMFT distribution

For now, distribute this versioned patch bundle and run the verification
script as an eDMFT preflight check.  This is preferable to immediately
copying complete WIEN2k programs into eDMFT because:

- the remaining sigma_zz discrepancy still needs localization;
- these changes are tied to the WIEN2k-24 source layout;
- distributing small patches minimizes duplication and licensing concerns.

After validation on cubic and tetragonal systems, for LDA and PBE, we can
decide whether to:

1. keep version-specific patches;
2. include only independent corrected stress kernels in eDMFT; or
3. maintain eDMFT versions of lapw0/lcore/mixer stress functionality.

The second option is likely the cleanest long-term direction if additional
upstream stress defects are confirmed.
