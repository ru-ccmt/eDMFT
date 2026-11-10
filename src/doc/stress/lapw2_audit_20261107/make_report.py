from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parent
r=json.loads((root/'comparison.json').read_text())
g=json.loads((root/'geometry_results.json').read_text())
rows=[]
for case,d in r.items():
    a=d['analytic']
    rows += [f'{case}:',
        f"  Numerical mean stress                 {d['numerical_mean_gpa']:.9f} GPa",
        f"  Numerical dE/dV (rounded volumes)      {d['numerical_volume_slope_gpa']:.9f} GPa"]
    for name,q in a.items():
        rows += [f"  {name:20s} diagonal = {q['diagonal']}",
                 f"  {'':20s} mean = {q['mean']:.9f}, analytic-numeric = {q['difference_from_fd']:.9f} GPa"]
    rows += [f"  k points: {d['nk']}; bands at every k: {d['bands']}",
             f"  Matching k, weights, matrix sizes, band counts and PW sets: {all(d['basis_consistency'].values())}",
             f"  E(-), E(0), E(+) [Ry]: {d['energy_after']}",
             f"  Original/updated lapw2 valence-density bytes identical: {d['clmval_identical']}",
             '  Final charge residuals:']
    rows += ['    '+x for x in d['DIS']]
    rows += ['']

text='''LAPW2 G2 AUDIT CORRECTIONS AND PURE-DFT REGRESSION
7 November 2026

Scope
-----
Four concrete defects shared with the recent dmft2 audit were corrected in
WIEN2k-24 SRC_lapw2/val_corr_stress.F. Both real lapw2 and complex lapw2c
were rebuilt. The existing TRACE production default is retained. FULL was
tested separately; an isotropic finite difference cannot validate its
deviatoric/shear tensor, even when the hydrostatic result agrees.

1. LO--LO radial derivative
   The stored table is J_ij = integral dr r^3 u_i(r) u'_j(r).
   Blaha G2 differentiates the bra, giving J_ji, whereas the ket retains
   its linearization energy E_j:
       2 Re sum_ij A_i* A_j (E_j-epsilon_n) J_ji.
   The use of jr3upu_22(jlo,jlop,l,lp) was changed to
   jr3upu_22(jlop,jlo,lp,l). Both LO and angular indices must be swapped.

2. Nonspherical radial angular contraction
   For the Gaunt ordering used by that block, n_alpha*n_beta requires
   c(alpha,t)*c(beta,u), without conjugation. The old expression reversed
   components having one Cartesian y index. The differently ordered
   spherical contractions were deliberately left unchanged.

3. Angular intermediate cutoff
   grad_Omega Y_l has an l+1 intermediate even at the maximum retained
   radial angular momentum. Removing that intermediate violates
   sum_alpha A_alpha,alpha Y_l = n dot grad_Omega Y_l = 0.
   Only negative intermediate angular momenta are excluded now.

4. Fortran bounds safety
   In both duplicate-channel checks, ih>1 is tested in a separate IF.
   A Fortran .AND. is not guaranteed to avoid evaluating array(ih-1).

Validation
----------
The actual corrected lapw2 Gaunt loops were extracted and compiled at -O0
with -fcheck=all. Independent numerical spherical-harmonic quadrature and
angular finite differences test six channels, including l=8 at the cutoff.
The first-record check is exercised with one-based arrays of length one.
An independent two-radial-channel integral tests the production LO indexing.
Results are in geometry_results.json. Both real and complex production builds
completed with the existing compiler flags. The old Makefile has a parallel
dependency omission for val_rela_kinstrmt -> TestCases; the build was completed
serially after the initial parallel attempt. That unrelated Makefile issue
is not changed by this patch.

Material runs
-------------
Prepared source directories:
  FeSe_again/FeSe_stress/DFT_RELA{,_m001,_p001}
  SrVO3/DFT_{center,m001,p001}
under /Users/haule/projects/materials.

FeSe: RELA, PBE, central RKMAX=8, real executable.
SrVO3: NREL, LDA, central RKMAX=7, complex executable.
Both use tetrahedron integration. RMT, fractional atomic positions, k mesh,
fixed band count and the prepared scaled RKMAX/GMAX values were preserved.
All calculations were done in independent copies, leaving the prepared
directories and their histories untouched. Stress was enabled in copied
in0 and in2[c] files for every endpoint. The chain was WIEN2k-24 lapw0,
lapw1[c], patched lapw2[c], lcore, merge_scf.py, mixer. Every x invocation
included -f CASE. Sources and binaries used are recorded by hash.

Leftover Eorb_imp.dat records were present even in these nominal DFT
directories. They were excluded from the isolated pure-DFT SCF assembly,
along with any stale optional DMFT/semicore records. The old EORB offsets
were +0.474566017926 Ry for FeSe and -0.183923905033 Ry for SrVO3. This
explains the large constant changes in absolute reported energy; it is
not an energy change caused by the stress patch. The impurity files are
archived with the report. Their nearly constant offsets did not explain
the stress discrepancies.

Switching the prepared dftKS results to WIEN2k lapw1 produced a small charge
residual. Further ordinary DFT cycles were therefore performed. Old MSR1
mixing histories were archived and restarted because they initially
suppressed this new residual. Relaxation logs and final SCF files are saved.
The first-cycle comparison is retained separately in initial_comparison.json.

Comparison convention
---------------------
For a(t)=a0(1+t), b(t)=b0(1+t), c(t)=c0(1+t), h=0.001,

  mean_sigma_FD = [E(+h)-E(-h)] / [2 h * 3 Omega0],
  mean_sigma_an = (sigma_xx+sigma_yy+sigma_zz)/3.

Conversion: 1 Ry/Bohr^3 = 14710.5164057 GPa, matching mixer.
These are mean normal stresses in the WIEN :PRESSURE: output convention;
the thermodynamic compressive pressure is minus this value. The factor
three has already been included: discrepancies below must NOT be divided
by three again. Rounded endpoint cell constants make dE/dV differ slightly
from the nominal h expression; both are reported.

Analytic assembly is sigma = (STR_DENS - STR_CORE - STR_VAL), followed by
the same conversion as mixer. This independently reconstructed result is
checked against the newly written mixer stress. Old/corrected lapw2 runs
use identical final central vectors and potentials.

'''+ '\n'.join(rows)+'''
What this establishes, and what it does not
------------------------------------------
The four identified implementation defects are repaired and have independent
algebraic tests. In these two high-symmetry examples their net hydrostatic
effect is negligible; the FeSe diagonal redistribution is appreciable.
Neither a good nor a bad hydrostatic result validates separate tensor
components. FULL and TRACE preserve the same trace by construction.

The residual analytic/numerical mismatch is NOT fixed by this patch.
Its origin is not proved to be radial-function response. The legacy LAPW2
G2 construction still lacks the complete primitive-matrix and matching
machinery newly added to dmft2 (notably general LO/APW matching and
strain-allowed nonspherical channels). This patch repairs the four shared
defects in that legacy construction; it is not a full port of the new
dmft2 G2 implementation. That larger port and separate directional-strain
tests remain relevant before certifying the full stress tensor.

No ordinary force routine, DMFT projector, self-energy, or impurity solver
was modified. Charge files are identical for old/corrected lapw2 on fixed
inputs; mode FULL likewise leaves charge unchanged.

Artifacts and installation
--------------------------
lapw2-g2-audit.patch is a standalone incremental patch against the pre-audit
SRC_lapw2/val_corr_stress.F (baseline copy included). It is not a replacement
for the earlier cumulative WIEN patch bundle or an update of its installer.
To apply elsewhere after checking the corresponding baseline:
  patch --dry-run -p1 -d "$WIENROOT" < lapw2-g2-audit.patch
  patch -p1 -d "$WIENROOT" < lapw2-g2-audit.patch
Then rebuild real and complex lapw2 and install both executables.

Installed sources and executables have .before_g2_audit_20261107 backups;
installation.json records before/after SHA256 hashes. Original material
directories remain unchanged. Exact input/output records are in results/;
full isolated calculation copies are under
/private/tmp/lapw2_stress_audit_20261107/runs/.
'''
(root/'REPORT.txt').write_text(text)
print(text)
