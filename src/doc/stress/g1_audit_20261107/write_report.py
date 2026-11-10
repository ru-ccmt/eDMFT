from pathlib import Path
import json,hashlib,re
import numpy as np
root=Path(__file__).resolve().parent
r=json.loads((root/'material_comparison.json').read_text())
assert len(r)==4
for item in r.values():
 assert item['charge_bytes_identical'] and item['ordinary_scf2_records_identical']
 assert np.isfinite(item['corrected']['tensor_gpa']).all()
 assert abs(np.trace(item['angular_delta_gpa']))<2e-7
lines=['G1 DERIVATION, REPAIR, AND MATERIAL TESTS','7 November 2026','',
'Implemented four source-file edits (exact before/after snapshots are retained):',
'  WIEN2k SRC_lapw2/val_kinstrmt.F and kin_str_integrals.F',
'  eDMFT src/dmft2/stress_kernels.f90 and stress_integrals.f90','',
'1. Restore the product-rule cross coefficient 2 R1 - R0. The previous',
'   g1_deriv_scale=0.5 change was ours and was incorrect.',
'2. Remove independent extra factors of two in the F3 LO--LO R1 and R0 sums.',
'3. Remove undefined IsEven from G1 radial differentiation. Differentiate the',
'   existing cubic spline on the positive radial mesh with its one-sided',
'   endpoint rule. Generic DERGL and density/force radial routines are unchanged.','',
'Both angular repairs multiply A3=sym(n_alpha D_beta), whose trace is zero.',
'They cannot by themselves remove a hydrostatic pressure discrepancy.',
'This is NOT the missing physical variation of radial functions under strain.','',
'Unit validation:',
'  Actual lapw2 and dmft2 G1 loops versus independent Cartesian Hessians: pass.',
'  u, udot, two LOs, complex coefficients, and signed complex density: pass.',
'  Both interstitial implementations versus independent Fourier matrix: pass',
'    (maximum error 4.5e-15).',
'  Positive-mesh radial helper: constants through cubics pass; weighted',
'    second-derivative error <1.5e-12. Cusp and fractional powers converge.',
'  Portable patches reproduce exact target source hashes.','',
'Material tests:',
'  Independent old / radial-only / corrected runs at IDENTICAL fixed vectors',
'  and potentials. DFT was completed before DMFT. Original material directories',
'  and their histories were not changed. WIEN2k real/complex DFT binaries;',
'  four-rank Homebrew MPI for DMFT, with all initial EF/normalization files',
'  restored before each mode.',
'  FeSe is RELA, RKMAX=8; SrVO3 is NREL, RKMAX=7, at both levels.',
'  Central vector headers: FeSe 140 k points x 66 bands; SrVO3 75 x 43.',
'  DFT retains existing TRACE G2 default; DMFT retains its existing FULL default.',
'  No G2 default or physical radial-response convention was changed.','',
'Comparison convention:',
'  mean = Tr(sigma)/3; finite difference = [K(+h)-K(-h)]/(2h*3*Omega0).',
'  h=0.001; +1/Omega energy derivative convention, matching WIEN output.',
'  The entries below are mean stresses in GPa. DO NOT divide them by three again.',
'  Thermodynamic compressive pressure has the opposite sign.','',
'  FD values are reused from the existing charge-converged isotropic triplets.',
'  DFT: prior isolated reconverged WIEN2k-24 triplets (lapw2_audit_20261107).',
'  DMFT: fresh READ-ONLY frozen_sigma_klein.py evaluation of prepared triplets.',
'  No strained SCF or impurity solve was rerun in this G1 audit. The changes',
'  affect stress only; unchanged charge and functional records confirm that',
'  reconverging endpoints for these particular code changes is unnecessary.',
'  These are comparisons at the existing cutoff, not newly certified limits',
'  in basis size, radial response, or impurity self-consistency.','']
for name,x in r.items():
 lines += [name]
 for mode in ['before','radial_only','corrected']:
  lines += [f"  {mode:13s} diagonal {np.diag(x[mode]['tensor_gpa']).tolist()}",f"                mean {x[mode]['mean_gpa']:.12f}"]
 lines += [f"  numerical mean {x['numerical_mean_gpa']:.12f}",f"  analytic - numerical {x['difference_gpa']:.12f}",f"  radial mean change {np.trace(x['radial_delta_gpa'])/3:.12f}",f"  angular mean change {np.trace(x['angular_delta_gpa'])/3:.12f}",'  Charge bytes and ordinary SCF2 records (including forces and functional sums): IDENTICAL.','']
lines += ['DFT mixer verification: fresh merge_scf.py and x mixer -f CASE on the',
'isolated corrected cases reproduce reconstructed tensors within printed rounding.',
'DMFT total tensors use matching saved STR_DENS/STR_CORE and fresh STR_VAL;',
'no user mixer or charge output was overwritten. Existing IEEE invalid/zero',
'flags at MPI shutdown occur with both before and corrected executables; no',
'nonfinite checked stress/SCF output was found. These flags are not claimed fixed.','',
'Installation:',
'  Sources and binaries backed up as *.before_g1_audit_20261107.',
'  lapw2/lapw2c installed in /Users/haule/wien2k_24/.',
'  dmft2 installed in /Users/haule/dbin/.',
'  Build used source ~/activate; Homebrew compiler/MPI paths were explicitly',
'  selected because conda mpif90 refers to an unavailable compiler.',
'  WIEN_STRESS_G1_PREAUDIT=1 is DEFAULT OFF and reproduces only the two old',
'  angular defects while keeping the defined radial helper. It is a diagnostic,',
'  not a recommended production mode. The older G1_LEGACY selector is superseded.','',
'Portable patches are lapw2-g1-audit-20261107.patch and dmft2-g1-audit-20261107.patch.',
'They record ONLY these changes, on top of the pre-existing work. The previous',
'G2 lapw2 patch is in ../lapw2_audit_20261107/lapw2-g2-audit.patch. These audit',
'patches must accompany an older cumulative WIEN patch distribution; this task',
'does not claim to have regenerated that entire distribution.','',
'What remains:',
'  The hydrostatic DMFT residual is essentially unchanged. No empirical term',
'  was fitted to make it disappear. The corrected anisotropy still needs its',
'  own controlled component-strain finite differences of the full functional.',
'  Next: a consistent RKMAX series (first DFT, then DMFT), with fixed G sets',
'  within each triplet, then an analytic radial-response audit if the residual',
'  persists. Increasing RKMAX is a hypothesis to test, not an established cause.',
'  No additional standalone M_Sigma term or pseudo-inverse change was introduced.']
(root/'REPORT.txt').write_text('\n'.join(lines)+'\n')
tex=[r'\section{Implementation and material results}',
'All three repairs described above were implemented in both programs. Real and complex lapw2 and MPI dmft2 were rebuilt and installed, retaining source and executable backups. The diagnostic \\texttt{WIEN\\_STRESS\\_G1\\_PREAUDIT=1} reproduces the two previous angular defects only; it is off by default. No original material directory was overwritten.',
r'\subsection{Mean stress comparison}',
'The following are mean stresses in GPa, already divided by three. The numerical values use existing charge-converged isotropic triplets: the prior isolated DFT tests and a fresh read of the prepared frozen-self-energy DMFT records. New strained SCF calculations were not required for this stress-only change; density and functional outputs are unchanged. These comparisons retain the existing convergence limitations.',
r'\begin{center}\small\begin{tabular}{lrrrr}\toprule',r'Case & Before & Corrected & Numerical & Difference\\\midrule']
for name in ['DFT_SrVO3','DFT_FeSe','DMFT_SrVO3','DMFT_FeSe']:
 x=r[name];tex.append(f"{name.replace('_',' ')} & {x['before']['mean_gpa']:.6f} & {x['corrected']['mean_gpa']:.6f} & {x['numerical_mean_gpa']:.6f} & {x['difference_gpa']:+.6f}"+r'\\')
tex += [r'\bottomrule\end{tabular}\end{center}',
'SrVO3 uses NREL with RKMAX 7; FeSe uses RELA with RKMAX 8, in both DFT and DMFT. The central vectors contain 75 k points and 43 bands for SrVO3, and 140 k points and 66 bands for FeSe. The earlier DFT G2 TRACE default and DMFT FULL default were preserved.',
r'\subsection{Tensor components and isolated effects}',
r'\begin{center}\small\begin{tabular}{lrrr}\toprule',r'Case & Corrected $\sigma_{xx}=\sigma_{yy}$ & Corrected $\sigma_{zz}$ & Radial change in mean\\\midrule']
for name in ['DFT_SrVO3','DFT_FeSe','DMFT_SrVO3','DMFT_FeSe']:
 x=r[name];d=np.diag(x['corrected']['tensor_gpa']);dr=np.trace(x['radial_delta_gpa'])/3
 tex.append(f"{name.replace('_',' ')} & {d[0]:.6f} & {d[2]:.6f} & {dr:+.6f}"+r'\\')
tex += [r'\bottomrule\end{tabular}\end{center}',
'The old/radial-only/corrected comparisons separate the angular repair from the radial derivative policy. The angular change has zero mean within printed precision. It substantially changes the FeSe anisotropy, so this repair matters even though it leaves the hydrostatic mismatch. The remaining anisotropic total stress has not been certified by a new component-strain finite difference in this audit.',
r'\subsection{Regression evidence and reproducibility}',
'Charge files are byte-identical before and after in all four tests. Ordinary SCF2 records, including available forces, chemical potential, band sums, and DMFT functional sums, are identical. Fresh DFT SCF merging and mixer runs reproduce the reconstructed corrected stresses within output rounding. Existing MPI floating-point shutdown flags occur in both old and new dmft2 runs; no nonfinite checked output was found, and those flags are not claimed fixed.',
'The accompanying \\texttt{g1\\_audit\\_20261107} directory contains before/after source snapshots, portable patches, unit checks, material inputs/outputs, build logs, and hashes. The changes were made on top of the pre-existing working tree, without reverting unrelated work. The current result is a corrected kinetic block, not a solution of the remaining hydrostatic error.']
(root/'g1_results_addendum.tex').write_text('\n\n'.join(tex)+'\n')
print('\n'.join(lines[-28:]))
