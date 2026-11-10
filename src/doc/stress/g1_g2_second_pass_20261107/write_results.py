from pathlib import Path
import json,hashlib
r=Path(__file__).resolve().parent
read=lambda n:json.loads((r/n).read_text())
g=read('g1_results.json');co=read('combined_results.json');mat=read('material_comparison.json');ma=read('matching_results.json');su=read('surface_results.json');ig=read('interstitial_results.json');g2=read('existing/complete_g2_results.json')
rows=[('G1, full complex density; WIEN and dmft2',max(g['production_complex_density_error'].values())),('G2, all six primitive contractions',g2['six_G2_contractions_max_error']),('Combined H/O, noncommuting DMFT',co['DMFT_noncommuting']['max_error']),('Combined H/O, diagonal DFT',co['DFT_diagonal']['max_error']),('Combined H/O, geometric cancellation',co['geometric_cancellation']['max_error']),('Matching, direct endpoints at $h=3\\times10^{-6}$',max(ma['3e-06'].values())),('APW boundary, complex density',max(x['complex_matrix_contraction_error'] for x in su.values())),('Interstitial, complex density',max(x['max_error'] for x in ig.values()))]
s=r'''\subsection{Measured regression results}
The absolute errors below are in the units of each synthetic fixture, not GPa. For the combined test the finite-difference step is $h=2\times10^{-5}$. The maximum is taken over all six tensor components. The $\sim10^{-7}$ combined error includes the production logarithmic-mesh radial quadrature; it is well below the $3\times10^{-5}$ test tolerance.
\begin{center}\begin{tabular}{p{0.73\linewidth}r}\toprule
Test & Maximum error\\\midrule
'''
for label,x in rows:s+=label+f' & ${x:.3e}$'+r'\\'+'\n'
s+=r'''\bottomrule\end{tabular}\end{center}
The existing normalization, window, direct-$C$, and basis-covariance regressions also pass. The positive-mesh G1 derivative reproduces polynomials through degree three to roundoff-weighted accuracy; cusp and fractional-power tests converge with mesh refinement. A deliberately restored pre-audit G1 mode fails the Cartesian-polynomial test, so the suite is sensitive to the previously corrected factors and not just checking that two current routines agree.

\subsection{Fresh SrVO$_3$ and FeSe reruns}
Installed binaries were rerun in isolated copies for both DFT and eDMFT. FeSe uses RELA and RKMAX 8; SrVO$_3$ uses NREL and RKMAX 7. The eDMFT vector data contain 66 bands at each of 140 $k$ points for FeSe, and 43 bands at each of 75 $k$ points for SrVO$_3$. The archived first-pass input provenance and endpoint references were reused unchanged.

The table compares the \textbf{mean stress} $\bar\sigma=\Tr\sigma/3$ in GPa. Thus the difference is already divided by three relative to an isotropic linear-strain trace; do not divide the displayed difference by three again. The finite-difference values are the existing converged endpoint reference, not fresh endpoint reconvergence in this pass.
\begin{center}\begin{tabular}{lrrr}\toprule
Case & Analytic & Existing FD & Analytic $-$ FD\\\midrule
'''
for key,x in mat.items():
 level,case=key.split('_');case=case.replace('SrVO3',r'SrVO$_3$')
 s+=f"{level} {case} & {x['mean_gpa']:.6f} & {x['existing_fd_mean_gpa']:.6f} & {x['analytic_minus_fd_gpa']:+.6f}"+r'\\'+'\n'
s+=r'''\bottomrule\end{tabular}\end{center}
The fresh total tensors exactly reproduce the first-pass results at the precision of the saved SCF component records. Valence charge files are byte-identical, and ordinary nonstress SCF2 records are unchanged in all four cases. DFT merge/mixer was also run on the isolated outputs and agrees with the component assembly to its printed precision. For eDMFT the total is assembled from the saved density, core, and valence stress records using the same convention.

The code therefore remains reproducible after the first-pass corrections. The mismatch has not decreased in this pass, because the independent tests found no additional justified production correction. In particular the DMFT residuals remain about $0.607$ GPa in SrVO$_3$ and $0.560$ GPa in FeSe. Agreement of the synthetic NREL G1/G2 tests does not explain or remove those residuals.
'''
# Turn engineering numeric exponents into readable LaTeX.
import re
s=re.sub(r'\$(\d+\.\d+)e([+-]\d+)\$',lambda m:'$'+m[1]+r'\times10^{'+str(int(m[2]))+'}$',s)
(r/'second_pass_results.tex').write_text(s)
files=[Path('/Users/haule/projects/eDMFT_new/src/dmft2')/f for f in ['stress_kernels.f90','stress_integrals.f90','stress_g2_matrices.f90','stress_interstitial.f90','alm_blm_clm.f90','lomain.f90','cmp_dmft_stress.f90','l2mainSubs.f90','cmp_dmft_weights.f90']]
files += [Path('/Users/haule/wien2k_24/SRC_lapw2')/f for f in ['val_kinstrmt.F','kin_str_integrals.F','stress_APW_surface.f','val_corr_stress.F']]
files += [Path('/Users/haule/dbin/dmft2'),Path('/Users/haule/wien2k_24/lapw2'),Path('/Users/haule/wien2k_24/lapw2c')]
(r/'SOURCE_AND_BINARY_SHA256.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2)+'\n')
print('Wrote report tables and provenance')
