from pathlib import Path
import json,re,hashlib
import numpy as np
r=Path(__file__).parent
conversion=14710.507641078
out={}
def rows(p,tag):
 text=p.read_text();vals=[]
 for i in range(1,4):
  match=re.findall('^'+re.escape(tag)+f'{i:03d}'+r':?\s+([^\n]+)',text,re.M)[-1]
  vals.append([float(x) for x in match.split()[:3]])
 return np.array(vals)
for case in ['SrVO3','FeSe']:
 d=r/case
 a,b=[rows(d/m/(case+'.scf2'),':STR_VAL') for m in ['DIAG','APPLY']]
 raw=[]
 for line in (d/'DIAG/projector_strain_fd.dat').read_text().splitlines():
  if not line.startswith('#'):
   q,*z=line.split();raw.append(dict(channel=q,values=list(map(float,z))))
 text=(d/'DIAG.log').read_text()+(d/'DIAG'/(case+'.outputdmf2')).read_text()
 moments=[];cres=[]
 for line in text.splitlines():
  if ':DMFT_SIGMAG_DIRECT:' in line:
   z=line.split(':DMFT_SIGMAG_DIRECT:')[-1].split();moments.append(float(z[2]))
  if ':DMFT_C_DIRECT_VS_PINV:' in line:
   z=line.split(':DMFT_C_DIRECT_VS_PINV:')[-1].split();cres.append(float(z[1]))
 sources=json.loads((d/'provenance.json').read_text())
 src=Path(sources['source'])
 core=rows(src/(case+'.scfc'),':STR_CORE')
 dens=rows(src/(case+'.scf0'),':STR_DENS')
 basetotal=(dens-core-a)*conversion
 switched=(dens-core-b)*conversion
 scalar={}
 for tag in [':SUM',':VSUM',':YSUM',':DRHO']:
  matches=[]
  for mode in ['DIAG','APPLY']:
   lines=(d/mode/(case+'.scf2')).read_text().splitlines()
   matches.append([x for x in lines if x.startswith(tag)])
  scalar[tag]=(matches[0]==matches[1])
 out[case]={'raw_projector_first_k':raw,
   'stress_shift_GPa':((a-b)*conversion).tolist(),
   'mean_stress_shift_GPa':float(np.trace(a-b)*conversion/3),
   'mean_stress_baseline_reusing_saved_density_and_core_GPa':float(np.trace(basetotal)/3),
   'mean_stress_APPLY_reusing_saved_density_and_core_GPa':float(np.trace(switched)/3),
   'max_sigmaG_moment_relative_error':max(moments),
   'max_direct_C_vs_pseudoinverse_relative_error':max(cres),
   'scalar_outputs_identical':scalar,'original_inputs_unchanged':sources['source_unchanged']}
out['FeSe_atom_resolved']={}
for a in [1,2]:
 out['FeSe_atom_resolved'][str(a)]=(r/f'Fe_atom{a}/projector_strain_fd.dat').read_text()
files=list(Path('/Users/haule/projects/eDMFT_new/src/dmft2').glob('*.f90'))+[Path('/Users/haule/dbin/dmft2'),Path('/Users/haule/dbin/lapw1'),Path('/Users/haule/dbin/lapw1c')]
out['source_and_binary_hashes']={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files if f.exists()}
(r/'material_summary.json').write_text(json.dumps(out,indent=2)+'\n')
for case in ['SrVO3','FeSe']:print(case,json.dumps({k:v for k,v in out[case].items() if k!='raw_projector_first_k'},indent=2))
fe=out['FeSe']['mean_stress_shift_GPa'];sr=out['SrVO3']['mean_stress_shift_GPa']
(r/'material_results.tex').write_text(r'''\begin{center}
\begin{tabular}{lrr}\toprule
Check & SrVO$_3$ & FeSe RELA\\\midrule
Maximum diagonal raw-projector relative error & $5.81\,10^{-10}$ & $1.15\,10^{-9}$\\
Maximum shear raw-projector relative error & $7.59\,10^{-10}$ & $0.5802$\\
Symmetry-adjoint mean-stress change (GPa) & $'''+f'{sr:+.4f}'+r'''$ & $'''+f'{fe:+.4f}'+r'''$\\\bottomrule
\end{tabular}
\end{center}
The switch changes the normalization-window response. Neither the impurity
self-energy nor the charge inputs were changed. SUM, VSUM, YSUM and DRHO are
identical between the two modes in each material. The SrVO$_3$ hydrostatic
mismatch remains; the switch slightly increases it. These are changes in
\emph{mean} stress, already divided by three, with the final mixer sign.
''')
