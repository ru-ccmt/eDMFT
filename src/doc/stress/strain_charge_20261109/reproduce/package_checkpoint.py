"""Package the reviewed-source candidate; no writes to repository or dbin."""
from pathlib import Path
import difflib,hashlib,json,shutil
import numpy as np
root=Path(__file__).resolve().parent
repo=Path('/Users/haule/projects/eDMFT_new');dbin=Path('/Users/haule/dbin')
old=Path('/private/tmp/edmft_band_geometry_20261108')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

report=(root/'REPORT_body.txt').read_text()
report+='Relative errors of the normalized material CHARGE source, not GPa.\n'
report+='Case/mesh      small-step source   Richardson source   Richardson PW    q difference\n'
timing=[]
for mesh in ['three_k','full_bz']:
 for case in ['SrVO3','FeSe']:
  r=json.loads((root/mesh/case/'REPORT.json').read_text())
  report+=f"{case:6s}/{mesh:7s} {r['steps']['0.000125']['comparisons']['source']['all']['relative']:14.6e} {r['richardson_source_relative']:19.6e} {r['richardson_PW_relative']:16.6e} {r['screened']['relative_to_numerical_source_solution']:14.6e}\n"
  timing.append(f"{case}/{mesh}: valence preparation {r['valence_preparation_seconds']:.2f} s; source {r['seconds']:.2f} s; warm-start solve {r['screened']['seconds']:.2f} s, {r['screened']['actions']} actions; source charge {r['source_charge']:.3e}; true relative residual {r['screened']['residual_norm']/r['source_norm']:.3e}; central fixed-point residual norm {r['metadata']['central_fixed_point_residual_norm']:.6g}.")
report+='\nAdditional three-k step 0.0000625:\n'
for case in ['SrVO3','FeSe']:
 r=json.loads((root/'three_k'/case/'REPORT.json').read_text())['third_step']
 report+=f"{case}: source relative {r['source_relative']:.6e}; Richardson {r['richardson_relative']:.6e}; Richardson PW {r['richardson_PW_relative']:.6e}.\n"
report+='\n'+ '\n'.join(timing)+'\n'
report+='\nThe source includes large radial-coordinate components; the separately listed\nPW errors are the relevant limits for that smaller sector. Their residual\norigin remains open; the test does not establish roundoff-level agreement.\n'
(root/'REPORT.txt').write_text(report)

road=(repo/'src/doc/stress/stress_response_roadmap_20261108.txt').read_text()
start=road.index('LATEST CHECKPOINT:');end=road.index('THE THREE LEVELS THAT CURRENTLY COEXIST')
road=road[:start]+'''LATEST CHECKPOINT: strain_charge_20261109/REPORT.txt

The HYDROSTATIC portion of point 3 is now implemented and checked: full native
J and J^T, explicit analytic potential and band-geometry sources, complete
analytic charge source s_nu, and screened q_nu from (I-J)q_nu=s_nu.
The full source combines electronic/normalized-projector response, explicit
spatial MT/PW valence geometry, radial/core response, normalization geometry,
inverse density transport and the reference's second central normalization.
Inverse transport acts on the central normalized OUTPUT, not input n0.

215 tests pass with repository and dbin dependencies. Independent native
source comparisons cover SrVO3 NREL (3 and 75 k points) and FeSe RELA (3 and
140 k points), with two strain sizes plus a third size on the small meshes.
The report records separate MT/PW errors, residuals and timings. Small
interstitial residuals remain: do not infer roundoff accuracy or pressure
agreement from the much smaller overall radial-dominated relative error.
These are local tangent tests at archived states with recorded nonzero
fixed-point residuals, not newly converged physical stress calculations.

NEXT: complete the functional contraction and baseline reconciliation (point
4), including forward/adjoint material readout and zero-Sigma controls.
The projector-only gradient exists; the full Klein/direct/radial gradient and
subtraction of terms already present in analytic stress remain unfinished.
Keep the interstitial potential mask in the explicit potential source and
the charge integration mask in normalization, each exactly once.

Six general physical strain sources still require symmetry/harmonic/site
layout expansion and a consistent anisotropic transport convention. The new
diagonal field adapter matches actual rounded hydrostatic native structures;
it is not general shear support. Production backend replacement, final tensor
publication, MPI/streaming and the deferred force audit remain later tasks.
No new total Klein derivative or pressure correction is claimed here.

'''+road[end:]
road=road.replace('   Spatial charge geometry and complete functional/baseline terms remain.',
 '   The full hydrostatic charge source is now analytic and tested. Complete\n   functional/baseline terms and six-source layout expansion remain.')
road=road.replace('   available. Complete analytic spatial-charge source assembly comes next.',
 '   available, as is complete analytic hydrostatic source assembly. Complete\n   functional/baseline reconciliation comes next.')
road=road.replace('   Complete analytic spatial-charge source assembly comes next, followed by\n   full functional/baseline reconciliation.',
 '   The complete hydrostatic charge source and its screened solution now pass\n   native checks too. Full functional/baseline reconciliation comes next.')
road=road.replace('   strain sources now include tested explicit H/O/raw geometry; spatial\n   charge geometry and complete functional reconciliation remain.',
 '   strain sources now include tested explicit H/O/raw geometry and full\n   hydrostatic spatial-charge assembly. General tensor layouts and complete\n   functional reconciliation remain.')
road=road.replace('   Numerical hydrostatic charge sources and an analytic hydrostatic potential\n   source exist. General strain sources with fixed-RMT surface/transport terms, local',
 '   The complete hydrostatic charge source is now analytic and validated.\n   General strain sources with fixed-RMT surface/transport terms, local')
start=road.index('THE IMMEDIATE DELIVERABLE');end=road.index('SOURCE AND EVIDENCE MAP',start)
road=road[:start]+'''THE IMMEDIATE DELIVERABLE

The analytic hydrostatic source and q are now available. Next assemble and
validate the COMPLETE functional contraction against the existing analytic
stress baseline: projector, radial and direct terms in shared conventions.
First verify the forward/adjoint observable contractions using the new source;
then zero-Sigma DFT and consistent frozen-Sigma functional finite differences.
Do not convert isolated projector or charge-source norms into pressure claims.
Six general strain layouts and production publication still follow.

'''+road[end:]
anchor='src/python/lapw_projector_functional.py\n'
road=road.replace(anchor,'''src/python/{lapw_strain_charge,lapw_strain_field}.py
  Complete hydrostatic unmixed charge source: explicit spatial valence geometry,
  fixed-RMT integration/normalization geometry, inverse output transport and
  post-transport normalization. Diagonal Coulomb/XC geometry matches rounded
  native probes with explicit scalar density transport. No stress publication.
'''+anchor,1)
road=road.replace('  band_geometry_20261108/REPORT.txt\n','  band_geometry_20261108/REPORT.txt\n  strain_charge_20261109/REPORT.txt\n',1)
(root/'stage/stress_response_roadmap_20261108.txt').write_text(road)

install=(old/'install_checkpoint.py').read_text().replace('explicit LAPW band-geometry patch','analytic hydrostatic charge-source patch').replace('stress/band_geometry_20261108','stress/strain_charge_20261109')
(root/'install_checkpoint.py').write_text(install)
record=(old/'record_installation.py').read_text().replace('stress/band_geometry_20261108','stress/strain_charge_20261109').replace('209','215').replace('all thirteen checked helpers and all seven checked tests','all fifteen checked helpers and all nine checked tests').replace('analytic_band_geometry_validated=True','analytic_hydrostatic_charge_source_validated=True').replace('Point 3: spatial charge geometry, normalization and inverse transport; compare complete analytic source and screened q, then functional/baseline contraction and six strains','Hydrostatic point 3 validated; next complete functional/baseline contraction, six strain layouts and production integration')
(root/'record_installation.py').write_text(record)
targets=[]
for name in ['lapw_strain_charge','lapw_strain_field']:
 targets.extend([(root/'stage'/(name+'.py'),repo/'src/python'/(name+'.py')),(root/'stage'/(name+'.py'),dbin/(name+'.py')),
                 (root/'stage/tests'/('test_'+name+'.py'),repo/'src/python/tests'/('test_'+name+'.py'))])
targets.append((root/'stage/stress_response_roadmap_20261108.txt',repo/'src/doc/stress/stress_response_roadmap_20261108.txt'))
manifest=[];patch=''
for source,dest in targets:
 manifest.append(dict(source=str(source),destination=str(dest),source_sha256=sha(source),before_sha256=sha(dest) if dest.exists() else None))
 if dest.is_relative_to(repo):
  relative=str(dest.relative_to(repo));before=dest.read_text().splitlines(True) if dest.exists() else []
  patch+=''.join(difflib.unified_diff(before,source.read_text().splitlines(True),fromfile='a/'+relative if dest.exists() else '/dev/null',tofile='b/'+relative))
(root/'strain_charge.patch').write_text(patch)
(root/'installation_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

evidence=root/'evidence';evidence.mkdir(exist_ok=True)
def copy(relative):
 dest=evidence/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(root/relative,dest)
for name in ['REPORT.txt','REPORT_body.txt','strain_charge.patch','installation_manifest.json','tests_repository.log','tests_dbin.log','refine_native.log','refinement.log']:
 copy(name)
for mesh in ['three_k','full_bz']:
 for case in ['SrVO3','FeSe']:
  prefix=Path(mesh)/case
  for name in ['REPORT.json','analytic_source.npy','analytic_density.npy','analytic_induced_potential.npy']:
   copy(str(prefix/name))
  copy(f'{case}_{mesh}.log')
  arrays=np.load(root/prefix/'analytic.npz')
  np.savez_compressed(evidence/prefix/'analytic.npz',**{k:arrays[k] for k in arrays.files if k not in ['band_geometry','joint','electronic_density']})
  for p in (root/prefix).glob('native_*.npz'):copy(str(p.relative_to(root)))
for case in ['SrVO3','FeSe']:
 for p in (root/'native_extra'/case).rglob('*'):
  if p.is_file() and p.name in ['done.json','timings.json','lapw0.log','lapw1.log','dmft2.log','dmft2_fixed_n.log','lcore.log']:
   copy(str(p.relative_to(root)))
for name in ['check_charge.py','refine_native.py','check_refinement.py','finalize_material_reports.py','run_tests.py','post_install_tests.py',
             'package_checkpoint.py','verify_checkpoint.py','install_checkpoint.py','record_installation.py']:
 dest=evidence/'reproduce'/name;dest.parent.mkdir(exist_ok=True);shutil.copy2(root/name,dest)
for path in (root/'stage').rglob('*'):
 if path.is_file() and path.suffix in ['.py','.txt']:copy(str(path.relative_to(root)))
provenance={}
for mesh in ['three_k','full_bz']:
 for case in ['SrVO3','FeSe']:
  ref=Path('/private/tmp/edmft_screened_native_20261108')/mesh/case
  source=Path('/private/tmp/edmft_full_k_response_20261108')
  if mesh=='full_bz':source=source/mesh
  source=source/case/'charge/0'
  paths=[source/(case+x) for x in ['.struct','.clmsum','.clmval','.vsp','.vns','.indmfl','.in1c' if case=='SrVO3' else '.in1']]
  paths += [source/'projectorw.dat',ref/'analytic_density.npy']
  paths += [p for p in ref.glob('source_*/*') if p.name in ['done.json','normalized.npy']]
  provenance[case+'/'+mesh]={str(p):sha(p) for p in paths}
(evidence/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n')
hashes={str(p.relative_to(evidence)):sha(p) for p in sorted(evidence.rglob('*')) if p.is_file() and p.name!='SHA256.json'}
(evidence/'SHA256.json').write_text(json.dumps(hashes,indent=2)+'\n')
print(f'Prepared {len(manifest)} targets and {len(hashes)} evidence files; nothing installed.')
