"""Package the staged band-geometry patch; never writes repository or dbin."""
from pathlib import Path
import difflib
import hashlib
import json
import shutil
import numpy as np

root = Path(__file__).resolve().parent
repo = Path('/Users/haule/projects/eDMFT_new')
dbin = Path('/Users/haule/dbin')
old = Path('/private/tmp/edmft_explicit_strain_20261108')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

road = (repo/'src/doc/stress/stress_response_roadmap_20261108.txt').read_text()
start = road.index('LATEST CHECKPOINT:')
end = road.index('THE THREE LEVELS THAT CURRENTLY COEXIST')
road = road[:start]+'''LATEST CHECKPOINT: band_geometry_20261108/REPORT.txt

We remain in point 3. The full native charge Jacobian J, numerical hydrostatic
charge source and screened q, projector-only functional readout, and analytic
hydrostatic potential source s_V have been implemented and checked.

The explicit BAND geometry source is now analytic and independently checked:
  b_total = b_geometry + B (s_V + K q),
  b = (delta H_RR, delta O_RR, delta Phi_R, delta H_QR, delta O_QR).
R means retained/exported bands, Q the remaining finite LAPW complement;
Phi is RAW projector overlap. The electronic map applies metric/window and
normalized-projector response afterwards, exactly once. b_geometry holds
native potential FILE coefficients fixed. The moving potential mask belongs
in s_V; this convention must agree with the later stress-baseline subtraction.

Native three-k SrVO3 NREL and FeSe RELA tests validate H, O and raw-projector
geometry. Errors decrease by four when the strain step is halved; joint
Richardson errors are 6.95e-11 and 3.88e-10. 209 tests pass with repository
and dbin dependencies. The material adapter is limited to symmetry-preserving
strains of the supported primitive orthogonal snapshots, fixed RMT, radial
energies and local axes. General symmetry-breaking/shear layouts are pending.

NEXT WITHIN POINT 3: finish the explicit spatial VALENCE geometry, native
output-normalization geometry and inverse density transport. Combine them
with b_geometry and s_V to form the full analytic s_charge. Compare this
source and its screened q against the existing numerical native references.

After that, reconcile the complete Klein/direct/radial readout with the
current analytic stress baseline, then six physical strain sources and
production publication. The isolated projector contractions are not total
Klein derivatives or pressure corrections. No new total-stress agreement is
claimed. The production wrapper still uses its earlier native reference.
MPI/streaming and the deferred force audit remain on the roadmap.

'''+road[end:]
road=road.replace('Current implementation roadmap, 8 November 2026','Current implementation roadmap, 9 November 2026')
road=road.replace('   strain sources remain to be connected.', '   strain sources now include tested explicit H/O/raw geometry; spatial\n   charge geometry and complete functional reconciliation remain.')
road=road.replace('   baseline reconciliation are the immediate next step.', '   baseline reconciliation remain; the tested explicit band geometry is now\n   available. Complete analytic spatial-charge source assembly comes next.')
road=road.replace('   checked. Remaining band geometry and functional/baseline terms must still\n   be assembled once.',
 '   checked. Explicit band geometry is now independently checked as well.\n   Spatial charge geometry and complete functional/baseline terms remain.')
road=road.replace('   Independent charge integration/gauge tests, reusable loading and production\n   MPI wiring remain pending.',
 '   Charge/gauge and broader-k native checks plus reusable loading are now\n   available; production MPI wiring remains pending.')
start=road.index('1. Finish the valence-potential adapter;')
end=road.index('3. Solve and validate screened response',start)
road=road[:start]+'''1. Valence-potential adapter: COMPLETED FOR THE TESTED NATIVE SCOPE.

   Spherical, nonspherical and interstitial input sectors, their radial and
   explicit basis-density contributions, the complete finite LAPW complement,
   and their exact transposes are connected. Reusable full-k native loading,
   RR/QR memory layout, number/gauge checks and independent native endpoints
   are documented in shared_gamma_response and screened_native checkpoints.
   Broader geometry/energy policies and production MPI are later extensions.

2. Potential/core response and full charge Jacobian: COMPLETED FOR TESTED SCOPE.

   Periodic Hartree, native nonmagnetic LDA/PBE XC, NREL/RELA core, native
   output-charge normalization, and reusable central archives are connected
   to the valence map. J and J^T pass native tests in SrVO3 and FeSe, including
   broader-k forward and screened-response checks. The native reference stays
   available for independent comparisons. Spin/SOC and general-cell coverage,
   automatic energy searches and production MPI remain outside this scope.

'''+road[end:]
road=road.replace('   independent native checks. Explicit band geometry and complete functional/\n   baseline reconciliation remain; the tested explicit band geometry is now\n   available. Complete analytic spatial-charge source assembly comes next.',
 '   independent native checks. Explicit band geometry is now tested too.\n   Complete analytic spatial-charge source assembly comes next, followed by\n   full functional/baseline reconciliation.')
start=road.index('THE IMMEDIATE DELIVERABLE')
end=road.index('SOURCE AND EVIDENCE MAP',start)
road=road[:start]+'''THE IMMEDIATE DELIVERABLE

An entirely analytic hydrostatic charge source s_charge in the SAME central
material coordinates used by the checked numerical native source. The full
native J and potential/radial response maps are already connected. Band
geometry and the hydrostatic potential source are now separately validated.
Explicit spatial valence geometry, native normalization geometry and inverse
transport must still be combined. At a nonzero fixed-point residual, inverse
transport acts on the central OUTPUT charge, not an assumed input equal to it.
Then solve (I-J)q=s_charge and compare with the existing native source/solution.
This is still point 3; total functional/baseline assembly and production tensor
publication follow. See dated checkpoints for evidence and precise scopes.

'''+road[end:]
anchor='src/python/lapw_projector_functional.py\n'
road=road.replace(anchor,'''src/python/lapw_geometry_response.py
  Fixed-native-potential H/O/raw geometry source: reciprocal/matching and LO
  derivatives, spherical surface and ordered Wronskian, nonspherical legs,
  interstitial overlap/kinetic and fixed local-frame raw projector overlaps.
  RR/QR columns are preserved. This is not yet a full spatial-charge source.
'''+anchor,1)
road=road.replace('  explicit_strain_20261108/REPORT.txt\n','  explicit_strain_20261108/REPORT.txt\n  band_geometry_20261108/REPORT.txt\n',1)
(root/'stage/stress_response_roadmap_20261108.txt').write_text(road)

install=(old/'install_checkpoint.py').read_text().replace('explicit hydrostatic potential-source patch','explicit LAPW band-geometry patch').replace('stress/explicit_strain_20261108','stress/band_geometry_20261108')
(root/'install_checkpoint.py').write_text(install)
record=(old/'record_installation.py').read_text().replace('stress/explicit_strain_20261108','stress/band_geometry_20261108').replace('201','209').replace('all twelve checked helpers and all six checked tests','all thirteen checked helpers and all seven checked tests').replace('analytic_hydrostatic_potential_validated=True','analytic_band_geometry_validated=True').replace('Point 3: explicit H/O/raw geometry and functional/baseline contraction; then fully analytic charge sources, six strains and production integration','Point 3: spatial charge geometry, normalization and inverse transport; compare complete analytic source and screened q, then functional/baseline contraction and six strains')
(root/'record_installation.py').write_text(record)

targets=[
 (root/'stage/lapw_geometry_response.py',repo/'src/python/lapw_geometry_response.py'),
 (root/'stage/lapw_geometry_response.py',dbin/'lapw_geometry_response.py'),
 (root/'stage/tests/test_lapw_geometry_response.py',repo/'src/python/tests/test_lapw_geometry_response.py'),
 (root/'stage/stress_response_roadmap_20261108.txt',repo/'src/doc/stress/stress_response_roadmap_20261108.txt'),
]
manifest=[];patch=''
for source,dest in targets:
 manifest.append(dict(source=str(source),destination=str(dest),source_sha256=sha(source),before_sha256=sha(dest) if dest.exists() else None))
 if dest.is_relative_to(repo):
  rel=str(dest.relative_to(repo));before=dest.read_text().splitlines(True) if dest.exists() else []
  patch+=''.join(difflib.unified_diff(before,source.read_text().splitlines(True),fromfile='a/'+rel if dest.exists() else '/dev/null',tofile='b/'+rel))
(root/'band_geometry.patch').write_text(patch)
(root/'installation_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

evidence=root/'evidence';evidence.mkdir(exist_ok=True)
def copy(relative):
 dest=evidence/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(root/relative,dest)
for name in ['REPORT.txt','band_geometry.patch','installation_manifest.json','tests_repository.log','tests_dbin.log']:
 copy(name)
provenance={}
for case in ['SrVO3','FeSe']:
 for name in [case+'/REPORT.json',case+'/timings.json',case+'_native.log',case+'_analytic.log']:
  copy(name)
 for step in ['0.00025','0.000125']:
  data=dict(np.load(root/case/f'comparison_{step}.npz'))
  np.savez_compressed(evidence/case/f'comparison_{step}.npz',**data)
  for label in ['m','p']:
   folder=root/case/f'{label}_{step}'
   for name in ['endpoint.json','lapw1.log','dmft2.log']:
    copy(str((folder/name).relative_to(root)))
   paths=list(folder.glob('native*.bin'))+list((folder/'response').glob('*'))
   provenance[str(folder)]={str(p):sha(p) for p in paths if p.is_file()}
for name in ['check_material.py','prepare_native.py','run_tests.py','post_install_tests.py',
             'package_checkpoint.py','verify_checkpoint.py','install_checkpoint.py','record_installation.py']:
 dest=evidence/'reproduce'/name;dest.parent.mkdir(exist_ok=True);shutil.copy2(root/name,dest)
for path in (root/'stage').rglob('*'):
 if path.is_file() and path.suffix in ['.py','.txt']:copy(str(path.relative_to(root)))
native=Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')
paths=[native/'lapw1',native/'lapw1c',dbin/'dmft2_response',dbin/'liblapw_response.so']
paths += [Path(p) for p in [
 '/private/tmp/edmft_charge_closure_20261108/native/libmixer_normalization.dylib',
 '/private/tmp/edmft_hartree_radial_20261108/native/libhartree_reference.dylib',
 '/private/tmp/edmft_hartree_periodic_20261108/native/libcharge_reference.dylib']]
provenance['native_binaries']={str(p):sha(p) for p in paths}
(evidence/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n')
hashes={str(p.relative_to(evidence)):sha(p) for p in sorted(evidence.rglob('*')) if p.is_file() and p.name!='SHA256.json'}
(evidence/'SHA256.json').write_text(json.dumps(hashes,indent=2)+'\n')
print(f'Prepared {len(manifest)} targets and {len(hashes)} evidence files; nothing installed.')
