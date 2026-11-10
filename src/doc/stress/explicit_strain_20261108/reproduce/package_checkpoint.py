"""Build reviewable source patch, guarded installer and compact evidence only."""
from pathlib import Path
import difflib,hashlib,json,shutil
import numpy as np

root=Path(__file__).resolve().parent
repo=Path('/Users/haule/projects/eDMFT_new');dbin=Path('/Users/haule/dbin')
old=Path('/private/tmp/edmft_projector_functional_20261108')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

road=(repo/'src/doc/stress/stress_response_roadmap_20261108.txt').read_text()
start=road.index('LATEST CHECKPOINT:');end=road.index('THE THREE LEVELS THAT CURRENTLY COEXIST')
road=road[:start]+'''LATEST CHECKPOINT: explicit_strain_20261108/REPORT.txt

We are in point 3. The complete native charge Jacobian J, numerical hydrostatic
charge source and q, and the projector-only functional readout are connected
and validated. The explicit hydrostatic POTENTIAL source is now analytic:
  n_eta = T_eta(n0 + eta*q),
  s_V = K*tau_n + Coulomb_geometry(n0,ions) + XC_geometry(n0),
  v_eta' = s_V + K*q.
Here tau_n is the derivative of the existing density transport; it must not
be added again to q. s_V is not the entire explicit CHARGE source s_charge.

The source includes fixed-RMT boundary/mask motion, periodic ionic Coulomb,
native pseudocharge/radial geometry, and LDA/PBE XC geometry. It uses mixed
native radial units (rV spherical, V_lm nonspherical) and the existing complex
VNS conjugation convention. A separate material-radial conversion is supplied
for diagnostics and is not silently inserted in the fixed-r response map.

201 tests pass with repository and dbin dependencies. Independent native
LAPW0 endpoint checks pass in SrVO3 and FeSe with step convergence; exact
errors and a small residual near the innermost FeSe spherical grid are
recorded in the checkpoint report. Source preparation costs about 0.2/0.5 s
after potential/core preparation, without native launches or SCF iteration.
An isolated native precision switch removed six-decimal lattice rounding
from the pure-hydrostatic validation. That switch is not installed in WIEN2k.

The projector readout uses
  g_V = B^T g_joint, c=K^T g_V, (I-J)^T lambda=c,
  dL_projector,potential = g_V^T s_V + lambda^T s_charge.
Explicit and induced pieces cancel substantially. This is still only the
projector potential contribution, not the full Klein derivative or total
stress. Diagnostic states were not newly reconverged and no new pressure
agreement or automatic stress publication is claimed.

NEXT WITHIN POINT 3: finish the explicit joint H/O/raw-projector geometry
source at fixed NATIVE potential coefficients. Reconcile its functional
contraction and direct radial terms with the existing analytic stress so that
the moving-mask term is counted exactly once. Validate the combined
hydrostatic derivative, then the fully analytic charge source, six strains,
and production publication. The wrapper still uses its native reference;
the prepared engine has not replaced that backend. MPI/streaming and the
deferred force audit remain as previously planned.

'''+road[end:]
road=road.replace('   the explicit geometric potential and remaining functional/baseline terms\n   must still be assembled once.',
    '   the hydrostatic explicit potential is now analytic and independently\n   checked. Remaining band geometry and functional/baseline terms must still\n   be assembled once.')
road=road.replace('   projector_functional_20261108/REPORT.txt for scope and numerical controls.',
    '   projector_functional_20261108/REPORT.txt and explicit_strain_20261108/REPORT.txt\n   for scope and numerical controls.')
road=road.replace('   Hydrostatic reference sources and a hydrostatic density-transport derivative\n   exist. General strain sources with fixed-RMT surface/transport terms, local',
    '   Numerical hydrostatic charge sources and an analytic hydrostatic potential\n   source exist. General strain sources with fixed-RMT surface/transport terms, local')
road=road.replace('   Complete total potential/readout validation is the immediate next step.',
    '   The analytic hydrostatic potential source and its sum with Kq now pass\n   independent native checks. Explicit band geometry and complete functional/\n   baseline reconciliation are the immediate next step.')
anchor='src/python/lapw_projector_functional.py\n'
assert anchor in road
road=road.replace(anchor,'''src/python/{lapw_strain_potential,lapw_hydrostatic_potential}.py
  Analytic hydrostatic potential source, native pseudocharge and ionic geometry,
  XC spatial scaling, moving interstitial mask, fixed-label potential files
  and explicit distinction between fixed-r and material-r potential derivatives.
'''+anchor,1)
road=road.replace('  projector_functional_20261108/REPORT.txt\n',
    '  projector_functional_20261108/REPORT.txt\n  explicit_strain_20261108/REPORT.txt\n',1)
(root/'stage/stress_response_roadmap_20261108.txt').write_text(road)

post=(old/'post_install_tests.py').read_text().replace('== 194','== 201')
post=post.replace("    'test_lapw_projector_functional',", "    'test_lapw_projector_functional', 'test_lapw_strain_potential',")
post=post.replace("helpers = ['lapw_projector_functional',", "helpers = ['lapw_strain_potential','lapw_hydrostatic_potential','lapw_projector_functional',")
post=post.replace("for name in ['test_lapw_projector_functional',", "for name in ['test_lapw_strain_potential','test_lapw_projector_functional',")
post=post.replace('all ten checked helpers and all five checked tests','all twelve checked helpers and all six checked tests')
(root/'post_install_tests.py').write_text(post)
install=(old/'install_checkpoint.py').read_text().replace('projector-functional readout patch','explicit hydrostatic potential-source patch').replace(
    'stress/projector_functional_20261108','stress/explicit_strain_20261108')
(root/'install_checkpoint.py').write_text(install)
record=(old/'record_installation.py').read_text().replace('stress/projector_functional_20261108','stress/explicit_strain_20261108').replace(
    '194','201').replace('all ten checked helpers and all five checked tests','all twelve checked helpers and all six checked tests').replace(
    'projector_only_readout_validated=True','analytic_hydrostatic_potential_validated=True').replace(
    'Point 3: explicit affine potential and complete functional/baseline contraction; then six analytic strain sources and production integration',
    'Point 3: explicit H/O/raw geometry and functional/baseline contraction; then fully analytic charge sources, six strains and production integration')
(root/'record_installation.py').write_text(record)

targets=[]
for name in ['lapw_strain_potential.py','lapw_hydrostatic_potential.py']:
    targets += [(root/'stage'/name,repo/'src/python'/name),(root/'stage'/name,dbin/name)]
targets += [(root/'stage/tests/test_lapw_strain_potential.py',repo/'src/python/tests/test_lapw_strain_potential.py'),
            (root/'stage/stress_response_roadmap_20261108.txt',repo/'src/doc/stress/stress_response_roadmap_20261108.txt')]
manifest=[];patch=''
for source,dest in targets:
    manifest.append(dict(source=str(source),destination=str(dest),source_sha256=sha(source),before_sha256=sha(dest) if dest.exists() else None))
    if dest.is_relative_to(repo):
        name=str(dest.relative_to(repo));before=dest.read_text().splitlines(True) if dest.exists() else []
        patch+=''.join(difflib.unified_diff(before,source.read_text().splitlines(True),fromfile='a/'+name if dest.exists() else '/dev/null',tofile='b/'+name))
(root/'explicit_strain.patch').write_text(patch)
(root/'installation_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

evidence=root/'evidence';evidence.mkdir(exist_ok=True)
def copy(relative):
    dest=evidence/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(root/relative,dest)
for name in ['REPORT.txt','explicit_strain.patch','installation_manifest.json','tests_repository.log','tests_dbin.log',
             'native_build.log','native_readstruct_precision.patch']:
    copy(name)
for case in ['SrVO3','FeSe']:
    for name in [case+'.log',case+'_analytic.log',case+'_total_analytic.log',case+'/REPORT.json']:
        copy(name)
    for kind in ['pure_hydro','total_hydro']:
        copy(case+'/'+kind+'/REPORT.json')
        for step in ['0.00025','0.000125']:
            relative=case+'/'+kind+'/comparison_'+step+'.npz'
            data=np.load(root/relative)
            dest=evidence/relative;dest.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(dest,**{k:data[k] for k in ['analytic','native','source','gradient','induced_from_previous_q']})
        for path in sorted((root/case/kind).glob('*/lapw0.log')):copy(str(path.relative_to(root)))
    data=np.load(root/case/'pure_hydro/comparison_0.000125.npz')
    np.savez_compressed(evidence/case/'hydrostatic_parts.npz',**{k:data[k] for k in
        ['transport','coulomb_geometry','xc_geometry','central_coulomb','central_xc','central_potential','mask_derivative']})
    for path in sorted((root/case).glob('potential_split_*.npz')):
        np.savez_compressed(evidence/case/path.name,**dict(np.load(path)))
for name in ['check_explicit_potential.py','check_analytic_source.py','run_tests.py','post_install_tests.py',
             'install_checkpoint.py','record_installation.py','package_checkpoint.py',
             'prepare_native_reference.py','verify_checkpoint.py']:
    dest=evidence/'reproduce'/name;dest.parent.mkdir(exist_ok=True);shutil.copy2(root/name,dest)
for name in ['native/SRC_lapw0/lapw0','native/SRC_lapw0/readstruct.f']:
    copy(name)
for path in (root/'stage').rglob('*'):
    if path.is_file():copy(str(path.relative_to(root)))
provenance={}
for case in ['SrVO3','FeSe']:
    previous=Path('/private/tmp/edmft_screened_native_20261108/full_bz')/case
    paths=[previous/'central'/(case+suffix) for suffix in ['.struct','.clmsum','.in0','.inm','.vsp','.vns']]
    paths += [previous/'analytic_density.npy',Path('/private/tmp/edmft_projector_functional_20261108/full_bz')/case/'covectors.npz',
              repo/'src/doc/stress/charge_closure_20261108'/(case+'_potential_core.npz')]
    provenance[case]={str(p):sha(p) for p in paths}
provenance['native']={str(p):sha(p) for p in [root/'native/SRC_lapw0/lapw0',root/'native/SRC_lapw0/readstruct.f']}
(evidence/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n')
hashes={str(p.relative_to(evidence)):sha(p) for p in sorted(evidence.rglob('*')) if p.is_file() and p.name!='SHA256.json'}
(evidence/'SHA256.json').write_text(json.dumps(hashes,indent=2)+'\n')
print('Prepared',len(manifest),'destinations and',len(hashes),'evidence files; nothing installed.')
