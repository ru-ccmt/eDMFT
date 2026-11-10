"""Prepare a reviewable checkpoint in scratch; do not install anything."""
from pathlib import Path
import difflib
import hashlib
import json
import shutil

root = Path(__file__).resolve().parent
repo = Path('/Users/haule/projects/eDMFT_new')
dbin = Path('/Users/haule/dbin')
old = Path('/private/tmp/edmft_screened_native_20261108')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

roadmap = (repo / 'src/doc/stress/stress_response_roadmap_20261108.txt').read_text()
start = roadmap.index('LATEST CHECKPOINT:')
end = roadmap.index('THE THREE LEVELS THAT CURRENTLY COEXIST')
roadmap = roadmap[:start] + '''LATEST CHECKPOINT: projector_functional_20261108/REPORT.txt

We are in point 3. NativeScreenedState connects the complete charge Jacobian
  J = N [R E B + D + C] K,   (I-J) q_nu = s_nu,
with native hydrostatic source and full-mesh screening validation recorded in
screened_native_20261108/. s_nu is the source called b_nu in earlier reports.

PreparedProjectorFunctional now adds the potential/charge readout of the extra
NORMALIZED-PROJECTOR dependence of the regularized lattice scalar:
  g_v = B^T g_joint,   c = K^T g_v,
  (I-J)^T lambda = c,
  dL_projector,potential,nu = g_v^T v_explicit,nu + lambda^T s_nu.
It includes the moving overlap metric, global X and actual orbital symmetry,
and the empty complement of the stored KS space using prepared shifted solves.
One adjoint serves multiple sources. Explicit potential inputs are required.

194 tests pass with repository and dbin dependencies. Independent three-k
native scalar endpoints validate the new gradient in FeSe and SrVO3; actual
projector-readout forward/adjoint contractions agree within true residual bounds.
Full configured meshes (140/75 k points) prepare the common gradient and pass
the potential-to-charge pullback check. Exact data and scope are in the report.

This is the PROJECTOR-ONLY scalar readout, not the full Klein-functional
gradient. It excludes direct density/Hamiltonian and baseline terms. The
induced potential Kq is only part of the total potential. Partial contractions,
especially FeSe's, can be large and must not be interpreted as total pressure
corrections. Diagnostic central states were not newly reconverged. No total
stress/GPa comparison or automatic stress publication is claimed.

NEXT WITHIN POINT 3: define the explicit affine/ionic potential in the same
native coordinates, assemble its projector contraction and complementary
direct radial H/O/geometry terms against the existing analytic stress, and
validate the COMBINED hydrostatic derivative without double counting.
Then complete six analytic geometry sources and production publication.
The production wrapper still uses its native reference; the prepared engine
has not replaced that backend. MPI/streaming remains necessary for large
meshes. The deferred force audit stays after strain validation.

''' + roadmap[end:]
old_h = '''   test the connection. A physical muffin-tin charge observable tests the
   adjoint and reuse for two sources on small meshes. This does not yet
   provide the functional/stress gradient. Kq is explicitly identified as
   the induced potential; an explicit geometric potential must be added once.
   See screened_native_20261108/REPORT.txt for scope and numerical controls.
'''
new_h = '''   test the connection. A muffin-tin charge observable first tested adjoint
   reuse. The projector-only functional readout now also passes native scalar
   endpoints and forward/adjoint checks; its metric and empty-complement
   pullbacks reuse the same native operators. Kq is the induced potential;
   the explicit geometric potential and remaining functional/baseline terms
   must still be assembled once. See screened_native_20261108/REPORT.txt and
   projector_functional_20261108/REPORT.txt for scope and numerical controls.
'''
assert old_h in roadmap
roadmap = roadmap.replace(old_h, new_h, 1)
roadmap = roadmap.replace('   Native export and the full-k valence loader are still separate work. Broader\n   number/gauge and k-point controls remain before physical screening.',
    '   Reusable full-k valence loading and physical screening are now connected\n   in point 3. Broader material/spin/geometry and number-policy controls remain.')
roadmap = roadmap.replace('WHAT THE MOST RECENT TESTS ESTABLISH\n',
    'CHECKPOINT HISTORY (OLDER LIMITATIONS MAY BE SUPERSEDED ABOVE)\n', 1)
old_step = '''   The software for this linear solve exists; the full physical J is what
   steps 1--2 provide. Start with hydrostatic FeSe and SrVO3 against the current
   reference delta n and delta V_KS. Record actual residuals and timings.
'''
new_step = '''   The complete native J, numerical hydrostatic sources and screened q now
   pass independent native checks in FeSe and SrVO3. True residuals, central
   fixed-point residuals and timings are recorded in screened_native_20261108.
   The projector-only readout passes its own native scalar/adjoint checks.
   Complete total potential/readout validation is the immediate next step.
'''
assert old_step in roadmap
roadmap = roadmap.replace(old_step,new_step,1)
roadmap = roadmap.replace('   part has already been tested; the full material gradient has not.',
    '   potential gradient now includes complete finite-basis metric and radial\n   projector dependence and passes native material checks. The complete\n   functional/baseline gradient has not yet been assembled.')
anchor = 'src/python/lapw_complement_response.py\n'
assert anchor in roadmap
roadmap = roadmap.replace(anchor, '''src/python/lapw_projector_functional.py
  Extra normalized-projector scalar readout, full finite-basis metric/QR
  pullbacks, one common native potential/charge gradient and adjoint reuse.
  Does not supply the complete stress or choose the baseline subtraction.
''' + anchor,1)
roadmap = roadmap.replace('  shared_gamma_response_20261108/REPORT.txt\n',
    '  shared_gamma_response_20261108/REPORT.txt\n'
    '  screened_native_20261108/REPORT.txt\n'
    '  projector_functional_20261108/REPORT.txt\n',1)
(root / 'stage/stress_response_roadmap_20261108.txt').write_text(roadmap)

# Prepare installed-only verification now, to avoid ambiguous stage imports.
post = (old / 'post_install_tests.py').read_text().replace('== 188','== 194')
post = post.replace("'test_lapw_core_response', 'test_lapw_xc_response', 'test_lapw_screened_state',",
    "'test_lapw_core_response', 'test_lapw_xc_response', 'test_lapw_screened_state',\n"
    "    'test_lapw_projector_functional',")
post = post.replace("helpers = ['lapw_screened_state',", "helpers = ['lapw_projector_functional', 'lapw_screened_state',")
post = post.replace("for name in ['test_lapw_screened_state',", "for name in ['test_lapw_projector_functional', 'test_lapw_screened_state',")
post = post.replace('all nine checked helpers and all four checked tests', 'all ten checked helpers and all five checked tests')
(root / 'post_install_tests.py').write_text(post)

installer = (old / 'install_checkpoint.py').read_text().replace(
    'native screened-charge connection patch', 'projector-functional readout patch').replace(
    'stress/screened_native_20261108', 'stress/projector_functional_20261108')
(root / 'install_checkpoint.py').write_text(installer)
record = (old / 'record_installation.py').read_text().replace(
    'stress/screened_native_20261108', 'stress/projector_functional_20261108').replace('188','194').replace(
    'all nine checked helpers and all four checked tests', 'all ten checked helpers and all five checked tests')
record = record.replace("physical_screened_solve_started=True,", "projector_only_readout_validated=True,")
record = record.replace('Point 3: physical functional covector and total potential/stress contraction; then general strain sources and production integration',
    'Point 3: explicit affine potential and complete functional/baseline contraction; then six analytic strain sources and production integration')
(root / 'record_installation.py').write_text(record)

targets = [
    ('stage/lapw_projector_functional.py', repo / 'src/python/lapw_projector_functional.py'),
    ('stage/lapw_projector_functional.py', dbin / 'lapw_projector_functional.py'),
    ('stage/tests/test_lapw_projector_functional.py', repo / 'src/python/tests/test_lapw_projector_functional.py'),
    ('stage/stress_response_roadmap_20261108.txt', repo / 'src/doc/stress/stress_response_roadmap_20261108.txt'),
]
manifest = []
patch = ''
for relative,destination in targets:
    source = root / relative
    item = dict(source=str(source), destination=str(destination),
                source_sha256=sha(source), before_sha256=sha(destination) if destination.exists() else None)
    manifest.append(item)
    if destination.is_relative_to(repo):
        name = str(destination.relative_to(repo))
        before = destination.read_text().splitlines(keepends=True) if destination.exists() else []
        patch += ''.join(difflib.unified_diff(before,source.read_text().splitlines(keepends=True),
            fromfile='a/'+name if destination.exists() else '/dev/null',tofile='b/'+name))
(root / 'projector_functional.patch').write_text(patch)
(root / 'installation_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

evidence = root / 'evidence'
evidence.mkdir(exist_ok=True)
files = ['REPORT.txt', 'projector_functional.patch', 'installation_manifest.json',
         'tests_repository.log', 'tests_dbin.log', 'SrVO3_material.log', 'FeSe_material.log',
         'SrVO3_full_mesh.log', 'FeSe_full_mesh.log', 'logdet_warning.json', 'logdet_warning_matrix.npy']
files += [case+'/'+name for case in ['SrVO3','FeSe'] for name in ['REPORT.json','projector_covectors.npz']]
files += ['full_bz/'+case+'/'+name for case in ['SrVO3','FeSe'] for name in ['REPORT.json','covectors.npz']]
for relative in files:
    dest = evidence / relative
    dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(root / relative,dest)
for relative in ['check_material.py','check_full_mesh.py','check_logdet_warning.py',
                 'run_tests.py','post_install_tests.py','install_checkpoint.py','record_installation.py',
                 'package_checkpoint.py','verify_checkpoint.py']:
    dest = evidence / 'reproduce' / relative
    dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(root / relative,dest)
for relative in ['lapw_projector_functional.py','tests/test_lapw_projector_functional.py',
                 'stress_response_roadmap_20261108.txt']:
    dest = evidence / 'stage' / relative
    dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(root / 'stage' / relative,dest)
hashes = {str(p.relative_to(evidence)):sha(p) for p in sorted(evidence.rglob('*'))
          if p.is_file() and p.name!='SHA256.json'}
(evidence / 'SHA256.json').write_text(json.dumps(hashes,indent=2)+'\n')
print('Prepared',len(manifest),'destinations and',len(hashes),'evidence files; nothing installed.')
