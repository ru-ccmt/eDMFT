"""Prepare the reviewable core/XC checkpoint without touching repository/dbin."""
from pathlib import Path
import difflib,hashlib,json,shutil
root=Path(__file__).resolve().parent
repo=Path('/Users/haule/projects/eDMFT_new');dbin=Path('/Users/haule/dbin')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
road=(repo/'src/doc/stress/stress_response_roadmap_20261108.txt').read_text()
start=road.index('LATEST CHECKPOINT:');end=road.index('THE THREE LEVELS THAT CURRENTLY COEXIST')
road=road[:start]+'''LATEST CHECKPOINT: core_functional_20261109/REPORT.txt

We are in POINT 4, assembling the complete functional readout. Hydrostatic
point 3 already supplies analytic s, J/J^T, and screened q. The lattice scalar
and its forward/adjoint readouts were installed in functional_readout_20261109.
This checkpoint adds the native core eigenvalue-sum gradient, spatial pairing,
and actual discrete XC energy gradient. 239 regression tests pass against
repository and dbin dependencies; default local XC potentials are unchanged.

Core potential and explicit hydrostatic scalar derivatives pass native tests.
The XC audit finds a discrete cancellation defect: grad E_xc differs from
the spatial pairing with the separately implemented V_xc. The new analytic
energy gradient reproduces the PBE-minus-LDA native defect in FeSe, both for
a probe and saved full-mesh screened q. This is a derived scalar derivative,
not a pressure fit. It is not yet added to production stress.

At fixed geometry the eventual charge readout is
 c = K^T[B^T g_lattice + g_core - W n_input] + r_C,
 r_C = grad C_native + K^T W n_input.
The XC part r_xc = grad E_xc - W V_xc is now implemented and validated.
The remaining electrostatic/gauge residual and explicit geometric partners
must be derived before the complete covector and Fortran-baseline subtraction
are assembled. One combined adjoint then serves all physical strain sources.

Matrix/spatial pairing can also differ because H uses LNSMAX=4 and density
l<=8 in the archived cases, with different radial banks. Preserve that
measured difference; do not silently identify P^T rho with spatial n_input.
Neither raw nor normalized output density can replace input density away
from a fixed point. Core exterior continuation is retained in its readout.

NEXT: complete the native electrostatic density scalar (nuclear constants
and gauge included), explicit geometric density-functional terms, and actual
baseline subtraction. Validate the complete frozen-Sigma tangent on converged
centers, then finish six physical strain sources and production publication.
No new total pressure agreement is claimed. Wrappers still report
stress_applied=false; their production behavior has not changed here.

'''+road[end:]
start=road.index('THE IMMEDIATE DELIVERABLE');end=road.index('SOURCE AND EVIDENCE MAP',start)
road=road[:start]+'''THE IMMEDIATE DELIVERABLE

Complete point 4 by differentiating the remaining native electrostatic scalar
and its geometry, combining it with the tested lattice/core/XC gradients,
then subtracting the existing Fortran baseline. Keep the finite matrix/spatial
pairing residual and charge fixed-point controls. Only the complete readout
should drive one shared adjoint and subsequent total-stress comparisons.

'''+road[end:]
anchor='src/python/lapw_complement_response.py\n'
road=road.replace(anchor,'''src/python/lapw_functional_cancellation.py
  Native spatial W pairing, matched-core eigenvalue gradient with exterior
  continuation, and exact discrete XC energy/density-term gradients. These
  are fixed-geometry components, not a published total stress correction.
'''+anchor,1)
road=road.replace('  functional_readout_20261109/REPORT.txt\n','  functional_readout_20261109/REPORT.txt\n  core_functional_20261109/REPORT.txt\n')
(root/'stage/stress_response_roadmap_20261108.txt').write_text(road)
changes=[]
for name in ['lapw_functional_cancellation.py','lapw_xc_response.py']:
    before=repo/'src/python'/name
    installed=dbin/name
    if before.exists() and installed.exists() and sha(before)!=sha(installed):
        raise RuntimeError('Repository/dbin preimages differ; inspect '+name)
    changes.append((root/'stage'/name,before))
    changes.append((root/'stage'/name,installed))
for p in sorted((root/'stage/tests').glob('test_*.py')):
    changes.append((p,repo/'src/python/tests'/p.name))
changes.append((root/'stage/stress_response_roadmap_20261108.txt',repo/'src/doc/stress/stress_response_roadmap_20261108.txt'))
manifest=[];patch=[]
for source,dest in changes:
    manifest.append(dict(source=str(source),destination=str(dest),source_sha256=sha(source),before_sha256=sha(dest) if dest.exists() else None))
    if dest.is_relative_to(repo):
        rel=dest.relative_to(repo)
        patch.extend(difflib.unified_diff(dest.read_text().splitlines(True) if dest.exists() else [],source.read_text().splitlines(True),fromfile='a/'+str(rel),tofile='b/'+str(rel)))
(root/'core_functional.patch').write_text(''.join(patch))
(root/'installation_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
evidence=root/'evidence';evidence.mkdir(exist_ok=True)
for name in ['REPORT.txt','tests_repository.log','tests_dbin.log','unit.log','potential_regression.json','material.log','core_strain.log','xc_scalar.log','density_energy_audit.log','core_functional.patch','installation_manifest.json']:
    shutil.copy2(root/name,evidence/name)
for case in ['SrVO3','FeSe']:
    dst=evidence/case;dst.mkdir(exist_ok=True)
    for name in ['REPORT.json','covectors.npz','core_strain.json','xc_scalar.json','xc_energy_covectors.npz','density_energy_audit.json']:
        shutil.copy2(root/case/name,dst/name)
    for base in ['native_energy','density_energy','xc_energy_central']:
        for source in (root/case/base).rglob(case+'.scf0'):
            if source.stat().st_size==0:continue # failed unsupported XC_NONE probe
            if '/Hartree/' in str(source):continue
            rel=source.relative_to(root/case)
            target=dst/'native'/rel;target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(source,target)
            shutil.copy2(source.parent/(case+'.in0'),target.parent/(case+'.in0'))
scripts=['check_material.py','check_native_energy.py','audit_density_energy.py','check_xc_scalar.py','check_core_strain.py','run_tests.py','package_checkpoint.py','install_checkpoint.py','verify_checkpoint.py']
for name in scripts:
    dst=evidence/'reproduce'/name;dst.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(root/name,dst)
for src in (root/'stage').rglob('*'):
    if src.is_file() and src.suffix in {'.py','.txt'}:
        dst=evidence/'stage'/src.relative_to(root/'stage');dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src,dst)
hashes={str(p.relative_to(evidence)):sha(p) for p in sorted(evidence.rglob('*')) if p.is_file() and p.name!='SHA256.json'}
(evidence/'SHA256.json').write_text(json.dumps(hashes,indent=2)+'\n')
print('Prepared',len(manifest),'installation targets and',len(hashes),'evidence files; repository/dbin unchanged.')
