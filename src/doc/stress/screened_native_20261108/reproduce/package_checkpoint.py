"""Package tested staged sources and compact, explicitly scoped evidence."""
from pathlib import Path
import difflib, hashlib, json, shutil
import numpy as np

root=Path(__file__).resolve().parent
repo=Path('/Users/haule/projects/eDMFT_new')
dbin=Path('/Users/haule/dbin')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
reports={}
for mesh in ['three_k','full_bz']:
    for case in ['SrVO3','FeSe']:
        p=root/mesh/case
        r=json.loads((p/'screening_report.json').read_text())
        assert r['forward']['converged'] and len(r['native_residuals'])==2
        assert max(v['relative_to_source'] for v in r['native_residuals'].values())<2e-6
        r['source']=json.loads((p/'source_report.json').read_text())
        if mesh=='three_k':
            r['native']=json.loads((p/'native_screening_report.json').read_text())
            assert r['native']['converged'] and r['native']['density_relative_difference']<2e-6
            r['adjoint_bound']=json.loads((p/'adjoint_residual_bound.json').read_text())
            assert r['adjoint_bound']['passed']
        reports[mesh,case]=r
for label in ['repository','dbin']:
    t=(root/f'tests_{label}.log').read_text()
    assert 'Ran 188 tests' in t and '\nOK\n' in t

lines=['''PHASE 3: COMPLETE NATIVE SCREENED CHARGE CONNECTION
8 November 2026

STATUS: tested staged patch; not yet installed.
No production stress tensor or new GPa agreement is claimed.

WHAT IS NEW

NativeScreenedState connects the independently prepared potential/core and
full-k valence states to ChargeResponseChain and the existing paired GMRES
solver. It checks density ordering/grids, core and valence reconstruction,
common structure/spherical-potential fingerprints, explicit number/normalization policy
and native potential coordinates. One electronic operator handles the global
k sum, normalized projector X, complement response and chemical potential.
The preparation is reused in forward and adjoint solves. No per-k screening
or extra native executable launch is hidden inside a response action.
Full execution-snapshot provenance remains the caller/loader's responsibility;
this adapter is not an automatic exporter or a complete execution-cache loader.

The only new production helper proposed is lapw_screened_state.py. There are
five new tests and a roadmap update. No existing numerical kernel, Fortran
executable, dmft2 wrapper, mixer or user material directory is changed.

EQUATIONS AND COORDINATES

q: total density variation, in the real central Density.vector coordinates
   (MT radial CLMSUM coefficients and complex Fourier star sums).
v: prepared potential/core field coordinates; MT entries include r*delta V,
   IR entries use internal REAN4 coefficients. The selector conjugates once
   at the complex VNS file boundary for the native band assembly.
K: Hartree+XC derivative, q -> v.
B: joint H/O/raw-projector derivative, including the separate WIEN24 H/O and
   ATPARN projector radial policies, and the potential-coordinate selector.
E: frozen-local-Sigma/DC electronic derivative, with fixed-index window motion,
   symmetrized normalization, overlap metric, empty complement and fixed-N mu.
R: fixed-basis density reconstruction; D: direct moving-basis density response.
C: spherical core response; N: MIXER output-charge normalization derivative,
   evaluated at the UNNORMALIZED total central valence+core output.

  J = N [R E B + D + C] K,
  (I-J) q_eta = b_eta.

b_eta contains all explicit strain/coordinate/normalization contributions.
For this validation it comes from independent native +/-hydrostatic one-step
charge evaluations, pulled back into the central coordinates. It is not yet
an analytic implementation of all six geometric sources.

The solver returns q_eta and v_induced=K q_eta. It deliberately does not call
v_induced the total potential. The latter requires the separately determined
explicit affine potential: v_total=v_affine+K q_eta, counted exactly once.
The existing chain.solve_perturbation accepts all those analytic sources when
they are available and reconstructs total band H/O/raw-projector inputs.

For ONE scalar observable with density covector c_q and potential covector c_v,
  c = c_q + K^T c_v,
  (I-J)^T lambda = c,
  dA/deta_nu = a_nu + lambda^T b_nu.
Transposes are exact real-coordinate transposes, including every quadrature,
normalization and Fourier convention. J is not assumed symmetric. One adjoint
can serve multiple sources of that SAME scalar observable. The physical
functional/stress covector has NOT yet been connected by this patch.

INDEPENDENT NATIVE REFERENCE

Each native evaluation runs LAPW0 -> LAPW1 -> dmft2 -> LCORE, normalizes total
charge and applies the same material-coordinate pullback/normalization.
Sigma/DC, RMT, local frames, radial energies, G labels, stored band counts,
k list and DMFT window indices stay fixed. SOlapm is recomputed. Endpoint mu
is obtained from an independent nonlinear resolvent trace/root, without the
analytic dmu. Native completion/error files and window boundary gaps are
checked. Two probe steps check Jq after the analytic solve; these native
actions are not used by GMRES. On the three-k meshes an independent cold-start
Anderson solve uses only native centered-difference actions.
The endpoint electronic-number target is the central quadrature N0, so the
derivative is evaluated at the identical central mu with dN=0. It is not
rounded to the nominal integer: full-mesh N0 is 41.0018262968246 for SrVO3
and 64.0011746495842 for FeSe. Outer total-charge normalization still targets
85 and 120 electrons respectively. A final material calculation must first
prepare a central state converged with its intended physical number policy.

Hydrostatic source steps are 0.00025 and 0.000125 in lattice scale, with actual
formatted structure lengths used in the denominator. FeSe's native fixed-width
cell rounding gives the near-hydrostatic direction
  (0.99987401, 0.99987401, 1.00025198)
at both steps; SrVO3's direction is (1,1,1). This small geometric qualification
is retained rather than declaring an exact analytic hydrostatic source.
Some FeSe strained outputs append zero high-G density stars. The reference
aligns by central labels and rejects any nonzero discarded channel.
The existing density transport uses cubic splines and a small extrapolation
at the fixed-RMT edge. Step halving tests this reference numerically; it does
not replace an analytic derivation of all moving-boundary geometry sources.

FULL CONFIGURED K-MESH RESULTS

Relative residual below means ||b+Jq-q||/||b||; it is not a pressure error.
The native residual uses Jq=[F(n0+t*q)-F(n0-t*q)]/(2t).
Material Nk stored bands analytic residual native(t=1e-5) native(t=5e-6)
''']
for case in ['SrVO3','FeSe']:
    r=reports['full_bz',case];f=r['forward']
    lines.append(f"{case:7s} {f['kpoints']:3d} {f['retained_bands']:3d} {f['residual_norm']/r['source_norm']:.9e} {r['native_residuals']['1e-05']['relative_to_source']:.9e} {r['native_residuals']['5e-06']['relative_to_source']:.9e}\n")
lines.append('''
SrVO3's smaller-step native probe is noisier; monotonic quadratic improvement
is not asserted. The analytic solve uses rtol=2e-8, atol=1e-9 and verifies its
true final residual. Failure raises and produces no screened-density result.

SMALL-MESH INDEPENDENT SOLVE AND ADJOINT CONTROLS

Material native-vs-analytic density relative error; adjoint contraction error
''')
for case in ['SrVO3','FeSe']:
    r=reports['three_k',case]
    lines.append(f"{case:7s} {r['native']['density_relative_difference']:.9e}; {r['adjoint']['absolute_difference']:.9e}\n")
lines.append('''
The native cold-start stopping threshold is 5e-4 in density-coefficient norm.
The tested adjoint observable is the integrated MT electron number on atom 1,
including multiplicity, NOT the Klein functional or stress. Its forward and
adjoint contractions agree within the true-residual error bound
  ||lambda|| ||r_forward|| + ||q|| ||r_adjoint||.
The same lambda also agrees with a separate forward solve for the second
native source step. All component transposes were validated previously.

TIMINGS AND CENTRAL-STATE LIMITATIONS

Material mesh preparation[s] forward actions/seconds peakRSS[GB] fixed-point residual norm
''')
for mesh in ['three_k','full_bz']:
    for case in ['SrVO3','FeSe']:
        r=reports[mesh,case];f=r['forward']
        lines.append(f"{case:7s} {mesh:7s} {r['preparation_seconds']:.3f} {f['actions']:3d}/{f['seconds']:.3f} {r['peak_process_rss_bytes']/1e9:.3f} {f['central_fixed_point_residual_norm']:.9e}\n")
lines.append('''
These are single-thread local measurements; small tests/native preparations
sometimes overlapped. They are not a controlled production speed benchmark.
Preparation is shared by subsequent directions and the adjoint, but resident
k data and repeated spatial/band actions remain costly. MPI/streaming and
further action optimization remain important before production use.

The central diagnostic states were NOT newly converged. The three-k grids
deliberately differ from the input-density mesh. The full-mesh test variants
retain a central charge mismatch. The adapter reports the norm of
F(n0)-n0 explicitly: passing a local linearized equation does not establish
the derivative of a converged physical free energy. Final stress validation
requires matched converged central states and the physical functional covector.
Relative to the central density-vector norm, the full-mesh fixed-point residual
is 4.724e-5 for SrVO3 and 7.930e-4 for FeSe; these are diagnostic snapshots,
not new statements about convergence in the user's live material directories.
The test fixtures use FeSe RELA/PBE and SrVO3 NREL/LDA, RKmax=8, fixed-index
projector 5, nonmagnetic/no SOC, primitive orthogonal geometry and fixed RMT.
The integrated production wrapper has NOT switched to this new backend.

SOURCE STEP CONTROL

Material mesh relative source step difference
''')
for mesh in ['three_k','full_bz']:
    for case in ['SrVO3','FeSe']:
        lines.append(f"{case:7s} {mesh:7s} {reports[mesh,case]['source']['step_relative_difference']:.9e}\n")
lines.append('''
These finite-source differences are larger than the screened solver residual.
They are tracked separately and are not hidden in the electronic tolerance.

TESTS AND NEXT STEP

188 tests pass, zero skips, with repository and dbin dependency overlays.
New controls cover an independent dense nonsymmetric screened solve, one
observable/two-source adjoint, layout and stale-state rejection, core closure,
explicit ensemble choice and failure without a returned result.

Point 2 remains complete for the tested scope. Point 3 is now underway with
the native screened charge connection validated. Next: connect the physical
functional/potential covector, combine the explicit affine potential exactly
once, validate the total potential/stress contraction, then supply six
analytic strain sources and integrate/publish the tensor in the executor.
General geometry, magnetic/SOC support and production MPI remain pending.
The force window/normalization audit remains deferred until strain is validated.
''')
(root/'REPORT.txt').write_text(''.join(lines))

mapping=[(root/'stage/lapw_screened_state.py',repo/'src/python/lapw_screened_state.py'),
         (root/'stage/lapw_screened_state.py',dbin/'lapw_screened_state.py'),
         (root/'stage/tests/test_lapw_screened_state.py',repo/'src/python/tests/test_lapw_screened_state.py'),
         (root/'roadmap.txt',repo/'src/doc/stress/stress_response_roadmap_20261108.txt')]
manifest=[];diff=[]
for src,dst in mapping:
    manifest.append(dict(source=str(src),destination=str(dst),source_sha256=sha(src),before_sha256=sha(dst) if dst.exists() else None))
    if dst.is_relative_to(repo):
        rel=str(dst.relative_to(repo));old=dst.read_text().splitlines(True) if dst.exists() else []
        diff.extend(difflib.unified_diff(old,src.read_text().splitlines(True),fromfile='a/'+rel if dst.exists() else '/dev/null',tofile='b/'+rel))
(root/'screened_native.patch').write_text(''.join(diff))
(root/'installation_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

evidence=root/'evidence';evidence.mkdir(exist_ok=True)
for name in ['REPORT.txt','NEXT_STEP.txt','screened_native.patch','installation_manifest.json','roadmap.txt',
             'tests_repository.log','tests_dbin.log']:
    shutil.copy2(root/name,evidence/name)
for name in ['native_reference.py','native_screening.py','check_screening.py','run_tests.py',
             'post_install_tests.py','package_checkpoint.py']:
    dest=evidence/'reproduce'/name;dest.parent.mkdir(exist_ok=True);shutil.copy2(root/name,dest)
shutil.copytree(root/'stage',evidence/'staged_source',dirs_exist_ok=True)
for mesh in ['three_k','full_bz']:
    for case in ['SrVO3','FeSe']:
        src=root/mesh/case;dst=evidence/mesh/case;dst.mkdir(parents=True,exist_ok=True)
        for p in src.iterdir():
            if p.is_file() and p.suffix in ['.json','.npy']: shutil.copy2(p,dst/p.name)
            elif p.is_dir() and p.name.startswith(('central','source_','verify_')):
                d=dst/p.name;d.mkdir(exist_ok=True)
                for n in ['done.json','normalized.npy']:
                    if (p/n).exists():shutil.copy2(p/n,d/n)
        for suffix in ['source','screening']:
            log=root/f'{case}_{suffix}_{mesh}.log'
            if log.exists():shutil.copy2(log,dst/log.name)
        if mesh=='three_k':shutil.copy2(root/f'{case}_native_screening.log',dst/'native_screening.log')
inputs={}
for case in ['SrVO3','FeSe']:
    src=Path('/private/tmp/edmft_full_k_response_20261108/full_bz')/case/'charge/0'
    for p in src.iterdir():
        if p.is_file() and (p.name.startswith(case+'.') or p.name.startswith(('response.','native.','sig.inp')) or p.name in ['projectorw.dat','EF.dat']):
            inputs[str(p)]=sha(p)
    p=repo/'src/doc/stress/charge_closure_20261108'/f'{case}_potential_core.npz';inputs[str(p)]=sha(p)
(evidence/'reference_fingerprints.json').write_text(json.dumps(inputs,indent=2)+'\n')
hashes={str(p.relative_to(evidence)):sha(p) for p in sorted(evidence.rglob('*')) if p.is_file() and p.name!='SHA256.json'}
(evidence/'SHA256.json').write_text(json.dumps(hashes,indent=2)+'\n')
print('Prepared patch:',len(manifest),'install targets;',len(hashes),'evidence files')
