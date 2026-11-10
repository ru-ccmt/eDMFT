"""Build a reviewable patch/evidence package; never edit repository or dbin."""

from pathlib import Path
import difflib
import hashlib
import json
import shutil
import numpy as np

root = Path(__file__).resolve().parent
repo = Path("/Users/haule/projects/eDMFT_new")
dbin = Path("/Users/haule/dbin")
previous = Path("/private/tmp/edmft_strain_charge_20261109")
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
report = (root / "REPORT_body.txt").read_text()
report += (
    "225 tests pass using repository dependencies and installed dbin dependencies.\n"
)
report += "Ten new tests cover nonlinear scalar derivatives, direct/projector\nseparation, pure-basis covariance, zero-Sigma DFT, particle-number/gauge,\nsix-source forward/adjoint algebra, fixed-basis physical-potential integrals/transpose/scalar, and rejected missing inputs. Zero skips.\n\n"
report += "Case/mesh       lattice explicit   Legendre extrap error   relative error\n"
rows = []
for mesh in ["three_k", "full_bz"]:
    for case in ["SrVO3", "FeSe"]:
        r = json.loads((root / mesh / case / "REPORT.json").read_text())
        f = json.loads((root / mesh / case / "scalar_refinement.json").read_text())
        best = f["richardson"][0]["legendre_fixed_number"]
        report += f"{case:6s}/{mesh:7s} {r['explicit_lattice_derivative_Ry']:17.9f} {best['error_Ry']:23.6e} {best['relative_error']:16.6e}\n"
        rows.append((case, mesh, r, f))
report += "\nAll values above are Ry per hydrostatic parameter; not total stress.\n"
report += "\nFixed-basis physical-potential material checks (3 k):\n"
for case in ["SrVO3", "FeSe"]:
    p = json.loads((root / "three_k" / case / "physical_potential.json").read_text())
    report += f"  {case}: extrapolated scalar error {p['richardson']['error_Ry']:.6e} Ry; transpose relative error {p['transpose_relative_error']:.6e}; constant-potential RR-I norm {p['constant_potential_H_RR_minus_I_norm']:.6e}, QR norm {p['constant_potential_H_QR_norm']:.6e}.\n"
    report += f"  Its induced-potential matrix decomposition in Ry: {p['induced_potential_split_Ry']}. This remainder is not a total stress correction.\n"
    pairing = json.loads((root / "three_k" / case / "spatial_pairing.json").read_text())
    angular = json.loads((root / "three_k" / case / "angular_cutoff.json").read_text())
    report += f"  Raw spatial-pairing audit: {pairing['comparison']['screened_induced']}; band-minus-spatial valence number {pairing['band_minus_spatial_number']:.6e}.\n"
    report += "  Nonspherical same-coefficient cutoff/bank experiment:\n"
    for label, a in angular["rows"].items():
        report += f"    {label}: norm difference {a['difference_norm']:.6e}; induced difference {a['induced_difference_Ry']:.6e} Ry.\n"
    report += "  Spherical same-coefficient cutoff/bank experiment:\n"
    for label, a in angular["spherical_rows"].items():
        report += f"    {label}: norm difference {a['difference_norm']:.6e}; induced difference {a['induced_difference_Ry']:.6e} Ry.\n"
report += (
    "\nIndependent fixed-geometry potential-direction lattice derivatives (3 k):\n"
)
for case in ["SrVO3", "FeSe"]:
    p = json.loads((root / "three_k" / case / "potential_scalar.json").read_text())
    report += f"  {case}: analytic {p['analytic_Ry']:.12g} Ry; extrapolated error {p['richardson']['error_Ry']:.6e} Ry, relative {p['richardson']['relative_error']:.6e}.\n"
report += "\nCase/mesh       component     forward-adjoint error     residual bound\n"
for case, mesh, r, f in rows:
    for component in ["projector", "lattice"]:
        d = r["components"][component]
        report += f"{case:6s}/{mesh:7s} {component:10s} {d['forward_adjoint_error_Ry']:24.6e} {d['residual_bound_Ry']:18.6e}\n"
report += (
    "\nSeparated contractions in Ry (these are NOT independently additive stresses):\n"
)
for case, mesh, r, f in rows:
    report += f"\n{case}/{mesh}; central charge fixed-point residual norm {r['central']['central_fixed_point_residual_norm']:.8g}:\n"
    for component in ["projector", "direct", "lattice"]:
        d = r["components"][component]
        report += f"  {component:10s}: geometry {d['geometry_derivatives'][0]: .10f}, explicit potential {d['explicit_potential_derivatives'][0]: .10f}, induced {d['induced_derivatives'][0]: .10f}, sum {d['forward_Ry'][0]: .10f}\n"
    report += f"  Valence preparation {r['valence_preparation_seconds']:.2f} s; gradient preparation {r['gradient_seconds']:.2f} s.\n"
    for component in ["projector", "lattice"]:
        a = r["components"][component]["adjoint"]["adjoint"]
        report += f"  {component} adjoint: {a['actions']} actions, {a['seconds']:.2f} s, true residual {a['residual_norm']:.6e} (threshold {a['threshold']:.6e}).\n"
    if mesh == "three_k":
        a = f["richardson"][-1]["legendre_fixed_number"]
        report += f"  Smaller-pair extrapolated scalar error: {a['error_Ry']:.6e} Ry ({a['relative_error']:.6e} relative).\n"
report += "\nPreparation/solve times are serial adapter diagnostics on this machine, not\nproduction timing promises. Saved analytic q and sources were reused; their\npreparation costs are in the previous checkpoint. No native process was rerun.\n"
(root / "REPORT.txt").write_text(report)

road = (repo / "src/doc/stress/stress_response_roadmap_20261108.txt").read_text()
start = road.index("LATEST CHECKPOINT:")
end = road.index("THE THREE LEVELS THAT CURRENTLY COEXIST")
road = (
    road[:start]
    + """LATEST CHECKPOINT: functional_readout_20261109/REPORT.txt

Hydrostatic point 3 is implemented and validated: full native J/J^T,
explicit analytic band/potential/charge sources s_nu, and screened q_nu.
We are now in POINT 4. Its first scalar-readout part is implemented:
the regularized lattice scalar, split into direct density/H and normalized-
projector covectors, including explicit geometry and potential plus Kq.
Both scalar readouts pass actual forward/adjoint contractions on SrVO3 NREL
and FeSe RELA, with 3-k and full-mesh native fixtures. Independent native
strain derivatives of the lattice scalar pass after using a consistent
fixed-N Legendre comparison; fixed-mu finite steps have stronger curvature.
Pure moving-basis cancellation, zero-Sigma DFT and six-source algebra pass.
225 tests pass with repository and dbin dependencies. The detailed report
records absolute errors, residual bounds, timings and scope restrictions.

The LATTICE scalar is NOT the complete Klein functional. Its direct physical
potential contribution must be reconciled with density/potential, core and
Hartree/XC terms, then the actual existing Fortran stress baseline subtracted.
The fixed-basis physical-potential map P and P^T rho are now implemented and
validated on the small material meshes. Their matrix-defined subtraction
separates radial-basis and normalized-projector response, but matching it to
the spatial Klein functional and production baseline remains unfinished.
Do not append an isolated M_Sigma O' to the new orthonormal-gauge contraction:
its complete projector derivative already includes the self-energy metric
dependence. Compare the two gauges as complete sums before changing baseline.

Controlled potential-pairing tests identify a further discrete cancellation
issue: LNSMAX=4 in the nonspherical H matrix versus l<=8 in the density.
At fixed coefficients, matching cutoff and radial bank removes the observed
nonspherical matrix/spatial residual. Spherical H uses l<=10 versus density
l<=8; its smaller residual is also recorded. This is a contraction diagnostic,
not a reconverged LNSMAX/stress test. Retain the matrix-minus-spatial pairing
in the full functional audit rather than silently assuming it vanishes.

NEXT: complete this functional/baseline reconciliation, including native
quadrature/normalization and fixed-point residual controls; validate the full
frozen-Sigma Klein-tangent derivative on consistent converged centers.
Six physical strain layouts and production integration still follow.
The archived centers have recorded nonzero fixed-point residuals; this stage
does not claim a new physical total stress or pressure improvement.

"""
    + road[end:]
)
start = road.index("THE IMMEDIATE DELIVERABLE")
end = road.index("SOURCE AND EVIDENCE MAP", start)
road = (
    road[:start]
    + """THE IMMEDIATE DELIVERABLE

Point 4 is underway. Explicit+induced projector and complete lattice-scalar
readouts now pass native scalar and material forward/adjoint tests. Next
use the tested fixed-basis physical-potential contraction from the same H
assembly to reconcile its cancellation with the implemented spatial density,
core and Hartree/XC functional, and subtract the actual Fortran stress baseline.
This produces the complete correction gradient for one shared adjoint.
Then check the full frozen-Sigma Klein tangent, preserving the selectable
radial correction, before six-source and production publication work.

"""
    + road[end:]
)
road = road.replace(
    "   functional/baseline gradient has not yet been assembled.",
    "   lattice-scalar gradient and explicit+screened contractions now pass tests;\n   full functional/baseline subtraction remains.",
)
anchor = "src/python/lapw_complement_response.py\n"
road = road.replace(
    anchor,
    """src/python/{lapw_strain_functional,lapw_physical_potential}.py
  Separated direct and normalized-projector joint lattice covectors; complete
  explicit geometry/potential plus screened readout, forward and one-adjoint.
  Fixed-basis physical-potential map/transpose and its weighted-density readout.
  Not the full Klein functional; does not choose or subtract the baseline.
"""
    + anchor,
    1,
)
road = road.replace(
    "  strain_charge_20261109/REPORT.txt\n",
    "  strain_charge_20261109/REPORT.txt\n  functional_readout_20261109/REPORT.txt\n",
    1,
)
(root / "stage/stress_response_roadmap_20261108.txt").write_text(road)

install = (
    (previous / "install_checkpoint.py")
    .read_text()
    .replace(
        "analytic hydrostatic charge-source patch", "lattice functional-readout patch"
    )
    .replace("stress/strain_charge_20261109", "stress/functional_readout_20261109")
)
(root / "install_checkpoint.py").write_text(install)
post = (previous / "post_install_tests.py").read_text().replace("215", "225")
post = post.replace(
    "'test_lapw_strain_field',",
    "'test_lapw_strain_field', 'test_lapw_strain_functional', 'test_lapw_physical_potential',",
    1,
)
post += """
for name in ['lapw_strain_functional','lapw_physical_potential']:
    assert Path(sys.modules[name].__file__).resolve() == source / (name+'.py')
for name in ['test_lapw_strain_functional','test_lapw_physical_potential']:
    assert Path(sys.modules[name].__file__).resolve() == repo / 'src/python/tests' / (name+'.py')
print('Verified all checked helpers and tests load from installed locations.', flush=True)
"""
(root / "post_install_tests.py").write_text(post)

targets = [
    (
        root / "stage/lapw_physical_potential.py",
        repo / "src/python/lapw_physical_potential.py",
    ),
    (root / "stage/lapw_physical_potential.py", dbin / "lapw_physical_potential.py"),
    (
        root / "stage/tests/test_lapw_physical_potential.py",
        repo / "src/python/tests/test_lapw_physical_potential.py",
    ),
    (
        root / "stage/lapw_strain_functional.py",
        repo / "src/python/lapw_strain_functional.py",
    ),
    (root / "stage/lapw_strain_functional.py", dbin / "lapw_strain_functional.py"),
    (
        root / "stage/tests/test_lapw_strain_functional.py",
        repo / "src/python/tests/test_lapw_strain_functional.py",
    ),
    (
        root / "stage/stress_response_roadmap_20261108.txt",
        repo / "src/doc/stress/stress_response_roadmap_20261108.txt",
    ),
]
manifest, patch = [], ""
for source, dest in targets:
    manifest.append(
        dict(
            source=str(source),
            destination=str(dest),
            source_sha256=sha(source),
            before_sha256=sha(dest) if dest.exists() else None,
        )
    )
    if dest.is_relative_to(repo):
        relative = str(dest.relative_to(repo))
        patch += "".join(
            difflib.unified_diff(
                dest.read_text().splitlines(True) if dest.exists() else [],
                source.read_text().splitlines(True),
                fromfile="a/" + relative if dest.exists() else "/dev/null",
                tofile="b/" + relative,
            )
        )
(root / "functional_readout.patch").write_text(patch)
(root / "installation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

evidence = root / "evidence"
evidence.mkdir(exist_ok=True)


def copy(relative, target=None):
    dest = evidence / (relative if target is None else target)
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(root / relative, dest)


for name in [
    "REPORT.txt",
    "REPORT_body.txt",
    "functional_readout.patch",
    "installation_manifest.json",
    "tests_repository.log",
    "tests_dbin.log",
    "FeSe_transpose_profile.txt",
    "FeSe_profile.log",
]:
    copy(name)
for mesh in ["three_k", "full_bz"]:
    for case in ["SrVO3", "FeSe"]:
        prefix = Path(mesh) / case
        for name in ["REPORT.json", "scalar_refinement.json"]:
            copy(prefix / name)
        if mesh == "three_k":
            copy(prefix / "potential_scalar.json")
            copy(prefix / "physical_potential.json")
            copy(prefix / "physical_potential_covectors.npz")
            copy(f"{case}_physical.log")
            for name in [
                "spatial_pairing.json",
                "spatial_pairing_covectors.npz",
                "angular_cutoff.json",
            ]:
                copy(prefix / name)
            copy(f"{case}_angular.log")
        copy(f"{case}_{mesh}.log")
        arrays = np.load(root / prefix / "covectors.npz")
        np.savez_compressed(
            evidence / prefix / "covectors.npz",
            **{
                k: arrays[k]
                for k in arrays.files
                if not (mesh == "full_bz" and k.endswith("_joint"))
            },
        )
for p in (root / "stage").rglob("*"):
    if p.is_file() and "__pycache__" not in str(p):
        copy(p.relative_to(root))
for p in root.glob("*.py"):
    copy(p.name, "reproduce/" + p.name)
hashes = {
    str(p.relative_to(evidence)): sha(p)
    for p in sorted(evidence.rglob("*"))
    if p.is_file() and p.name != "SHA256.json"
}
(evidence / "SHA256.json").write_text(json.dumps(hashes, indent=2) + "\n")
print(
    "Staged",
    len(targets),
    "targets and",
    len(hashes),
    "evidence files; no source installed.",
)
