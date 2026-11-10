"""Check actual material results and guarded patch preimages."""

from pathlib import Path
import ast
import hashlib
import json

root = Path(__file__).resolve().parent
if root.name == "reproduce":
    root = root.parent
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
for label in ["repository", "dbin"]:
    log = (root / f"tests_{label}.log").read_text()
    assert "Ran 225 tests" in log and "\nOK\n" in log
    assert "Verified staged readout helper/test and installed dependencies." in log
for mesh, nks in [("three_k", [3, 3]), ("full_bz", [75, 140])]:
    for case, nk in zip(["SrVO3", "FeSe"], nks):
        r = json.loads((root / mesh / case / "REPORT.json").read_text())
        f = json.loads((root / mesh / case / "scalar_refinement.json").read_text())
        assert r["kpoints"] == nk
        for rr in f["richardson"]:
            error = rr["legendre_fixed_number"]
            assert abs(error["error_Ry"]) < 3e-6, (case, mesh, rr)
            assert error["relative_error"] < 2e-8, (case, mesh, rr)
        for component in ["projector", "lattice"]:
            c = r["components"][component]
            assert not c["full_klein_functional"] and not c["stress_applied"]
            assert (
                c["forward_adjoint_error_Ry"] <= 1.001 * c["residual_bound_Ry"] + 2e-10
            )
            a = c["adjoint"]["adjoint"]
            assert a["converged"] and a["residual_norm"] <= a["threshold"]
        if mesh == "three_k":
            p = json.loads((root / mesh / case / "potential_scalar.json").read_text())
            assert p["richardson"]["relative_error"] < 3e-8
            p = json.loads((root / mesh / case / "physical_potential.json").read_text())
            assert p["richardson"]["relative_error"] < 3e-8
            assert p["transpose_relative_error"] < 1e-11
            assert p["constant_potential_H_RR_minus_I_norm"] < 1e-9
            assert p["constant_potential_H_QR_norm"] < 1e-9
            assert abs(p["constant_potential_number_derivative"] - p["number"]) < 1e-9
            a = json.loads((root / mesh / case / "angular_cutoff.json").read_text())
            assert a["rows"]["matrix_l4"]["versus_native_matrix_norm"] < 1e-12
            assert a["rows"]["density_l8"]["difference_norm"] < 1e-11
            assert (
                a["rows"]["matrix_l8"]["difference_norm"]
                < 1e-5 * a["rows"]["matrix_l4"]["difference_norm"]
            )
            # The smaller spherical FeSe residual remains documented; do not
            # demand or claim roundoff agreement for that separate sector.
            assert (
                a["spherical_rows"]["matrix_l10"]["versus_native_matrix_norm"] < 1e-12
            )
            assert a["spherical_rows"]["density_l8"]["difference_norm"] < 1e-8
for path in (root / "stage").rglob("*.py"):
    ast.parse(path.read_text(), filename=str(path))
for item in json.loads((root / "installation_manifest.json").read_text()):
    source, dest = Path(item["source"]), Path(item["destination"])
    assert sha(source) == item["source_sha256"]
    assert (sha(dest) if dest.exists() else None) == item["before_sha256"]
print(
    "225-test suites, native scalar derivatives, material adjoints and unchanged targets verified."
)
