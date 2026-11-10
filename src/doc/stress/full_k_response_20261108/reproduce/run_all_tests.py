from pathlib import Path
import sys, re, unittest, os

root = Path(__file__).resolve().parent
source = os.environ.get(
    "EDMFT_TEST_SOURCE", "/Users/haule/projects/eDMFT_new/src/python"
)
sys.path[:0] = (
    [str(root / "stage")] if os.environ.get("EDMFT_TEST_INSTALLED") != "1" else []
) + [
    str(root / "stage/tests"),
    source,
    "/Users/haule/projects/eDMFT_new/src/python/tests",
]
text = Path(
    "/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/post_install_repository.log"
).read_text()
modules = sorted(
    set(re.findall(r"\((test_\w+)\.", text))
    | {
        "test_lapw_compact_response",
        "test_lapw_valence_state",
        "test_lapw_charge_response",
        "test_lapw_potential_core_state",
        "test_lapw_core_response",
        "test_lapw_xc_response",
    }
)
print(
    "Source:",
    source,
    "Installed:",
    os.environ.get("EDMFT_TEST_INSTALLED", "0"),
    flush=True,
)
print("Modules:", ", ".join(modules), flush=True)
suite = unittest.defaultTestLoader.loadTestsFromNames(modules)
result = unittest.TextTestRunner(verbosity=2).run(suite)
raise SystemExit(not result.wasSuccessful() or bool(result.skipped))
