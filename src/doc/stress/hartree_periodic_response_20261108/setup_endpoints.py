from pathlib import Path
import sys, os, shutil, subprocess
import numpy as np

sys.path.insert(0, "/Users/haule/projects/eDMFT_new/src/python")
from wien_response_density import Density

root = Path(__file__).resolve().parent
case = sys.argv[1]
original = (
    Path("/Users/haule/projects/materials")
    / dict(FeSe="FeSe_again/FeSe_stress/FeSe_RELA", SrVO3="SrVO3/V_plus_force")[case]
)
fixture = (
    Path("/private/tmp/edmft_coupled_valence_response_20261108")
    / dict(FeSe="FeSe_b509_one", SrVO3="SrVO3_b703_one")[case]
    / "central"
)
central = Density(original / (case + ".clmsum"), original / (case + ".struct"))
induced = Density(fixture / (case + ".clmval"), fixture / (case + ".struct"))
induced = induced.from_vector(
    np.load(
        Path("/private/tmp/edmft_valence_charge_checks_20261108")
        / (case + "_joint_fixedN.npz")
    )["density"]
)
assert (original / (case + ".struct")).read_bytes() == (
    fixture / (case + ".struct")
).read_bytes()
dn = central.zeros().add(induced)
env = os.environ.copy()
env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
for amplitude in [0, 0.01, -0.01, 0.005, -0.005]:
    folder = root / case / str(amplitude)
    folder.mkdir(parents=True, exist_ok=True)
    for suffix in ["struct", "in0", "inm"]:
        shutil.copy2(original / (case + "." + suffix), folder / (case + "." + suffix))
    text = (folder / (case + ".in0")).read_text()
    text = "TOT " + text[4:]
    (folder / (case + ".in0")).write_text(text)
    # Relative .def outputs and a private copy of the read-only XC table.
    text = (
        (original / "lapw0.def")
        .read_text()
        .replace("/Users/haule/wien2k_24/SRC_lapw0/xc_funcs.h", "xc_funcs.h")
    )
    (folder / "lapw0.def").write_text(text)
    shutil.copy2("/Users/haule/wien2k_24/SRC_lapw0/xc_funcs.h", folder / "xc_funcs.h")
    central.copy().add(dn, amplitude).write(folder / (case + ".clmsum"))
    with (folder / "lapw0.log").open("w") as log:
        p = subprocess.run(
            [str(root / "native/SRC_lapw0/lapw0"), "lapw0.def"],
            cwd=folder,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    p.check_returncode()
    if "HARTREE REFERENCE COMPLETE" not in (folder / "lapw0.log").read_text():
        raise RuntimeError(str(folder))
    print(case, amplitude, "complete", flush=True)
np.savez(root / (case + "_direction.npz"), density=dn.vector())
