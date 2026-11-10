"""Connect the new Hartree K to the existing complete valence map.

This is a Hartree-only component of the unmixed charge Jacobian. It is NOT
the full Hartree/XC/core map and is not used to solve physical screening.
"""

from pathlib import Path
import sys, time, json
import numpy as np

root = Path(__file__).resolve().parent
sys.path[:0] = [str(root / "stage"), "/Users/haule/projects/eDMFT_new/src/python"]
from lapw_hartree_periodic import PeriodicHartreeResponse
from lapw_potential_response import NativeInterstitialPotential, NativePotentialBands
from lapw_complement_response import complement_electronic_operator
from stress_response_adjoint import PairedOperator
from wien_response_density import Density
from check_material import read_native

sys.path.insert(0, "/private/tmp/edmft_complement_response_20261108")
from build_reduced_material import build

case = sys.argv[1]
nr = {"FeSe": 66, "SrVO3": 43}[case]
source = (
    Path("/private/tmp/edmft_coupled_valence_response_20261108")
    / dict(FeSe="FeSe_b509_one", SrVO3="SrVO3_b703_one")[case]
    / "central"
)
t = time.perf_counter()
f = build(case, source, nr)
s = f["mat"].s
b = f["bands"]
lines = (source / (case + ".vns")).read_text().splitlines()
j = next(i for i, l in enumerate(lines) if "NUMBER OF PW" in l)
n = int(lines[j].split()[0])
reps = np.array(
    [[int(l[a : a + 5]) for a in [3, 8, 13]] for l in lines[j + 1 : j + 1 + n]]
)
irs = []
for record in f["records"]:
    npw = record.nbasis - record.nlo
    irs.append(
        NativeInterstitialPotential(
            record.coefficients[:npw, :, 0],
            record.labels[:, :npw].T,
            reps,
            s.timat.transpose(0, 2, 1),
            s.tau,
            complex_values=f["mat"].complex,
        )
    )
bands = NativePotentialBands(b, irs, retained=nr)
e, inputs, outputs = complement_electronic_operator(f["complement"], fixed_number=True)
shape = b.coordinates.band.shape


def rforward(x):
    rr, qr = outputs.unpack(x)
    rho = np.zeros(shape, complex)
    rho[:, :nr, :nr] = rr
    rho[:, nr:, :nr] = qr
    rho[:, :nr, nr:] = np.array(qr).swapaxes(1, 2).conj()
    return f["spatial"].action(rho)


def rreverse(x):
    rho, _, _ = f["spatial"].pullback(x)
    return outputs.pack(rho[:, :nr, :nr], rho[:, nr:, :nr])


r = PairedOperator(
    outputs.size,
    len(f["reference"].vector()),
    rforward,
    rreverse,
    "complete spatial valence",
)
valence = r @ e @ bands.operator + f["direct"].operator @ bands.spherical_selector
native = read_native(root / case / "0")
k = PeriodicHartreeResponse(
    f["reference"],
    s,
    native["shape"],
    ncon=native["ncon"],
    inversion=not f["mat"].complex,
)
select, coords = k.native_potential_output(
    tuple(b.slices),
    {a: part.layout.harmonics for a, (part, ks) in bands.ns.items()},
    reps,
    complex_interstitial=f["mat"].complex,
)
assert coords.nonspherical_shapes == bands.input.nonspherical_shapes
assert coords.size == bands.input.size
kh = select @ k.operator
combined = valence @ kh
dn = np.load(root / (case + "_direction.npz"))["density"]
assert len(dn) == k.size
result = {
    "scope": "fixed-N first-k valence map composed with Hartree only; XC/core absent",
    "preparation_seconds": time.perf_counter() - t,
    "input_size": k.size,
    "potential_size": coords.size,
}
t = time.perf_counter()
y = combined.matvec(dn)
result["forward_seconds"] = time.perf_counter() - t
t = time.perf_counter()
result["transpose_error"] = combined.check_transpose(samples=1, tolerance=1e-10)
result["paired_check_seconds"] = time.perf_counter() - t
np.savez(root / (case + "_hartree_valence.npz"), density=y, potential=kh.matvec(dn))
(root / (case + "_composition.json")).write_text(json.dumps(result, indent=2) + "\n")
print(case, json.dumps(result, indent=2), flush=True)
