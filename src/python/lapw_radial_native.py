"""In-memory access to the native discrete scalar radial response.

Build liblapw_response.so in src/dmft2/atpar and supply that path, or install
it alongside this helper. No external executable, files, energy search or
finite differences occur in action/pullback. One context is shared across
all k points and equivalent sites of an atom/l channel.

These are derivatives of OUTWIN/ATPARN's discrete scalar NREL/RELA algorithm.
DIRACOUT p1/2 is unsupported. Input potential is r*V_spherical in Ry and
linearization/LO energies are Ry. The grid, RMT and nuclear charge are fixed.
"""

import ctypes as ct
from pathlib import Path
import threading
import weakref
import numpy as np


# The native handle registry is process-local and not thread-safe. Serialize
# calls, including creation/destruction, but do not serialize MPI processes.
_LOCK = threading.RLock()
_LIBRARIES = {}
_ARRAY = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")


def _real(value, shape, name):
    a = np.asarray(value)
    if a.shape != shape or np.iscomplexobj(a) or not np.all(np.isfinite(a)):
        raise ValueError("invalid " + name)
    return np.array(a, dtype=np.float64, order="F", copy=True)


def _status(code):
    if code:
        raise RuntimeError(
            {
                1: "invalid or closed native radial context",
                2: "invalid native radial input",
                3: "native radial context registry is full",
                4: "nonfinite native radial result",
            }.get(code, f"native radial error {code}")
        )


def _library(path):
    if path is None:
        directory = Path(__file__).resolve().parent
        choices = [
            directory / "liblapw_response.so",
            directory.parent / "dmft2/atpar/liblapw_response.so",
        ]
        path = next((p for p in choices if p.is_file()), None)
        if path is None:
            raise FileNotFoundError(
                "build/install liblapw_response.so or pass library=path"
            )
    path = str(Path(path).resolve())
    with _LOCK:
        if path in _LIBRARIES:
            return _LIBRARIES[path]
        lib = ct.CDLL(path)
        i, d, a = ct.c_int, ct.c_double, _ARRAY
        lib.edmft_radial_api_version.argtypes = []
        lib.edmft_radial_api_version.restype = i
        if lib.edmft_radial_api_version() != 2:
            raise RuntimeError(
                "incompatible native radial API; rebuild/install the response library"
            )
        signatures = {
            "api_version": [],
            "create": [i, i, i, i, i, d, d, d, d, d, a, a, ct.POINTER(i)],
            "destroy": [i],
            "value": [i, a, a],
            "metrics": [i, a, a, a],
            "action": [i, a, d, a, a, a],
            "pullback": [i, a, a, a, ct.POINTER(d), a],
        }
        if hasattr(lib, "edmft_radial_create_policy"):
            signatures["create_policy"] = signatures["create"][:-1] + [i, ct.POINTER(i)]
        for name, args in signatures.items():
            func = getattr(lib, "edmft_radial_" + name)
            func.argtypes, func.restype = args, i
        _LIBRARIES[path] = lib
        return lib


def _release(lib, handle):
    with _LOCK:
        _status(lib.edmft_radial_destroy(handle))


class NativeRadialContext:
    """Prepared scalar channel, usable as a context manager or with close().

    primitive[r,component,slot] has components (A,B,AP,BP) in the stored
    OUTWIN convention and slots (u,udot,ordinary LOs). lo_coefficients maps
    these primitives to ordinary LOs; APW+lo's first LO uses u and udot.
    quadrature="atparn" preserves dmft2's RINT13 rule. "wien2k24" uses
    native LAPW1's default CHARGE rule (extended Simpson, not a cubic spline).
    This policy affects radial normalization, udot, ABC and all derivatives.
    Returned values are read-only. Forward and reverse use real Euclidean
    pairings; the native radial quadratures are already inside the derivative.
    """

    def __init__(
        self,
        potential,
        *,
        r0,
        dx,
        rmt,
        energy,
        lo_energies,
        ell,
        charge,
        relativistic,
        lapw,
        dirac_lo=False,
        quadrature="atparn",
        library=None,
    ):
        if quadrature not in ("atparn", "wien2k24"):
            raise ValueError("quadrature must be atparn or wien2k24")
        self.quadrature = quadrature
        if np.any(dirac_lo):
            raise NotImplementedError(
                "DIRACOUT p1/2 radial response is not implemented"
            )
        v, elo = np.asarray(potential), np.asarray(lo_energies)
        if v.ndim != 1 or len(v) < 6 or elo.ndim != 1:
            raise ValueError("invalid radial dimensions")
        if not isinstance(ell, (int, np.integer)) or ell < 0:
            raise ValueError("ell must be a nonnegative integer")
        if not isinstance(relativistic, (bool, np.bool_)) or not isinstance(
            lapw, (bool, np.bool_)
        ):
            raise ValueError("relativistic and lapw must be booleans")
        scalars = _real([r0, dx, rmt, energy, charge], (5,), "radial geometry/energy")
        r0, dx, rmt, energy, charge = scalars
        if r0 <= 0 or dx <= 0 or rmt <= r0 or charge <= 0:
            raise ValueError("invalid logarithmic grid or nuclear charge")
        if abs(np.log(rmt / r0) - dx * (len(v) - 1)) > 1e-10:
            raise ValueError("RMT is inconsistent with the logarithmic grid")
        v = _real(v, v.shape, "radial potential")
        elo = _real(elo, elo.shape, "LO energies")
        self.n, self.nlo, self.ell = len(v), len(elo), int(ell)
        self.lapw, self.relativistic = bool(lapw), bool(relativistic)
        self.rmt = float(rmt)
        self.energy = float(energy)
        self.lo_energies = elo.copy()
        self.lo_energies.setflags(write=False)
        self._lib = _library(library)
        if quadrature == "wien2k24" and not hasattr(
            self._lib, "edmft_radial_create_policy"
        ):
            raise RuntimeError(
                "WIEN2k24 quadrature requires the updated response library"
            )
        create = self._lib.edmft_radial_create
        policy_args = []
        if quadrature == "wien2k24":
            create = self._lib.edmft_radial_create_policy
            policy_args = [1]
        handle = ct.c_int()
        with _LOCK:
            _status(
                create(
                    self.n,
                    self.nlo,
                    int(relativistic),
                    ell,
                    int(lapw),
                    r0,
                    dx,
                    rmt,
                    energy,
                    charge,
                    v,
                    elo,
                    *policy_args,
                    ct.byref(handle),
                )
            )
        self._handle = handle.value
        self._finalizer = weakref.finalize(self, _release, self._lib, self._handle)
        self.primitive = np.empty((self.n, 4, 2 + self.nlo), order="F")
        self.lo_coefficients = np.empty((2 + self.nlo, self.nlo), order="F")
        self.radius = np.empty(self.n)
        self.overlap_metric = np.empty((self.n, 4), order="F")
        self.projector_metric = np.empty((self.n, 4), order="F")
        try:
            with _LOCK:
                _status(
                    self._lib.edmft_radial_value(
                        self._handle, self.primitive, self.lo_coefficients
                    )
                )
                _status(
                    self._lib.edmft_radial_metrics(
                        self._handle,
                        self.radius,
                        self.overlap_metric,
                        self.projector_metric,
                    )
                )
            if not all(
                np.all(np.isfinite(a)) for a in (self.primitive, self.lo_coefficients)
            ):
                raise RuntimeError("nonfinite central radial solution")
        except Exception:
            self.close()
            raise
        self.primitive.setflags(write=False)
        self.lo_coefficients.setflags(write=False)
        for value in (self.radius, self.overlap_metric, self.projector_metric):
            value.setflags(write=False)

    def _check_open(self):
        if not self._finalizer.alive:
            raise RuntimeError("native radial context is closed")

    def close(self):
        """Release the prepared native context; repeated close is harmless."""
        with _LOCK:
            self._finalizer()

    def __enter__(self):
        self._check_open()
        return self

    def __exit__(self, *args):
        self.close()

    def action(self, potential_dot, energy_dot=0.0, lo_energies_dot=None):
        """Return primitive and ABC tangents; default keeps all energies fixed."""
        dv = _real(potential_dot, (self.n,), "radial potential derivative")
        de = float(_real(energy_dot, (), "energy derivative"))
        delo = _real(
            np.zeros(self.nlo) if lo_energies_dot is None else lo_energies_dot,
            (self.nlo,),
            "LO energy derivative",
        )
        dp, dm = np.empty_like(self.primitive), np.empty_like(self.lo_coefficients)
        with _LOCK:
            self._check_open()
            _status(self._lib.edmft_radial_action(self._handle, dv, de, delo, dp, dm))
        if not all(np.all(np.isfinite(a)) for a in (dp, dm)):
            raise RuntimeError("nonfinite radial tangent")
        return dp, dm

    def pullback(self, primitive_bar, lo_coefficients_bar):
        """Return covectors of r*V_Ry, E_Ry and E_LO_Ry, in that order."""
        gp = _real(primitive_bar, self.primitive.shape, "radial primitive covector")
        gm = _real(
            lo_coefficients_bar, self.lo_coefficients.shape, "LO coefficient covector"
        )
        dv, de, delo = np.empty(self.n), ct.c_double(), np.empty(self.nlo)
        with _LOCK:
            self._check_open()
            _status(
                self._lib.edmft_radial_pullback(
                    self._handle, gp, gm, dv, ct.byref(de), delo
                )
            )
        if not np.all(np.isfinite(np.r_[dv, de.value, delo])):
            raise RuntimeError("nonfinite radial pullback")
        return dv, de.value, delo
