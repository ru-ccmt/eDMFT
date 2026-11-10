"""Versioned central data for preparing native potential/core response in memory.

The archive contains numeric arrays and JSON only, never pickled Python objects
or commands. It supplies native angular quadrature and central core matching
metadata explicitly. Loading runs no WIEN executable, reads no other state
implicitly, and requires no material-specific paths. It is the potential/core
part of preparation; KS matrices/vectors and k-point valence preparation remain
separate. This is not a general production exporter or a screened-stress driver.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
import hashlib, json, os, tempfile
from lapw_lattice import NativeLattice
import numpy as np
from wien_response_density import Density
from lapw_density_response import checked
from lapw_hartree_periodic import PeriodicHartreeResponse
from lapw_xc_response import (
    MuffinTinXCResponse,
    InterstitialXCResponse,
    PeriodicXCResponse,
)
from lapw_core_response import (
    CorePotentialExtension,
    MatchedCoreState,
    AtomCoreResponse,
    PeriodicCoreResponse,
)
from lapw_charge_response import (
    NativeChargeNormalization,
    native_mixer_inversion,
    native_charge_covectors,
)

FORMAT = "EDMFT_POTENTIAL_CORE_STATE"
VERSION = 1
STRUCTURE_ARRAYS = (
    "mult",
    "jrj",
    "r0",
    "rmt",
    "pos",
    "timat",
    "tau",
    "rotloc",
    "iatnr",
    "Znuc",
)
INTEGER_ARRAYS = {"mult", "jrj", "timat", "iatnr", "density_keys", "pw_representatives"}
POLICY = dict(
    spin_channels=1,
    spin_orbit=False,
    radial_energy="fixed",
    rmt="fixed",
    local_axes="fixed",
    charge_rule="native_default",
    xc_parameters="native_default",
    core_nucleus="point",
)


def fingerprint_files(files):
    """Fingerprint caller-named physical inputs; store logical names, not paths."""
    return {
        str(key): hashlib.sha256(Path(path).read_bytes()).hexdigest()
        for key, path in sorted(files.items())
    }


@dataclass
class PreparedPotentialCore:
    input_density: object
    structure: object
    hartree: object
    xc: object
    core: object
    metadata: dict

    @property
    def potential_operator(self):
        return self.xc.hxc_operator

    def normalization(self, raw_output, target_charge, *, enabled=True):
        if tuple(raw_output.rad) != self.hartree.keys or tuple(raw_output.pw) != tuple(
            self.input_density.pw
        ):
            raise ValueError("raw output density coordinates differ from central state")
        return NativeChargeNormalization(
            raw_output,
            self.structure,
            target_charge,
            inversion=self.metadata["inversion"],
            enabled=enabled,
        )


def _prepare(metadata, arrays):
    """Validate all physical conventions before constructing any operator."""
    m = metadata
    if not isinstance(m, dict):
        raise ValueError("state metadata must be a JSON object")
    if m.get("format") != FORMAT or m.get("version") != VERSION:
        raise ValueError("unsupported potential/core state version")
    if m.get("policy") != POLICY:
        raise NotImplementedError("unsupported central response policy")
    if m.get("functional") not in ("LDA", "PBE"):
        raise NotImplementedError("only native unpolarized LDA/PBE supported")
    if not isinstance(m.get("inversion"), bool):
        raise ValueError("explicit potential reality flag required")
    if m.get("mode") not in ("NREL", "RELA"):
        raise NotImplementedError("unsupported radial mode")
    shape = m.get("fft_shape")
    if (
        not isinstance(shape, list)
        or len(shape) != 3
        or any(type(n) is not int or n < 1 for n in shape)
    ):
        raise ValueError("invalid FFT shape")
    if type(m.get("ncon")) is not int or m["ncon"] < 1:
        raise ValueError("invalid pseudocharge order")
    hashes = m.get("source_fingerprints")
    if not isinstance(hashes, dict) or any(
        not isinstance(k, str)
        or not isinstance(v, str)
        or len(v) != 64
        or any(c not in "0123456789abcdef" for c in v)
        for k, v in hashes.items()
    ):
        raise ValueError("invalid source fingerprints")
    nat = m.get("nat")
    if not isinstance(nat, int) or nat < 1 or len(m.get("core", [])) != nat:
        raise ValueError("invalid atom/core metadata")
    required = set(STRUCTURE_ARRAYS) | {
        "cell",
        "density_keys",
        "pw_representatives",
        "density",
    }
    for a in range(nat):
        required.update(f"{name}_{a}" for name in ["angular", "weights", "vsp"])
        if m["functional"] == "PBE":
            required.update(
                f"{name}_{a}" for name in ["sine", "cosine", "angular_derivatives"]
            )
    if set(arrays) != required:
        raise ValueError("missing or unrecognized state arrays")
    for key, value in arrays.items():
        if not np.issubdtype(value.dtype, np.number) or not np.all(np.isfinite(value)):
            raise ValueError("invalid numeric array: " + key)
        if key in INTEGER_ARRAYS and not np.issubdtype(value.dtype, np.integer):
            raise ValueError("integer array required: " + key)
        if np.iscomplexobj(value):
            raise ValueError("state arrays use explicit real coordinates: " + key)
    cell = checked(arrays["cell"], (6,), "cell", real=True)
    s = SimpleNamespace(
        nat=nat,
        lattyp=m["lattyp"],
        a=cell[0],
        b=cell[1],
        c=cell[2],
        alpha=cell[3],
        beta=cell[4],
        gamma=cell[5],
        mode=m["mode"],
    )
    for name in STRUCTURE_ARRAYS:
        setattr(s, name, np.array(arrays[name], copy=True))
    for name in ["mult", "jrj", "r0", "rmt", "iatnr", "Znuc"]:
        if getattr(s, name).shape != (nat,):
            raise ValueError("invalid structure array: " + name)
    if (
        np.any(s.mult < 1)
        or np.any(s.jrj < 8)
        or np.any(s.r0 <= 0)
        or np.any(s.rmt <= s.r0)
        or np.any(s.Znuc <= 0)
    ):
        raise ValueError("invalid atom radial data")
    if s.pos.shape != (sum(s.mult), 3) or s.rotloc.shape != (nat, 3, 3):
        raise ValueError("invalid positions or local axes")
    if native_mixer_inversion(s) != m["inversion"]:
        raise ValueError("reality flag disagrees with native symmetry convention")
    keys = arrays["density_keys"]
    reps = arrays["pw_representatives"]
    if keys.ndim != 2 or keys.shape[1] != 3 or len(set(map(tuple, keys))) != len(keys):
        raise ValueError("invalid density harmonic labels")
    if any(not 0 <= a < nat or not 0 <= mm <= abs(l) for a, l, mm in keys):
        raise ValueError("invalid atom/harmonic index")
    if reps.ndim != 2 or reps.shape[1] != 3 or len(set(map(tuple, reps))) != len(reps):
        raise ValueError("invalid Fourier representatives")
    d = Density.__new__(Density)
    d.valence_format = False
    d.lines = []
    d.spans = {}
    d.pw_lines = {}
    d.mult = list(s.mult)
    d.grids = [
        np.linspace(np.log(s.r0[a]), np.log(s.rmt[a]), s.jrj[a]) for a in range(nat)
    ]
    d.rad = {tuple(key): np.zeros(s.jrj[key[0]]) for key in keys}
    d.pw = {tuple(g): 0j for g in reps}
    vector = checked(
        arrays["density"], (len(d.vector()),), "central density", real=True
    )
    d = d.from_vector(vector)
    h = PeriodicHartreeResponse(
        d, s, m["fft_shape"], ncon=m["ncon"], inversion=m["inversion"]
    )
    mt = []
    atoms = {}
    empty = []
    for a in range(nat):
        lm = [(int(l), int(mm)) for aa, l, mm in keys if aa == a]
        r = np.exp(d.grids[a])
        rho = np.array([d.rad[a, l, mm] for l, mm in lm]).T
        extra = {}
        if m["functional"] == "PBE":
            extra = dict(
                sine=arrays[f"sine_{a}"],
                cosine=arrays[f"cosine_{a}"],
                angular_derivatives=arrays[f"angular_derivatives_{a}"],
            )
        mt.append(
            MuffinTinXCResponse(
                r,
                [l for l, mm in lm],
                rho,
                arrays[f"angular_{a}"],
                arrays[f"weights_{a}"],
                functional=m["functional"],
                harmonics=lm,
                cubic=s.iatnr[a] > 0,
                **extra,
            )
        )
        info = m["core"][a]
        if set(info) != {"shift_Ry", "states"} or not isinstance(info["states"], list):
            raise ValueError("invalid atom core metadata")
        ext = CorePotentialExtension(r, arrays[f"vsp_{a}"], shift_ry=info["shift_Ry"])
        states = []
        for state in info["states"]:
            if set(state) != {"energy_Ha", "kappa", "occupation", "match", "end"}:
                raise ValueError("invalid core-state metadata")
            states.append(
                MatchedCoreState(
                    ext.r,
                    ext.central,
                    state["energy_Ha"],
                    z=s.Znuc[a],
                    kappa=state["kappa"],
                    match=state["match"],
                    end=state["end"],
                    occupation=state["occupation"],
                    relativistic=m["mode"] == "RELA",
                    output_points=len(r),
                )
            )
        if states:
            atoms[a] = AtomCoreResponse(ext, states)
        else:
            empty.append(a)
    ir = InterstitialXCResponse(
        h.reanalysis,
        np.array(list(d.pw.values())),
        NativeLattice.from_structure(s).cell,
        lattice_type=s.lattyp,
        functional=m["functional"],
        rho_min=m["rho_min"],
    )
    xc = PeriodicXCResponse(h, mt, ir)
    core = PeriodicCoreResponse(h, atoms, empty_atoms=empty)
    return PreparedPotentialCore(d, s, h, xc, core, m)


def load_potential_core_state(path, *, expected_fingerprints=None):
    """Reconstruct cached maps without native programs; reject stale inputs.

    expected_fingerprints is a logical-name -> SHA256 map supplied by the
    execution snapshot. If supplied, it must match exactly. Omitting it means
    explicitly loading the self-contained recorded state, as in a diagnostic.
    """
    with np.load(path, allow_pickle=False) as archive:
        if "metadata" not in archive:
            raise ValueError("missing state metadata")
        raw = archive["metadata"]
        if raw.shape != () or raw.dtype.kind != "U":
            raise ValueError("metadata must be scalar JSON text")
        metadata = json.loads(str(raw))
        arrays = {k: archive[k].copy() for k in archive.files if k != "metadata"}
    if not isinstance(metadata, dict):
        raise ValueError("state metadata must be a JSON object")
    if (
        expected_fingerprints is not None
        and metadata.get("source_fingerprints") != expected_fingerprints
    ):
        raise ValueError("central-state input fingerprints differ")
    return _prepare(metadata, arrays)


def save_potential_core_state(
    path,
    *,
    density,
    structure,
    fft_shape,
    inversion,
    functional,
    angular_data,
    core_data,
    source_files=None,
    ncon=9,
    rho_min=(2e-8, 1e-9),
):
    """Validate, prepare, and atomically save caller-provided native central data.

    angular_data[a] contains angular, weights and (for PBE) sine, cosine,
    angular_derivatives, in the exact density harmonic order. core_data[a]
    contains vsp (r*V Ry), shift_Ry and states with energy_Ha, kappa,
    occupation and zero-based match/end indices. An empty state list is an
    explicit no-core atom. Sources must already represent one common central
    calculation; a generic native exporter remains a separate integration task.
    """
    s = structure
    if len(angular_data) != s.nat or len(core_data) != s.nat:
        raise ValueError("missing atom preparation data")
    shape = np.asarray(fft_shape)
    if shape.shape != (3,) or not np.issubdtype(shape.dtype, np.integer):
        raise ValueError("integer FFT shape required")
    if not isinstance(inversion, (bool, np.bool_)):
        raise ValueError("explicit potential reality flag required")
    inversion = bool(inversion)
    # The archive reconstructs logarithmic radial grids from STRUCT. Refuse a
    # density sampled on other grids instead of reinterpreting its coefficients.
    native_charge_covectors(density, s, inversion=inversion)
    metadata = dict(
        format=FORMAT,
        version=VERSION,
        policy=POLICY.copy(),
        nat=int(s.nat),
        lattyp=s.lattyp.strip(),
        mode=s.mode.strip(),
        functional=functional,
        inversion=inversion,
        fft_shape=list(map(int, fft_shape)),
        ncon=ncon,
        rho_min=list(rho_min),
        source_fingerprints=fingerprint_files(source_files or {}),
        core=[],
    )
    arrays = {name: np.array(getattr(s, name)) for name in STRUCTURE_ARRAYS}
    arrays.update(
        cell=np.array([s.a, s.b, s.c, s.alpha, s.beta, s.gamma]),
        density_keys=np.array(list(density.rad), int),
        pw_representatives=np.array(list(density.pw), int),
        density=density.vector(),
    )
    for a, (angular, core) in enumerate(zip(angular_data, core_data)):
        for key, value in angular.items():
            arrays[f"{key}_{a}"] = np.array(value)
        arrays[f"vsp_{a}"] = np.array(core["vsp"], float)
        metadata["core"].append(
            {key: value for key, value in core.items() if key != "vsp"}
        )
    # JSON also rejects accidental ndarray/nonserializable metadata before a
    # file is published. No NaN/Inf scalar is allowed in the state header.
    encoded = json.dumps(metadata, sort_keys=True, allow_nan=False)
    prepared = _prepare(json.loads(encoded), arrays)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as out:
            np.savez_compressed(out, metadata=np.array(encoded), **arrays)
        os.replace(temp, path)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)
    return prepared
