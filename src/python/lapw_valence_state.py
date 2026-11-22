"""Prepare the valence part of the screened response from a native snapshot.

Read the central H/O, vector and radial exports prepared by the response
driver. NativeLattice supplies the stored-label reciprocal basis, primitive
volume and site frames. For monoclinic CXZ the matrix and vector files label
k differently, even though they use the same G labels.

For eDMFT each spin branch uses scalar KS vectors, normalized single-site
projector 5 and a fixed band-index window. Retained KS counts agree across
k, while the full LAPW dimension can vary. Complement states must be
thermally empty. The driver attaches the radial-energy search response to
these fixed-energy maps. Sphere radii and local frames are held fixed.

Preparation reduces each rank's complement once. Subsequent actions use
retained/retained and complement/retained blocks. Close the prepared state
or use it as a context manager to release its native radial contexts and any
rank-private disk storage for matching coefficients.
"""

from contextlib import ExitStack
from functools import partial
from response_parallel import SERIAL, parallel_of
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
import copy
import re
import time
from lapw_lattice import NativeLattice, native_site_frames, native_structure
from lapw_cubic_harmonics import cubic_potential_gaunt
import numpy as np
from scipy.io import FortranFile
from scipy.special import spherical_jn, expit
from scipy.fft import next_fast_len
from wstruct import Struct
from indmffile import Indmfl
from frozen_sigma_tangent import CachedFrozenSigmaResponse
from wien_response_density import Density
from lapw_radial_native import NativeRadialContext
from lapw_augmentation_response import PreparedAugmentationChannel
from lapw_spherical_response import PreparedSphericalBand, prepare_wronskian_band
from lapw_nonspherical_response import PreparedNonspherical
from lapw_projector_response import PreparedProjector5
from lapw_compact_response import (
    CompactSphericalPotentialBands,
    CompactSpatialValenceDensity,
)
from lapw_density_response import (
    MuffinTinDensitySum,
    SpatialValenceDensity,
    InterstitialDensity,
    FourierStars,
)
from lapw_basis_density_response import BasisDensityResponse
from lapw_complement_response import (
    LAPWComplementBasis,
    EmptyComplementResponse,
    complement_electronic_operator,
)
from lapw_potential_response import NativePotentialBands, NativeInterstitialPotential
from lapw_native_snapshot import read_basis_matrix, vector_k_from_matrix_header
from lapw_potential_core_state import fingerprint_files


@dataclass
class PreparedValenceState:
    """Prepared maps; contexts is the H/O bank, projector_contexts is ATPARN.

    Both banks share potential coordinates, not radial normalization rules.
    Their native contexts are owned by this state and closed together.
    """

    model: object
    complement: object
    bases: list
    records: list
    structure: object
    contexts: dict
    bands: object
    spatial: object
    direct: object
    reference: object
    _stack: object
    fingerprints: dict = field(default_factory=dict)
    preparation_seconds: float = 0.0
    matrix_quadrature: str = ""
    matrix_overlap_error: float = 0.0
    matching_scratch_bytes: int = 0

    @property
    def projector_contexts(self):
        return self.bands.spherical.projector_contexts

    def close(self):
        self._stack.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def valence_operator(self, *, fixed_number=True):
        # One E for the whole k sum: X and dmu must never be determined per k.
        e, inputs, outputs = complement_electronic_operator(
            self.complement, fixed_number=fixed_number
        )
        if (
            inputs.size != self.bands.output.size
            or outputs.size != self.spatial.coordinates.size
        ):
            raise ValueError("valence coordinate mismatch")
        parallel = parallel_of(self.model)
        return parallel.local_to_global(
            self.spatial.operator
        ) @ e @ parallel.global_to_local(self.bands.operator) + parallel.sum_operator(
            self.direct.operator @ self.bands.spherical_selector
        )


def _retained_channel(channel, count):
    # Prepared geometry sums are small and already in the common band gauge.
    # D needs only retained columns because central rho_QQ=rho_QR=0.
    result = copy.copy(channel)
    result.fj = channel.fj[:, :count].copy()
    result.dfj = channel.dfj[:, :count].copy()
    result.local = channel.local[:, :, :count].copy()
    result.nband = count
    return result


def _radial_fields(path, structure):
    text = path.read_text()
    blocks = re.split(r"ATOMNUMBER\s*=\s*\d+", text)[1:]
    if len(blocks) != structure.nat:
        raise ValueError("potential atom count mismatch")
    out = []
    pattern = r"[-+]?\d+\.\d+[DE][-+]\d+"
    for a, block in enumerate(blocks):
        parts = re.split(r"VLM\(R\) FOR L\s*(-?\d+)\s+M=\s*(\d+)", block)
        channels = {}
        for j in range(1, len(parts), 3):
            key = (int(parts[j]), int(parts[j + 1]))
            values = np.array(
                [float(x.replace("D", "E")) for x in re.findall(pattern, parts[j + 2])]
            )
            if key in channels or len(values) < structure.jrj[a]:
                raise ValueError("invalid radial potential layout")
            channels[key] = values[: structure.jrj[a]]
        out.append(channels)
    return out


def _snapshot_inputs(case, source, complex_vectors, spin_channel="", *, dft=False):
    s = Struct()
    s.ReadStruct(str(source / (case + ".struct")))
    s = native_structure(s)
    NativeLattice.from_structure(s)  # Validate both stored and primitive cell maps.
    if not dft:
        control = (
            (source / (case + ".indmfl" + spin_channel)).read_text().splitlines()[0].split("#")[0].split()
        )
        if len(control) < 4 or int(control[2]) != 1 or abs(int(control[3])) != 5:
            raise NotImplementedError("normalized projector 5 is required")
    v = _radial_fields(source / (case + ".vsp"), s)
    vns = _radial_fields(source / (case + ".vns"), s)
    header = []
    lo = []
    with FortranFile(source / (case + ".vector")) as handle:
        for a in range(s.nat):
            header.append(handle.read_record("f8"))
            lo.append(handle.read_record("f8").reshape(-1, 4).T)
    ilo = [np.sum(e < 995, axis=1) for e in lo]
    indices = {}
    index = 0
    for a in range(s.nat):
        for l in range(4):
            for j in range(ilo[a][l]):
                count = s.mult[a] * (2 * l + 1)
                indices[a, l, j] = np.arange(index, index + count)
                index += count
    return SimpleNamespace(
        s=s,
        complex=complex_vectors,
        v=v,
        vns=vns,
        header=header,
        lo=lo,
        ilo=ilo,
        lo_indices=indices,
        nlo_total=index,
    )


def load_native_valence_state(
    directory,
    case,
    *,
    complex_vectors,
    energy_policy,
    retained=None,
    library=None,
    cache_wavefunctions=False,
    expected_fingerprints=None,
    matrix_quadrature=None,
    fft_grid="fast",
    parallel=SERIAL,
    cache_mb=256.0,
    basis_cache_mb=128.0,
    disk_matching=False,
    spin_channel="",
    diagnostic_gap_cutoff=None,
    electronic_model="dmft",
):
    """Prepare central native data; optionally reject a stale execution snapshot.

    matrix_quadrature explicitly identifies the native H/O assembler:
    "wien2k24" selects default LAPW1 CHARGE; "atparn" selects the historical
    rule. Projector and density reconstruction always follows dmft2 ATPARN.
    The input filenames alone do not identify a native quadrature policy.

    complex_vectors is explicit, because the native VNS phase convention must
    match the executable used to export these files. Omitting fingerprints is
    for an independently inspected diagnostic snapshot. Production callers
    must pass the expected file fingerprints from their execution snapshot.
    """
    if spin_channel not in ("", "dn"):
        raise ValueError("spin_channel must be empty or dn")
    if electronic_model not in ("dmft", "dft"):
        raise ValueError("electronic_model must be dmft or dft")
    if electronic_model == "dft" and spin_channel:
        raise NotImplementedError("initial DFT response requires unpolarized KS potential")
    response_name = "response" + ("_dn" if spin_channel else "")
    if matrix_quadrature not in ("atparn", "wien2k24"):
        raise ValueError("explicit matrix_quadrature=atparn or wien2k24 is required")
    if fft_grid not in ("fast", "minimal"):
        raise ValueError("fft_grid must be fast or minimal")
    if energy_policy != "fixed":
        raise NotImplementedError("only fixed radial energies are supported")
    if not isinstance(complex_vectors, (bool, np.bool_)):
        raise ValueError("explicit complex-vector policy required")
    directory = Path(directory)
    names = [
        case + suffix
        for suffix in [
            ".struct",
            ".vector",
            ".vsp",
            ".vns",
            ".clmval" + spin_channel,
            (".in1c" if complex_vectors else ".in1"),
        ]
    ]
    if electronic_model == "dmft":
        names += [case + ".indmfl" + spin_channel, "projectorw.dat", response_name + ".meta"]
    else:
        names += [case + (".in2c" if complex_vectors else ".in2"), case + ".weight"]
    paths = {name: directory / name for name in names}
    patterns = ["native.k????????.bin", "native.ham.k????????.bin"]
    if electronic_model == "dmft":
        patterns.insert(0, response_name + ".k????????.bin")
    for pattern in patterns:
        found = list(directory.glob(pattern))
        if not found:
            raise ValueError("missing native snapshot records: " + pattern)
        paths.update({p.name: p for p in found})
    if (directory / "EF.dat").is_file():
        paths["EF.dat"] = directory / "EF.dat"
    mixer_policy = None
    mixer_path = directory / (case + ".inm")
    if mixer_path.is_file():
        fields = mixer_path.read_text().splitlines()[0].split()
        if len(fields) < 3 or fields[2].upper() not in ("YES", "NO"):
            raise ValueError("explicit YES/NO charge-normalization policy required in INM")
        mixer_policy = fields[2].upper() == "YES"
        paths[mixer_path.name] = mixer_path
    fingerprints = parallel.fingerprints(paths)
    if expected_fingerprints is not None and fingerprints != expected_fingerprints:
        raise ValueError("native valence snapshot fingerprints differ")
    stack = ExitStack()
    started = time.perf_counter()
    try:
        result = _prepare_snapshot(
            case,
            directory,
            retained,
            complex_vectors,
            library,
            cache_wavefunctions,
            stack,
            matrix_quadrature=matrix_quadrature,
            fft_grid=fft_grid,
            parallel=parallel,
            cache_mb=cache_mb,
            basis_cache_mb=basis_cache_mb, disk_matching=disk_matching,
            spin_channel=spin_channel,
            diagnostic_gap_cutoff=diagnostic_gap_cutoff,
            electronic_model=electronic_model,
        )
    except BaseException:
        stack.close()
        raise
    result.spin_channel = spin_channel
    result.fingerprints = fingerprints
    result.fingerprint_paths = paths
    result.mixer_normalization_enabled = mixer_policy
    result.preparation_seconds = time.perf_counter() - started
    return result


def _prepare_snapshot(
    case,
    source,
    nret,
    complex_vectors,
    library,
    cache_wavefunctions,
    stack,
    *,
    matrix_quadrature="atparn",
    fft_grid="fast",
    parallel=SERIAL,
    cache_mb=256.0,
    basis_cache_mb=128.0,
    disk_matching=False,
    spin_channel="",
    diagnostic_gap_cutoff=None,
    electronic_model="dmft",
):
    central = Path(source)
    dft = electronic_model == "dft"
    mat = _snapshot_inputs(case, central, complex_vectors, spin_channel, dft=dft)
    s = mat.s
    if dft:
        from lapw_dft_vector import load_vectors
        from lapw_dft_fermi import (
            DFTFermiModel, DFTComplementResponse, native_fermi_parameters,
        )

        records, kw = load_vectors(central / (case + ".vector"), s.nat, mat.nlo_total, parallel)
        available = records[0].nstored
        nret = available if nret is None else nret
        if int(nret) != nret or not 0 < nret <= available:
            raise ValueError("invalid DFT retained count")
        nret = int(nret)
        mu, temperature = native_fermi_parameters(
            central, case, complex_vectors=complex_vectors
        )
        model = DFTFermiModel(
            np.array([r.energies[:nret] for r in records]), kw, mu, temperature,
            parallel=parallel,
        )
    else:
        exported, records = CachedFrozenSigmaResponse.from_export(
            central / ("response_dn" if spin_channel else "response"),
            parallel=parallel, cache_mb=cache_mb, diagnostic_gap_cutoff=diagnostic_gap_cutoff,
        )
        if nret is None:
            nret = exported.eps.shape[1]
        if int(nret) != nret or not 0 < nret <= exported.eps.shape[1]:
            raise ValueError("invalid retained count")
        nret = int(nret)
    mat.cv = [
        (r.energies[:nret].copy(), r.coefficients[:, :nret, 0].copy(), r.labels)
        for r in records
    ]
    if any(r.nlo != mat.nlo_total for r in records):
        raise ValueError("native LO count differs from vector header")
    # Only retained eigenvectors enter response preparation. Exported high-state
    # data, if present, is not used to construct the complement.
    if not dft and nret <= max(exported.W):
        raise ValueError("retained space misses DMFT window")
    if not dft:
        normalized = np.array([r.projector[:nret] for r in records])
        raw_native = np.zeros_like(normalized)
        for orb, ic in enumerate(exported.meta.cix):
            if ic >= 0:
                for j in range(exported.meta.nindo[orb]):
                    raw_native[:, :, j, orb] = exported.raw[
                        :, :nret, exported.meta.isx[j, orb], ic
                    ]
        model = CachedFrozenSigmaResponse(
            exported.meta, exported.eps[:, :nret], raw_native, normalized,
            exported.kw, exported.W, parallel=parallel, cache_mb=cache_mb,
            diagnostic_gap_cutoff=diagnostic_gap_cutoff,
        )
        del exported
        overlap = model.overlap(model.raw[:, model.W])
        for ic, n in enumerate(model.meta.dimensions):
            active = np.flatnonzero(model.meta.active[:n, ic])
            ix = np.ix_(active, active)
            x = model.X[:, :, ic][ix]
            ss = overlap[:, :, ic][ix]
            if not np.allclose(x @ ss @ x, np.eye(len(active)), atol=2e-9, rtol=0):
                raise ValueError("normalization is stale for the exported global k sum")
    mat.cv = [
        (ev[:nret].copy(), c[:, :nret].copy(), labels) for ev, c, labels in mat.cv
    ]
    bases = []
    k_coordinates = []
    for k, (ev, c, labels) in enumerate(mat.cv):
        hh = read_basis_matrix(
            central / f"native.ham.k{records[k].kpoint:08d}.bin", b"EDMFTHV1"
        )
        oo = read_basis_matrix(
            central / f"native.k{records[k].kpoint:08d}.bin", b"EDMFTOV1"
        )
        if (
            hh.kpoint != records[k].kpoint
            or oo.kpoint != records[k].kpoint
            or hh.nlo != mat.nlo_total
            or oo.nlo != mat.nlo_total
            or not np.array_equal(hh.reciprocal_k, oo.reciprocal_k)
            or not np.array_equal(hh.labels, labels)
            or not np.array_equal(oo.labels, labels)
        ):
            raise ValueError("H/O and vector coefficient labels or headers differ")
        base = LAPWComplementBasis.prepare(
            hh.matrix,
            oo.matrix,
            c,
            ev,
            tridiagonal_basis=True,
            mu=model.meta.mu,
            temperature=model.meta.temperature,
            active_columns=(np.isin(np.arange(nret), model.W) |
                (2*nret*expit((model.meta.mu-ev)/model.meta.temperature)
                 > min(1e-12, model.meta.temperature*1e-12))),
        )
        bases.append(base)
        vector_k = vector_k_from_matrix_header(oo.reciprocal_k, s)
        if dft and not np.allclose(records[k].reciprocal_k, vector_k, rtol=0, atol=1e-12):
            raise ValueError("vector and matrix k-point headers disagree after native conversion")
        records[k].reciprocal_k = vector_k
        k_coordinates.append(vector_k)
        mat.cv[k] = (ev, base.coefficients, labels)

    return _prepare_representation(case, central, mat, records, k_coordinates, bases,
        model, nret, stack, complex_vectors=complex_vectors, library=library,
        cache_wavefunctions=cache_wavefunctions, matrix_quadrature=matrix_quadrature,
        fft_grid=fft_grid, parallel=parallel, basis_cache_mb=basis_cache_mb,
        disk_matching=disk_matching, spin_channel=spin_channel, dft=dft)


def _prepare_representation(case, central, mat, records, k_coordinates, bases,
                            model, nret, stack, *, complex_vectors, library,
                            cache_wavefunctions, matrix_quadrature, fft_grid,
                            parallel, basis_cache_mb, disk_matching, spin_channel,
                            dft, scalar_spin_frame=False, geometry_only=False):
    """Scalar LAPW H/O and spatial charge in the given retained/complement frame.

    For SOC, the central scalar occupation is the partial spin trace of the
    spinor DMFT density matrix. The spinor projector and resolvent are prepared
    separately; the scalar LAPW terms use the same native kernels as DMFT2.
    """
    s = mat.s
    has_projector = not (dft or scalar_spin_frame)
    if not np.isfinite(basis_cache_mb) or basis_cache_mb < 0:
        raise ValueError("basis cache budget must be nonnegative")
    basis_cache_limit = int(basis_cache_mb * 1024**2)
    basis_cache_bytes = 0
    nk = len(records)
    for pattern in ["native.k????????.bin", "native.ham.k????????.bin"]:
        if not scalar_spin_frame and len(list(central.glob(pattern))) != parallel.sum(nk):
            raise ValueError("extra/missing native matrix k records")
    rho = model.central_density()
    reference = (None if geometry_only else
                 Density(central / (case + ".clmval" + spin_channel), central / (case + ".struct")))
    if has_projector:
        inl = Indmfl(str(central / case), extn="indmfl" + spin_channel)
        inl.read()
    ctx = {}
    matrix_ctx = {}
    from lapw_response_scratch import MatchingScratch
    scratch = MatchingScratch(stack) if disk_matching else None
    sphere = []
    nonsphere = []
    project = []
    ltop, lns = map(
        int,
        (central / (case + (".in1c" if mat.complex else ".in1")))
        .read_text()
        .splitlines()[1]
        .split()[1:3],
    )
    if ltop < 8 or not 0 <= lns <= ltop:
        raise NotImplementedError(
            "native density adapter needs LMAX>=8 and LNSMAX<=LMAX"
        )
    # LAPW1 can retain angular channels above l=8 (FeSe uses l=10), but the
    # native dmft2 density still uses the compile-time LMAX2=8 in modules.f90.
    # Reconstructing its output from all LAPW1 channels adds high-l products
    # which dmft2 never wrote to clmval. The discrepancy is most visible in
    # the outer radial samples of high-l nonspherical density harmonics and
    # invalidates the strict native charge-map closure needed by the adjoint.
    # DFT's native LAPW2 path has its own angular cutoff and retains ltop.
    density_limit = ltop if dft else 8
    first = np.r_[0, np.cumsum(s.mult)[:-1]]
    site_atom = np.repeat(np.arange(s.nat), s.mult)
    axes, locations = native_site_frames(s)
    geometry = []
    for k, (ev, c, labels) in enumerate(mat.cv):
        kf = k_coordinates[k]
        frac = labels.T + kf
        q = NativeLattice.from_structure(s).momenta(frac)
        qn = np.linalg.norm(q, axis=1)
        geometry.append((c, labels, frac, q, qn, len(c) - mat.nlo_total))
    density_atoms = {}
    density_matching = {}
    density_complement = {}
    density_ctx = matrix_ctx if dft else ctx
    volume = NativeLattice.from_structure(s).volume
    for a in range(s.nat):
        n = s.jrj[a]
        ctx[a] = {}
        for l in range(ltop + 1):
            e = mat.header[a][l]
            lapw = bool(e < 150)
            e = e if lapw else e - 200
            elo = mat.lo[a][l, : mat.ilo[a][l]] if l < 4 else np.array([])
            ctx[a][l] = stack.enter_context(
                NativeRadialContext(
                    mat.v[a][0, 0],
                    r0=s.r0[a],
                    dx=np.log(s.rmt[a] / s.r0[a]) / (n - 1),
                    rmt=s.rmt[a],
                    energy=e,
                    lo_energies=elo,
                    ell=l,
                    charge=s.Znuc[a],
                    relativistic=bool("RELA" in s.mode),
                    lapw=lapw,
                    library=library,
                )
            )
        matrix_ctx[a] = ctx[a]
        if matrix_quadrature != "atparn":
            matrix_ctx[a] = {}
            for l, c in ctx[a].items():
                matrix_ctx[a][l] = stack.enter_context(
                    NativeRadialContext(
                        mat.v[a][0, 0],
                        r0=s.r0[a],
                        dx=np.log(s.rmt[a] / s.r0[a]) / (n - 1),
                        rmt=s.rmt[a],
                        energy=c.energy,
                        lo_energies=c.lo_energies,
                        ell=l,
                        charge=s.Znuc[a],
                        relativistic=c.relativistic,
                        lapw=c.lapw,
                        library=library,
                        quadrature=matrix_quadrature,
                    )
                )
        dens = []
        density_matching[a] = []
        ns = []
        ks = []
        for site in range(first[a], first[a] + s.mult[a]):
            for k, (c, labels, frac, q, qn, npw) in enumerate(geometry):
                nb = c.shape[1]
                phase = np.exp(2j * np.pi * (labels.T @ s.pos[site]))
                loc = (s.rotloc[a] @ axes[site] @ q.T).T
                orbit_phases = None
                if site == first[a] and s.mult[a] > 1:
                    positions = s.pos[first[a]:first[a] + s.mult[a]]
                    orbit_phases = np.exp(2j * np.pi * (positions @ labels[:, :npw]))
                channels = {}
                z = {}
                zks = {}
                for l, rc in density_ctx[a].items():
                    groups = [npw + mat.lo_indices[a, l, j] for j in range(rc.nlo)]
                    matching = np.array(
                        [
                            spherical_jn(l, qn[:npw] * rc.rmt),
                            qn[:npw]
                            * spherical_jn(l, qn[:npw] * rc.rmt, derivative=True),
                        ]
                    )
                    ch = PreparedAugmentationChannel(
                        l,
                        rc.rmt,
                        volume,
                        q,
                        qn,
                        phase,
                        c,
                        npw,
                        groups,
                        matching_values=matching,
                    )
                    rebuild_triangle = partial(
                        prepare_wronskian_band, l, rc.rmt, volume, q[:npw],
                        phase[:npw], c[:npw], matching, retained=nret,
                    )
                    site_triangle = None
                    if not rc.lapw:
                        triangle = np.broadcast_to(np.array(0j), (nb, nret))
                    elif s.mult[a] == 1:
                        triangle = rebuild_triangle()
                    else:
                        site_triangle = rebuild_triangle
                        triangle = None
                        if site == first[a]:
                            triangle = prepare_wronskian_band(
                                l, rc.rmt, volume, q[:npw], orbit_phases,
                                c[:npw], matching, retained=nret,
                            )
                    spherical_band = PreparedSphericalBand(
                        ch, matrix_ctx[a][l].primitive, matrix_ctx[a][l].lo_coefficients,
                        matrix_ctx[a][l].overlap_metric, matrix_ctx[a][l].radius,
                        matrix_ctx[a][l].energy, matrix_ctx[a][l].lo_energies, triangle,
                        lapw=rc.lapw, cache_products=False, share_radial=True,
                        site_triangle=site_triangle, retained=nret,
                    )
                    if scratch is not None:
                        spherical_band.store_matching(scratch)
                    sphere.append((a, l, k, spherical_band))
                    if l <= (lns if geometry_only else density_limit):
                        ch = PreparedAugmentationChannel(
                            l,
                            rc.rmt,
                            volume,
                            loc,
                            qn,
                            phase,
                            c,
                            npw,
                            groups,
                            matching_values=matching,
                        )
                        channels[l] = ch
                        if not geometry_only:
                            z[l] = ch.value(
                                rc.primitive[-1], rc.lo_coefficients, lapw=rc.lapw
                            ).reshape(-1, nb)
                        mc = matrix_ctx[a][l]
                        zks[l] = ch.value(
                            mc.primitive[-1], mc.lo_coefficients, lapw=mc.lapw
                        ).reshape(-1, nb)
                if not geometry_only:
                    dz = np.concatenate([z[l] for l in range(density_limit + 1)])
                    # In scratch mode a retained-column slice must not pin the
                    # full complement buffer after its explicit spill.
                    retained_z = dz[:, :nret]
                    complement_z = dz[:, nret:]
                    if scratch is not None:
                        retained_z = np.ascontiguousarray(retained_z)
                        complement_z = scratch.store(complement_z)
                    else:
                        complement_z.setflags(write=False)
                    retained_z.setflags(write=False)
                    dens.append((k, retained_z, 1 / s.mult[a]))
                    density_matching[a].append(
                        {l: _retained_channel(ch, nret) for l, ch in channels.items()}
                    )
                    density_complement.setdefault(a, []).append(complement_z)
                if scratch is not None:
                    for l in range(lns+1):
                        channels[l].store_matching(scratch)
                zks_local = np.concatenate([zks[l] for l in range(lns + 1)])
                zks_local.setflags(write=False)
                ns.append(
                    (
                        zks_local,
                        {l: channels[l] for l in range(lns + 1)},
                        {l: matrix_ctx[a][l].lo_coefficients for l in range(lns + 1)},
                        {l: ctx[a][l].lapw for l in range(lns + 1)},
                    )
                )
                ks.append(k)
        if not geometry_only:
            hs = [key[1:] for key in reference.rad if key[0] == a]
            small = 1 / 137.0359895**2 if "RELA" in s.mode else 4e-22
            density_atoms[a] = MuffinTinDensitySum(
                dens, rho, {l: density_ctx[a][l].primitive[:, :2] for l in range(density_limit + 1)}, hs, small, share_matching=True
            )
        nonsphere.append(
            (
                a,
                PreparedNonspherical(
                    {l: matrix_ctx[a][l].primitive for l in range(lns + 1)},
                    matrix_ctx[a][0].overlap_metric,
                    list(mat.vns[a]),
                    np.array(list(mat.vns[a].values())).T,
                    ns,
                    share_matching=True,
                    gaunt_matrices=(cubic_potential_gaunt(
                        {l: matrix_ctx[a][l].primitive.shape[2] for l in range(lns+1)},
                        list(mat.vns[a])) if s.iatnr[a] > 0 else None),
                ),
                ks,
            )
        )
    for a, l, k, op in sphere:
        for name, compute in [
            ("_z", lambda op=op: op.z),
            ("_zpw", lambda op=op: op.zpw),
            ("_kz", lambda op=op: op.kz),
            ("_sz", lambda op=op: op.sz),
        ]:
            if (
                getattr(op, name) is None
                and basis_cache_bytes + ((2*l+1)*(2+op.channel.nlo)*op.channel.nband*16) <= basis_cache_limit
            ):
                array = compute()
                array.setflags(write=False)
                setattr(op, name, array)
                basis_cache_bytes += array.nbytes
    model.basis_cache_bytes = basis_cache_bytes
    model.basis_cache_limit = basis_cache_limit
    radial_projector = {}
    if has_projector:
        with (central / "projectorw.dat").open() as f:
            for _ in range(int(f.readline().split()[1])):
                nr, nr0, a, l = map(int, f.readline().split()[1:])
                rows = np.array([list(map(float, f.readline().split())) for _ in range(nr)])
                if nr != nr0:
                    raise NotImplementedError(
                        "nonzero outside-MT projector tails require an explicit adapter"
                    )
                radial_projector[a - 1, l] = rows[:, 1:3]
    raw = [
        np.zeros((b.retained.shape[0], *model.raw.shape[2:]), complex) for b in bases
    ]
    for ic, orbs in (inl.cix if has_projector else {}).items():
        if len(orbs) != 1:
            raise NotImplementedError("only single-site impurity blocks are supported")
        site, l, _ = orbs[0]
        site -= 1
        a = int(site_atom[site])
        rc = ctx[a][l]
        locrot, xyz, shift, rmt2 = inl.atoms[site + 1]
        if shift is not None:
            raise NotImplementedError("shifted local orbitals are not supported")
        xyz = np.array(xyz) if len(xyz) else s.rotloc[a]
        if xyz.shape != (3, 3):
            raise NotImplementedError("local frame must have three explicit axes")
        for k, (c, labels, frac, q, qn, npw) in enumerate(geometry):
            phase = np.exp(2j * np.pi * (frac @ locations[site]))
            qloc = (xyz @ axes[site] @ q.T).T
            groups = [npw + mat.lo_indices[a, l, j] for j in range(rc.nlo)]
            ch = PreparedAugmentationChannel(
                l, rc.rmt, volume, qloc, qn, phase, c, npw, groups
            )
            op = PreparedProjector5(
                ch,
                rc.primitive,
                rc.lo_coefficients,
                rc.projector_metric,
                radial_projector[a, l],
                inl.cftrans[ic].T,
                lapw=rc.lapw,
            )
            cols = slice(0, 2 * l + 1)
            raw[k][:, cols, ic - 1] = op.value()
            project.append((a, l, k, int(ic) - 1, cols, op))
    if has_projector and not np.allclose(
        np.array([r[:nret] for r in raw]), model.raw, atol=1e-12, rtol=1e-10
    ):
        difference=np.linalg.norm(np.array([r[:nret] for r in raw])-model.raw)/np.linalg.norm(model.raw)
        raise ValueError(f"assembled raw projectors differ from native snapshot: relative error {difference:g}")
    complement = None
    if not scalar_spin_frame:
        if dft:
            from lapw_dft_fermi import DFTComplementResponse
        cls = DFTComplementResponse if dft else EmptyComplementResponse
        complement = cls(model, [b.solver for b in bases], [r[nret:] for r in raw])
    bands = CompactSphericalPotentialBands(
        matrix_ctx,
        sphere,
        nonsphere,
        project,
        [b.solver.nq for b in bases],
        model.raw.shape,
        energy_policy="fixed",
        projector_contexts=ctx,
    )
    allg = [r.labels[:, : r.nbasis - r.nlo].T for r in records]
    direct = compact_spatial = None
    if not geometry_only:
        reps = np.array(list(reference.pw))
        extent = np.max([np.ptp(g, axis=0) for g in allg], axis=0)
        extent = parallel.reduce(extent, "max")
        star_extent = np.max(
            abs(np.einsum("sij,qj->sqi", s.timat.transpose(0, 2, 1), reps)), axis=(0, 1)
        )
        minimum = 2 * np.maximum(extent, star_extent) + 1
        # These FFTs evaluate finite plane-wave products exactly. Their grid is
        # NOT LAPW0's XC integration grid. Any larger alias-free grid gives the
        # same Fourier coefficients. The minimal odd lengths are often prime;
        # round up to efficient complex FFT lengths with only modest extra storage.
        shape = (
            tuple(next_fast_len(int(n)) for n in minimum)
            if fft_grid == "fast"
            else tuple(minimum)
        )
        stars = FourierStars(
            shape, reps, s.timat.transpose(0, 2, 1), s.tau, complex_phases=mat.complex
        )
        ir = [
            InterstitialDensity(
                c[: len(g), :nret],
                g,
                volume,
                stars,
                cache_wavefunctions=cache_wavefunctions, share_coefficients=True,
            )
            for (c, *_), g in zip(geometry, allg)
        ]
        spatial = SpatialValenceDensity(reference, density_atoms, ir)
        direct = BasisDensityResponse(
            spatial,
            {a: {l: density_ctx[a][l] for l in range(density_limit + 1)} for a in ctx},
            density_matching,
            energy_policy="fixed",
        )
        compact_spatial = CompactSpatialValenceDensity(
            spatial,
            [b.solver.nq for b in bases],
            density_complement,
            [c[: len(g), nret:] for (c, *_), g in zip(geometry, allg)],
            share_coefficients=True,
        )
    lines = (central / (case + ".vns")).read_text().splitlines()
    j = next(i for i, line in enumerate(lines) if "NUMBER OF PW" in line)
    nstar = int(lines[j].split()[0])
    reps = np.array(
        [
            [int(line[a : a + 5]) for a in [3, 8, 13]]
            for line in lines[j + 1 : j + 1 + nstar]
        ]
    )
    ir_potential = [
        NativeInterstitialPotential(
            c[: len(g)],
            g,
            reps,
            s.timat.transpose(0, 2, 1),
            s.tau,
            complex_values=complex_vectors, share_coefficients=True,
        )
        for (c, *_), g in zip(geometry, allg)
    ]
    native_bands = NativePotentialBands(bands, ir_potential, retained=nret)
    # A representative-k overlap closure detects an incorrectly declared
    # H/O quadrature before it can silently contaminate a response solve.
    # Use the exact geometric interstitial mask, independently of VNS input.
    c, labels, frac, q, qn, npw = geometry[0]
    difference = np.linalg.norm(q[:npw, None] - q[None, :npw], axis=-1)
    theta = np.eye(npw, dtype=complex)
    for site, a in enumerate(site_atom):
        x = difference * s.rmt[a]
        mask = (4 * np.pi * s.rmt[a] ** 3 / volume) * np.divide(
            spherical_jn(1, x), x, out=np.full_like(x, 1 / 3), where=x != 0
        )
        phase = np.exp(2j * np.pi * (labels[:, :npw].T @ s.pos[site]))
        theta -= mask * phase.conj()[:, None] * phase[None, :]
    assembled = c[:npw].conj().T @ (theta @ c[:npw, :nret])
    for a, l, k, op in sphere:
        if k == 0:
            assembled += np.einsum(
                "min,mip->np", op.z.conj(), op.sz[..., :nret], optimize=True
            )
    assembled[:nret] -= np.eye(nret)
    overlap_error = float(np.linalg.norm(assembled) / np.sqrt(nret))
    if overlap_error > 5e-9:
        raise ValueError(
            f"assembled overlap disagrees with native H/O ({overlap_error:g}); check matrix_quadrature and snapshot consistency"
        )
    return PreparedValenceState(
        model,
        complement,
        bases,
        records,
        s,
        matrix_ctx,
        native_bands,
        compact_spatial,
        direct,
        reference,
        stack,
        matrix_quadrature=matrix_quadrature,
        matrix_overlap_error=overlap_error,
        matching_scratch_bytes=scratch.file.tell() if scratch is not None else 0,
    )
