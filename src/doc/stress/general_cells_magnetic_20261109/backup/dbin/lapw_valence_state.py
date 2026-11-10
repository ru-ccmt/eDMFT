"""Prepare the compact full-k valence map from a consistent native snapshot.

Required files: case.struct, case.in1[c], case.indmfl, case.vector, case.vsp,
case.vns, case.clmval, projectorw.dat, response.meta, response.k########.bin,
and native[.ham].k########.bin. The latter come from the existing default-off
native exporters. This loader launches no program and modifies no input file.

Supported: nonmagnetic, no SOC, primitive orthogonal cells, projector 5,
fixed-index DMFT window, common retained KS band count and thermally empty
remaining finite LAPW space. Total LAPW sizes MAY differ at each k point.
Fixed radial energies, RMT and local frames are explicit policies. Ordinary
LOs with l<=3 follow the exported vector header. Projector tails outside RMT,
shifted/multisite impurity orbitals and higher-l LOs require separate adapters.

Preparation performs dense complement reductions once per k point. It is not
an iterative electronic solve. Prepared actions carry only RR/QR blocks;
radial contexts and global electronic/normalization data are shared across k.
This is a native snapshot loader, not yet the production execution/export
workflow or a corrected stress tensor. Call close(), or use a with statement,
to release native radial contexts, including on preparation failures.
"""

from contextlib import ExitStack
from response_parallel import SERIAL, parallel_of
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
import copy
import re
import time
import numpy as np
from scipy.io import FortranFile
from scipy.special import spherical_jn
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
from lapw_native_snapshot import read_basis_matrix
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


def _snapshot_inputs(case, source, complex_vectors):
    s = Struct()
    s.ReadStruct(str(source / (case + ".struct")))
    if s.lattyp.strip() != "P" or not np.allclose(
        [s.alpha, s.beta, s.gamma], 90, atol=1e-9, rtol=0
    ):
        raise NotImplementedError("only primitive orthogonal cells are supported")
    control = (
        (source / (case + ".indmfl")).read_text().splitlines()[0].split("#")[0].split()
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
            ".indmfl",
            ".vector",
            ".vsp",
            ".vns",
            ".clmval",
            (".in1c" if complex_vectors else ".in1"),
        ]
    ] + ["projectorw.dat", "response.meta"]
    paths = {name: directory / name for name in names}
    for pattern in [
        "response.k????????.bin",
        "native.k????????.bin",
        "native.ham.k????????.bin",
    ]:
        found = list(directory.glob(pattern))
        if not found:
            raise ValueError("missing native snapshot records: " + pattern)
        paths.update({p.name: p for p in found})
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
            basis_cache_mb=basis_cache_mb,
        )
    except BaseException:
        stack.close()
        raise
    result.fingerprints = fingerprints
    result.fingerprint_paths = paths
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
):
    central = Path(source)
    mat = _snapshot_inputs(case, central, complex_vectors)
    s = mat.s
    exported, records = CachedFrozenSigmaResponse.from_export(
        central / "response", parallel=parallel, cache_mb=cache_mb
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
    if nret <= max(exported.W):
        raise ValueError("retained space misses DMFT window")
    normalized = np.array([r.projector[:nret] for r in records])
    raw_native = np.zeros_like(normalized)
    for orb, ic in enumerate(exported.meta.cix):
        if ic >= 0:
            for j in range(exported.meta.nindo[orb]):
                raw_native[:, :, j, orb] = exported.raw[
                    :, :nret, exported.meta.isx[j, orb], ic
                ]
    model = CachedFrozenSigmaResponse(
        exported.meta,
        exported.eps[:, :nret],
        raw_native,
        normalized,
        exported.kw,
        exported.W,
        parallel=parallel,
        cache_mb=cache_mb,
    )
    del exported
    # The derivative of X assumes X is the inverse square root for THIS k sum.
    # Saved normalization from another mesh/iteration is not a valid snapshot.
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
            mu=model.meta.mu,
            temperature=model.meta.temperature,
        )
        bases.append(base)
        records[k].reciprocal_k = oo.reciprocal_k.copy()
        k_coordinates.append(oo.reciprocal_k)
        mixed = np.column_stack((base.retained, base.complement))
        mat.cv[k] = (ev, mixed, labels)

    if not np.isfinite(basis_cache_mb) or basis_cache_mb < 0:
        raise ValueError("basis cache budget must be nonnegative")
    basis_cache_limit = int(basis_cache_mb * 1024**2)
    basis_cache_bytes = 0
    nk = len(records)
    for pattern in ["native.k????????.bin", "native.ham.k????????.bin"]:
        if len(list(central.glob(pattern))) != parallel.sum(nk):
            raise ValueError("extra/missing native matrix k records")
    rho = model.central_density()
    reference = Density(central / (case + ".clmval"), central / (case + ".struct"))
    inl = Indmfl(str(central / case))
    inl.read()
    ctx = {}
    matrix_ctx = {}
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
    if ltop < 8 or not 0 <= lns <= 8:
        raise NotImplementedError(
            "native density l<=8 adapter needs LMAX>=8 and LNSMAX<=8"
        )
    axes = []
    locations = []
    first = np.r_[0, np.cumsum(s.mult)[:-1]]
    site_atom = np.repeat(np.arange(s.nat), s.mult)
    for site, a in enumerate(site_atom):
        for rotation, tau in zip(s.timat, s.tau):
            diff = rotation @ s.pos[first[a]] + tau - s.pos[site]
            diff -= np.rint(diff)
            if max(abs(diff)) < 2e-7:
                break
        else:
            raise RuntimeError("missing site rotation")
        axes.append(rotation.T)
        locations.append(rotation @ s.pos[first[a]] + tau)
    geometry = []
    for k, (ev, c, labels) in enumerate(mat.cv):
        kf = k_coordinates[k]
        frac = labels.T + kf
        q = 2 * np.pi * frac / np.array([s.a, s.b, s.c])
        qn = np.linalg.norm(q, axis=1)
        geometry.append((c, labels, frac, q, qn, len(c) - mat.nlo_total))
    density_atoms = {}
    density_matching = {}
    density_complement = {}
    volume = s.a * s.b * s.c
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
                channels = {}
                z = {}
                zks = {}
                for l, rc in ctx[a].items():
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
                    triangle = (
                        prepare_wronskian_band(
                            l,
                            rc.rmt,
                            volume,
                            q[:npw],
                            phase[:npw],
                            c[:npw],
                            matching,
                            retained=nret,
                        )
                        if rc.lapw
                        else np.broadcast_to(np.array(0j), (nb, nret))
                    )
                    sphere.append(
                        (
                            a,
                            l,
                            k,
                            PreparedSphericalBand(
                                ch,
                                matrix_ctx[a][l].primitive,
                                matrix_ctx[a][l].lo_coefficients,
                                matrix_ctx[a][l].overlap_metric,
                                matrix_ctx[a][l].radius,
                                matrix_ctx[a][l].energy,
                                matrix_ctx[a][l].lo_energies,
                                triangle,
                                lapw=rc.lapw,
                                cache_products=False,
                                share_radial=True,
                            ),
                        )
                    )
                    if l <= 8:
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
                        z[l] = ch.value(
                            rc.primitive[-1], rc.lo_coefficients, lapw=rc.lapw
                        ).reshape(-1, nb)
                        mc = matrix_ctx[a][l]
                        zks[l] = ch.value(
                            mc.primitive[-1], mc.lo_coefficients, lapw=mc.lapw
                        ).reshape(-1, nb)
                dz = np.concatenate([z[l] for l in range(9)])
                dens.append((k, dz[:, :nret], 1 / s.mult[a]))
                density_matching[a].append(
                    {l: _retained_channel(ch, nret) for l, ch in channels.items()}
                )
                density_complement.setdefault(a, []).append(dz[:, nret:])
                ns.append(
                    (
                        np.concatenate([zks[l] for l in range(lns + 1)]),
                        {l: channels[l] for l in range(lns + 1)},
                        {l: matrix_ctx[a][l].lo_coefficients for l in range(lns + 1)},
                        {l: ctx[a][l].lapw for l in range(lns + 1)},
                    )
                )
                ks.append(k)
        hs = [key[1:] for key in reference.rad if key[0] == a]
        small = 1 / 137.0359895**2 if "RELA" in s.mode else 4e-22
        density_atoms[a] = MuffinTinDensitySum(
            dens, rho, {l: ctx[a][l].primitive[:, :2] for l in range(9)}, hs, small
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
                ),
                ks,
            )
        )
    for a, l, k, op in sphere:
        for name, compute in [
            ("_zpw", lambda op=op: op.zpw),
            ("_kz", lambda op=op: op.kz),
            ("_sz", lambda op=op: op.sz),
        ]:
            if (
                getattr(op, name) is None
                and basis_cache_bytes + op.z.nbytes <= basis_cache_limit
            ):
                array = compute()
                array.setflags(write=False)
                setattr(op, name, array)
                basis_cache_bytes += array.nbytes
    model.basis_cache_bytes = basis_cache_bytes
    model.basis_cache_limit = basis_cache_limit
    radial_projector = {}
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
    for ic, orbs in inl.cix.items():
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
    if not np.allclose(
        np.array([r[:nret] for r in raw]), model.raw, atol=1e-12, rtol=1e-10
    ):
        raise ValueError("assembled raw projectors differ from native snapshot")
    complement = EmptyComplementResponse(
        model, [b.solver for b in bases], [r[nret:] for r in raw]
    )
    reps = np.array(list(reference.pw))
    allg = [r.labels[:, : r.nbasis - r.nlo].T for r in records]
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
            cache_wavefunctions=cache_wavefunctions,
        )
        for (c, *_), g in zip(geometry, allg)
    ]
    spatial = SpatialValenceDensity(reference, density_atoms, ir)
    direct = BasisDensityResponse(
        spatial,
        {a: {l: ctx[a][l] for l in range(9)} for a in ctx},
        density_matching,
        energy_policy="fixed",
    )
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
    compact_spatial = CompactSpatialValenceDensity(
        spatial,
        [b.solver.nq for b in bases],
        density_complement,
        [c[: len(g), nret:] for (c, *_), g in zip(geometry, allg)],
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
            complex_values=complex_vectors,
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
            "assembled overlap disagrees with native H/O; check matrix_quadrature and snapshot consistency"
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
    )
