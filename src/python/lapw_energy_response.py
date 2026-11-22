"""Small-rank completion for successful native radial-energy searches.

The ordinary response holds radial energies fixed. SELECT's searched energies
are extra scalar inputs e=A V_sph. Compute each independent energy column once:
P_i = N [R E B_energy,i + D_energy,i], g_i = g_lattice.B_energy,i.
Then J_total=J_fixed+P A K and c_total=c_fixed+K^T A^T g. The source and
explicit scalar get the matching P A V_exp and g.A V_exp partners. This
retains their cancellations and adds no independent fitted stress term.
"""

import numpy as np
from response_parallel import parallel_of
from lapw_density_response import checked
from lapw_radial_energy_policy import searched_energy_derivatives
from stress_response_adjoint import PairedOperator, ScreenedResponse


def energy_direction(valence, energies):
    """Joint band and direct spatial derivatives at fixed physical potential."""
    if hasattr(valence, "branches"):
        pieces = [energy_direction(v, energies) for v in valence.branches]
        return valence.bands.output.pack(*(x[0] for x in pieces)), .5 * sum(x[1] for x in pieces)
    op = valence.bands.spherical
    c = op.coordinates
    nk, nr = c.raw_shape[:2]
    t = {}
    pt = {}
    for (a, l), de in energies.items():
        ctx = op.contexts[a][l]
        de = checked(de, (1 + ctx.nlo,), "radial energy direction", real=True)
        t[a, l] = ctx.action(np.zeros(ctx.n), de[0], de[1:])
        ctx = op.projector_contexts[a][l]
        pt[a, l] = ctx.action(np.zeros(ctx.n), de[0], de[1:])
    h = np.zeros((nk, nr, nr), complex)
    o = np.zeros_like(h)
    hq = [np.zeros(shape, complex) for shape in c.matrix.shapes]
    oq = [np.zeros_like(v) for v in hq]
    raw = np.zeros(c.raw_shape, complex)
    for a, l, k, s in op.spherical:
        if (a, l) not in energies:
            continue
        de = energies[a, l]
        dh, do = s.action(*t[a, l], energy_dot=de[0], lo_energies_dot=de[1:], orbit_sum=True)
        h[k] += dh[:nr]
        hq[k] += dh[nr:]
        o[k] += do[:nr]
        oq[k] += do[nr:]
    for a, ns, ks in op.nonspherical:
        if not any(aa == a for aa, l in energies):
            continue
        dp = {
            l: t[a, l][0] if (a, l) in t else np.zeros_like(op.contexts[a][l].primitive)
            for l in ns.phi
        }
        dm = {
            l: (
                t[a, l][1]
                if (a, l) in t
                else np.zeros_like(op.contexts[a][l].lo_coefficients)
            )
            for l in ns.phi
        }
        for k, dh in zip(ks, ns.action(dp, dm, retained=nr)):
            h[k] += dh[:nr]
            hq[k] += dh[nr:]
    for a, l, k, ic, cols, p in op.projectors:
        if (a, l) in pt:
            raw[k, :, cols, ic] += p.action(*pt[a, l])[:nr]
    direct = valence.direct
    dp = {a: None for a in direct.contexts}
    dz = {a: None for a in direct.contexts}
    for a in direct.contexts:
        if not any(aa == a for aa, l in energies):
            continue
        dp[a] = {
            l: (
                pt[a, l][0][:, :2]
                if (a, l) in pt
                else np.zeros_like(ctx.primitive[:, :2])
            )
            for l, ctx in direct.contexts[a].items()
        }
        dz[a] = []
        for channels in direct.matching[a]:
            blocks = []
            for l, ch in sorted(channels.items()):
                ctx = direct.contexts[a][l]
                if (a, l) in pt:
                    p, m = pt[a, l]
                    block = ch.action(ctx.primitive[-1], p[-1], m, lapw=ctx.lapw)
                else:
                    block = np.zeros((2 * l + 1, 2 + ctx.nlo, ch.nband), complex)
                blocks.append(block.reshape(-1, ch.nband))
            dz[a].append(np.concatenate(blocks))
    density = direct.spatial.action(matching_dots=dz, primitive_dots=dp)
    return c.pack(h, o, raw, hq, oq), density


class SearchedEnergyResponse:
    """Attach searched radial energies, holding the IN1 reference EF fixed.

    Fixed RMT means no separate explicit geometry of the radial root. A
    changing IN1 reference EF, a root/fallback switch, HDLO, or moving RMT
    requires a different policy and is not silently supplied here.
    """

    def __init__(
        self,
        state,
        observable,
        directory,
        case,
        *,
        complex_vectors,
        library,
        progress=None
    ):
        parallel = parallel_of(state)
        maps, details = searched_energy_derivatives(
            directory, case, complex_vectors=complex_vectors, library=library
        )
        self.state = state
        self.observable = observable
        self.details = details
        self.A = np.zeros((len(details), state.chain.potential.output_size))
        columns = []
        direct = []
        self.search_bands = []
        for i, item in enumerate(details):
            a, l = item["atom"], item["l"]
            slots = item["slots"]
            mat = maps[a, l]
            sl = next(
                sl
                for key, sl in zip(
                    state.potential_core.hartree.keys,
                    state.potential_core.hartree.slices,
                )
                if key == (a, 0, 0)
            )
            self.A[i, sl] = mat[slots[0]]
            de = np.zeros(len(mat))
            de[slots] = 1
            band, density = energy_direction(state.valence, {(a, l): de})
            if getattr(state.valence.model, "electronic_model", "dmft") == "dft":
                self.search_bands.append(band)
            ch = state.chain
            columns.append(
                ch.normalization.matvec(
                    ch.density.matvec(ch.electronic.matvec(band))
                    + parallel.sum(density)
                )
            )
            direct.append(parallel.sum(observable.lattice.joint_gradient @ band))
            if progress is not None:
                progress(i, item)
        self.P = (
            np.column_stack(columns)
            if columns
            else np.zeros((state.chain.normalization.output_size, 0))
        )
        self.g = np.array(direct)
        self.K = state.chain.potential
        self.J = state.jacobian
        # Factor on the common potential input: K (especially PBE) is
        # expensive and must be applied only once per forward/reverse action.
        fixed = state.chain.charge_from_potential
        self.charge_from_potential = PairedOperator(
            fixed.input_size,
            fixed.output_size,
            lambda v: fixed.matvec(v) + self.P @ (self.A @ v),
            lambda y: fixed.rmatvec(y) + self.A.T @ (self.P.T @ y),
            "charge response including searched radial energies",
        )
        self.jacobian = self.charge_from_potential @ self.K
        self.screened = ScreenedResponse(self.jacobian)
        self.charge_gradient = observable.charge_gradient + self.K.rmatvec(
            self.A.T @ self.g
        )

    def complete_source(self, fixed_source, explicit_potential):
        return fixed_source + self.P @ (self.A @ explicit_potential)

    def explicit_addition(self, explicit_potential):
        return self.g @ (self.A @ explicit_potential)

    def total_energy_directions(self, explicit_potential, density_response):
        return self.A @ (explicit_potential + self.K.matvec(density_response))

def in1_reference_fermi_slots(path, snapshot):
    """Return ATPAR's EF-dependent radial-energy slots from literal IN1 input.

    Native ATPAR maps each 0.30 input to EF +/- 0.20 Ry.  APW+lo creates
    an LO at the first occurrence of a channel; subsequent same-l entries
    create further LOs.  HDLO (iapw=2) inherits its parent energy.  SELECT
    search roots are removed by the caller after this input-level traversal.
    This parser follows the channel and LO ordering, not a proximity guess
    based on the final energies.
    """
    from pathlib import Path
    lines = Path(path).read_text().splitlines()
    pos = 2  # WFFIL/EF header and RKMAX line
    result = {}
    for a in range(snapshot.s.nat):
        words = lines[pos].split()
        if len(words) < 3:
            raise ValueError('incomplete IN1 atom header')
        default = abs(float(words[0]) - 0.30) < 1e-6
        count = int(words[1])
        pos += 1
        main = [default] * len(snapshot.header[a])
        local = {}
        previous = -1
        for unused in range(count):
            words = lines[pos].split()
            pos += 1
            if len(words) < 2:
                raise ValueError('incomplete IN1 channel')
            l, energy = int(words[0]), float(words[1])
            if l < 0 or l >= len(main):
                raise ValueError('IN1 channel outside radial header')
            apw = int(words[4]) if len(words) > 4 else 0
            flag = abs(energy - 0.30) < 1e-6
            if apw == 2:
                flag = local.get(l, [main[l]])[0]
            if l == previous:
                local.setdefault(l, []).append(flag)
            else:
                main[l] = flag
                if apw == 1:
                    local[l] = [flag]
            previous = l
        for l in range(len(main)):
            result[a, l] = np.array([main[l], *local.get(l, [])], dtype=float)
    if not lines[pos].lstrip().startswith('K-VECTORS'):
        raise ValueError('IN1 atom block parse did not reach K-VECTORS')
    return result

class DFTReferenceEnergyResponse(SearchedEnergyResponse):
    """Differentiate the IN1 reference EF as a second SCF variable.

    WIEN LAPW1 replaces an IN1 energy of 0.30 by EF +/- 0.20 Ry.  The
    energy-search roots above have their own potential derivative, while these
    unsearched energies move with the converged Fermi level.  Let q be the
    material charge and e the change of the IN1 EF.  Their joint equations are

        q = s_q + J q + p e,
        e = s_mu + m_q q + m_e e.

    The same scalar e changes every eligible radial channel.  The second row
    is the fixed-N Mermin chemical-potential derivative, not an independent
    fit to strained endpoint calculations.  One augmented adjoint still serves
    all strain components.
    """

    def __init__(self, state, observable, directory, case, *, complex_vectors,
                 library, progress=None):
        super().__init__(state, observable, directory, case,
                         complex_vectors=complex_vectors, library=library,
                         progress=progress)
        import re
        from pathlib import Path
        from lapw_valence_state import _snapshot_inputs
        from response_parallel import parallel_of

        if getattr(state.valence.model, 'electronic_model', None) != 'dft':
            raise ValueError('reference-EF response requires Mermin DFT')
        if not state.fixed_number:
            raise ValueError('reference-EF response requires fixed particle number')
        directory = Path(directory)
        suffix = '.in1c' if complex_vectors else '.in1'
        header = (directory / (case + suffix)).read_text().splitlines()[0]
        match = re.search(r'EF=\s*([-+0-9.Ee]+)', header)
        if match is None:
            raise ValueError('IN1 reference EF is missing')
        reference = float(match[1])
        mu = float(state.valence.model.meta.mu)
        if abs(reference - mu) > 5e-5:
            raise ValueError('IN1 reference EF differs from the converged Fermi energy')
        m = _snapshot_inputs(case, directory, complex_vectors, dft=True)
        selected = {(d['atom'], d['l']): set(d['slots']) for d in self.details}
        input_flags = in1_reference_fermi_slots(directory / (case + suffix), m)
        energies = {}
        slots_report = []
        for a, channels in state.valence.bands.spherical.contexts.items():
            for l, ctx in channels.items():
                main = float(m.header[a][l])
                if main > 150:
                    main -= 200
                lo = m.lo[a][l, :m.ilo[a][l]] if l < 4 else np.empty(0)
                values = np.r_[main, lo]
                flags = input_flags[a, l].copy()
                if len(flags) != 1 + ctx.nlo or len(values) != len(flags):
                    raise ValueError('IN1/reference radial energy layout differs')
                for j in selected.get((a, l), ()):
                    flags[j] = 0
                for j in np.flatnonzero(flags):
                    v = values[j]
                    if min(abs(v-reference-.2), abs(v-reference+.2)) >= 3e-5:
                        raise ValueError('IN1 EF-relative energy differs from native vector')
                if np.any(flags):
                    energies[a, l] = flags
                    slots_report.append(dict(atom=a, l=l,
                                             slots=np.flatnonzero(flags).tolist()))
        if not energies:
            raise ValueError('no EF-relative native radial channels found')
        # Validate every literal 0.30-derived input against the native vector.
        self.reference_slots = slots_report
        self.reference_EF_Ry = reference
        band, direct_density = energy_direction(state.valence, energies)
        ch = state.chain
        parallel = parallel_of(state)
        self.reference_charge = ch.normalization.matvec(
            ch.density.matvec(ch.electronic.matvec(band))
            + parallel.sum(direct_density)
        )
        self.reference_gradient = float(parallel.sum(observable.lattice.joint_gradient @ band))
        model = state.valence.model
        weights = 2 * model.kw[:, None] * model._slope / model._charge_slope
        hbar = np.array([np.diag(row) for row in weights], complex)
        obar = -hbar * model.eps[:, :, None]
        output = state.valence.bands.output
        qr = [np.zeros(shape, complex) for shape in output.matrix.shapes]
        self.mu_band_covector = output.pack(
            hbar, obar, np.zeros(output.raw_shape, complex), qr, qr
        )
        self.mu_reference = self.mu_from_band(band)
        mu_search = np.array([
            self.mu_from_band(v) for v in self.search_bands
        ])
        del self.search_bands
        self.mu_potential_covector = (
            ch.basis.rmatvec(self.mu_band_covector) + self.A.T @ mu_search
        )
        self.mu_charge_covector = self.K.rmatvec(self.mu_potential_covector)
        base_jacobian = self.jacobian
        base_gradient = self.charge_gradient
        p = self.reference_charge
        m_q = self.mu_charge_covector
        m_e = self.mu_reference
        size = len(p)
        self.jacobian = PairedOperator(
            size+1, size+1,
            lambda z: np.r_[base_jacobian.matvec(z[:-1])+p*z[-1],
                            m_q@z[:-1]+m_e*z[-1]],
            lambda y: np.r_[base_jacobian.rmatvec(y[:-1])+m_q*y[-1],
                            p@y[:-1]+m_e*y[-1]],
            'charge and IN1 reference-EF fixed point',
        )
        self.screened = ScreenedResponse(self.jacobian)
        self.charge_gradient = np.r_[base_gradient, self.reference_gradient]
        self.reference_report = dict(
            in1_reference_Ry=reference, central_mu_Ry=mu,
            radial_slots=slots_report,
            direct_energy_gradient_Ry_per_Ry=self.reference_gradient,
            charge_source_norm=float(np.linalg.norm(p)),
            mu_reference_derivative=m_e,
            mu_charge_covector_norm=float(np.linalg.norm(m_q)),
        )

    def mu_from_band(self, band):
        from response_parallel import parallel_of
        return float(parallel_of(self.state).sum(self.mu_band_covector @ band))

    def complete_source(self, fixed_source, explicit_potential, band_geometry=None):
        if band_geometry is None:
            raise ValueError('DFT reference-EF source needs band geometry')
        density = super().complete_source(fixed_source, explicit_potential)
        chem = (self.mu_from_band(band_geometry)
                + self.mu_potential_covector @ explicit_potential)
        return np.r_[density, chem]

    def source_mu(self, band_geometry, explicit_potential):
        return (self.mu_from_band(band_geometry)
                + self.mu_potential_covector @ explicit_potential)
