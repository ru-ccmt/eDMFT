"""Prepared spherical H/O in central bands, including native triangle order.

All radial arrays use the discrete ATPARN convention. The symmetric small
primitive kernel can be contracted directly in bands. The only ordered part
is a PW Wronskian term; its geometry/band matrix is prepared once. Thus no
Npw-by-Npw matrix is built inside a forward or reverse response action.
"""

import numpy as np
from scipy.special import eval_legendre
from lapw_response_scratch import StoredMatrix
from lapw_density_response import checked, herm
from lapw_response_assembly import NativeBandProjection
from lapw_compact_blocks import (
    checked_columns,
    hermitian_columns,
    right_multiply_columns,
    pair_columns,
)


def prepare_wronskian_band(
    ell, rmt, volume, q, phase, coefficients, matching, *, retained=None
):
    """One site's ordered PW matrix T, with the native upper-triangle copy.

    q is GLOBAL Cartesian k+G; phase is exp(i G.R_atom). Matching contains
    j_l(qR) and q*j'_l(qR). Supply the same matching values as the native
    assembler for exact finite-precision comparisons. Only PW rows appear.
    T=C^dagger F[angular*(dj_i*j_j-j_i*dj_j)] C.
    With retained supplied, return only [T_RR; T_QR]. Native coefficient-space
    triangle ordering is applied BEFORE projection, exactly as in the full map.
    """
    q = np.asarray(q)
    c = np.asarray(coefficients)
    if q.ndim != 2 or q.shape[1] != 3 or c.ndim != 2 or len(q) != len(c):
        raise ValueError("invalid PW geometry")
    checked(q, q.shape, "reciprocal vectors", real=True)
    # Site phases enter only as phase_i^* phase_j. For a radial response
    # the same scalar multiplies every site in an orbit, so the projector
    # C and the native triangle operation may act after the site sum.
    phase = np.atleast_2d(np.asarray(phase))
    if phase.ndim != 2 or not len(phase):
        raise ValueError("empty or invalid site phases")
    checked(phase, (len(phase), len(q)), "site phases")
    b = checked(matching, (2, len(q)), "Bessel matching", real=True)
    if (
        int(ell) != ell
        or ell < 0
        or not np.isfinite(rmt)
        or rmt <= 0
        or not np.isfinite(volume)
        or volume <= 0
    ):
        raise ValueError("invalid spherical geometry")
    length = np.linalg.norm(q, axis=1)
    unit = np.divide(
        q, length[:, None], out=np.zeros_like(q, dtype=float), where=length[:, None] > 0
    )
    angular = (
        4
        * np.pi
        / volume
        * (2 * ell + 1)
        * eval_legendre(ell, np.clip(unit @ unit.T, -1, 1))
    )
    angular = angular * (phase.conj().T @ phase)
    radial = b[1, :, None] * b[0, None, :] - b[0, :, None] * b[1, None, :]
    if retained is None:
        return herm(NativeBandProjection(angular, c).project(radial))
    if int(retained) != retained or not 0 < retained <= c.shape[1]:
        raise ValueError("invalid retained Wronskian count")
    upper = np.triu(angular * radial, 1)
    native = upper + upper.conj().T
    np.fill_diagonal(native, np.diag(angular * radial).real)
    return hermitian_columns(c.conj().T @ (native @ c[:, : int(retained)]))


class PreparedSphericalBand:
    """One atom/site/l/k spherical contribution and its exact transpose.

    A square Wronskian kernel preserves the original full-matrix interface.
    A (Nband,Nretained) kernel makes value/action return only RR/QR columns;
    pullback then uses the physical QR+RQ pairing, with no QQ allocation.

    cache_products=False recomputes small primitive-by-band products during
    an action instead of storing KZ, SZ and a second Z for every k/site/l.
    share_radial=True shares immutable native radial arrays across k and sites.
    Neither option changes the forward map or its exact real transpose.

    channel is PreparedAugmentationChannel in the GLOBAL spherical frame;
    local rotations are unnecessary for a spherical kernel. primitive and
    metric are native (radius,4,radial_slot) and (radius,4), respectively.
    The physical-minus-basis potential is r*(Vphysical-Vbasis) in Ry. It is
    zero for a common physical/basis potential, including its tangent; a
    frozen-physical-H test instead supplies difference_dot=-basis_potential_dot.

    Native H sets the LO boundary VALUE to zero but retains its boundary
    SLOPE. The ordered PW correction is 0.5*R^2*W*(W-1)*T for LAPW and zero
    for APW. W=AP_u*A_udot-A_u*AP_udot. It must not be discarded or averaged
    before the native triangle operation, even though W is close to one.
    """

    def __init__(
        self,
        channel,
        primitive,
        lo_coefficients,
        metric,
        radius,
        energy,
        lo_energies,
        wronskian_band,
        *,
        lapw,
        difference=None,
        cache_products=True,
        share_radial=False,
        site_triangle=None,
        retained=None,
    ):
        self.cache_products = bool(cache_products)

        def stored(value, shape, name):
            a = checked(value, shape, name, real=True)
            if share_radial:
                # Only immutable native-context arrays may be shared. Never
                # freeze caller-owned writable inputs as an accidental side effect.
                if a.flags.writeable or a.dtype != np.float64:
                    raise ValueError("shared radial arrays must be read-only float64")
                return a
            return np.array(a, dtype=float, copy=True)

        self.channel = channel
        self.lapw = bool(lapw)
        self.r = channel.rmt
        array = np.asarray(primitive)
        if array.ndim != 3:
            raise ValueError("invalid primitive dimensions")
        nr = len(array)
        m = 2 + channel.nlo
        self.phi = stored(array, (nr, 4, m), "primitive")
        self.metric = stored(metric, (nr, 4), "radial metric")
        self.radius = stored(radius, (nr,), "radial mesh")
        if np.any(self.radius <= 0):
            raise ValueError("nonpositive radial mesh")
        self.lo = stored(lo_coefficients, (m, channel.nlo), "LO factors")
        energies = np.r_[
            energy,
            energy,
            checked(lo_energies, (channel.nlo,), "LO energies", real=True),
        ]
        checked(energies, (m,), "radial energies", real=True)
        self.t = np.diag(energies)
        self.t[0, 1] = 1.0
        self.difference = (
            np.zeros(nr)
            if difference is None
            else np.array(
                checked(difference, (nr,), "potential difference", real=True), copy=True
            )
        )
        # A rectangular kernel requests only RR and QR output. This retains
        # the native triangle correction but never stores its unused QQ part.
        self._site_triangle = site_triangle
        if wronskian_band is None:
            # The first site carries the whole orbit matrix. The other
            # sites reconstruct their own triangle only for a site query.
            if site_triangle is None or retained is None or not 0 < retained <= channel.nband:
                raise ValueError("missing Wronskian matrix or orbit metadata")
            self._orbit_triangle = None
            self.band_shape = (channel.nband, retained)
        else:
            triangle = checked_columns(wronskian_band, channel.nband, "Wronskian band kernel")
            self.band_shape = triangle.shape
            self._orbit_triangle = (np.array(triangle, complex) if self.lapw else
                                    np.broadcast_to(np.array(0j), self.band_shape))
        self.retained = self.band_shape[1]
        self.s = self._gram(self.phi, self.phi)
        self.s[0, 0] = 1.0
        self.s[0, 1] = self.s[1, 0] = 0.0
        self.k = (self.s @ self.t + self.t.T @ self.s) / 2 + self._gram(
            self.phi, self.phi, self.difference / self.radius
        )
        # Thousands of (site,l,k) objects can exist. Z, p and d are
        # inexpensive radial combinations of F,D and the LO matching sums.
        # In the bounded-cache path form them once per call and release them
        # afterwards. The immutable radial arrays and matching sums suffice.
        self._z = channel.value(self.phi[-1], self.lo, lapw=lapw) if cache_products else None
        self._zpw = None
        if cache_products:
            self._zpw = (self._z if channel.nlo == 0 else
                         channel.value(self.phi[-1], np.zeros_like(self.lo), lapw=lapw))
        self.w = (
            self.phi[-1, 2, 0] * self.phi[-1, 0, 1]
            - self.phi[-1, 0, 0] * self.phi[-1, 2, 1]
        )
        self._kz = (
            np.einsum("ij,mjn->min", self.k, self.z) if self.cache_products else None
        )
        self._sz = (
            np.einsum("ij,mjn->min", self.s, self.z) if self.cache_products else None
        )

    def _wronskian_triangle(self, orbit_sum):
        """Radial directions and their adjoints share one coefficient per orbit.

        For those contractions the first site carries the preprojected sum.
        A site-resolved value uses the original individual PW phases.
        Explicit geometry derivatives use each site's own phase derivative.
        """
        if self._site_triangle is not None and not orbit_sum:
            return self._site_triangle()
        if isinstance(self._orbit_triangle, StoredMatrix):
            return self._orbit_triangle.read()
        return self._orbit_triangle

    def store_matching(self, scratch):
        self.channel.store_matching(scratch)
        if self.lapw and self._orbit_triangle is not None:
            self._orbit_triangle = scratch.store(self._orbit_triangle)

    @property
    def triangle(self):
        return self._wronskian_triangle(False)

    @property
    def z(self):
        return self._z if self._z is not None else self.channel.value(
            self.phi[-1], self.lo, lapw=self.lapw
        )

    def _boundary_value(self, zpw):
        return np.einsum("i,min->mn", self.phi[-1, 0] / self.r, zpw)

    def _boundary_slope(self, z):
        return np.einsum("i,min->mn", self.phi[-1, 2] / self.r, z)

    @property
    def p(self):
        return self._boundary_value(self.zpw)

    @property
    def d(self):
        return self._boundary_slope(self.z)

    def _matching(self):
        z = self.z
        return z, z if self.channel.nlo == 0 else self.zpw

    def _products(self, z):
        kz = self._kz if self._kz is not None else np.einsum("ij,mjn->min", self.k, z)
        sz = self._sz if self._sz is not None else np.einsum("ij,mjn->min", self.s, z)
        return kz, sz

    @property
    def kz(self):
        return (
            self._kz
            if self._kz is not None
            else np.einsum("ij,mjn->min", self.k, self.z)
        )

    @property
    def sz(self):
        return (
            self._sz
            if self._sz is not None
            else np.einsum("ij,mjn->min", self.s, self.z)
        )

    @property
    def zpw(self):
        return (
            self._zpw
            if self._zpw is not None
            else self.channel.value(
                self.phi[-1], np.zeros_like(self.lo), lapw=self.lapw
            )
        )

    def _gram(self, a, b, potential=None):
        weight = self.metric if potential is None else self.metric * potential[:, None]
        return np.einsum("rci,rc,rcj->ij", a, weight, b, optimize=True)

    def value(self, *, orbit_sum=False):
        n = self.retained
        z, zpw = self._matching()
        p, d = self._boundary_value(zpw), self._boundary_slope(z)
        kz, sz = self._products(z)
        h = np.einsum("min,mip->np", z.conj(), kz[..., :n], optimize=True)
        o = np.einsum("min,mip->np", z.conj(), sz[..., :n], optimize=True)
        h += (
            self.r**2
            / 2
            * (p.conj().T @ d[:, :n] + d.conj().T @ p[:, :n])
        )
        triangle = self._wronskian_triangle(orbit_sum) if self.lapw else None
        if triangle is not None:
            h += 0.5 * self.r**2 * self.w * (self.w - 1) * triangle
        return hermitian_columns(h), hermitian_columns(o)

    def action(
        self,
        primitive_dot,
        lo_dot,
        energy_dot=0.0,
        lo_energies_dot=None,
        difference_dot=None,
        *, orbit_sum=False,
    ):
        u = checked(primitive_dot, self.phi.shape, "primitive derivative", real=True)
        dm = checked(lo_dot, self.lo.shape, "LO factor derivative", real=True)
        de = float(checked(energy_dot, (), "energy derivative", real=True))
        delocal = (
            np.zeros(self.channel.nlo)
            if lo_energies_dot is None
            else checked(
                lo_energies_dot, (self.channel.nlo,), "LO energy derivative", real=True
            )
        )
        dv = (
            np.zeros(len(self.phi))
            if difference_dot is None
            else checked(
                difference_dot,
                (len(self.phi),),
                "potential difference derivative",
                real=True,
            )
        )
        ds = self._gram(u, self.phi) + self._gram(self.phi, u)
        ds[0, 0] = ds[0, 1] = ds[1, 0] = 0.0
        dt = np.diag(np.r_[de, de, delocal])
        dk = (ds @ self.t + self.s @ dt + dt.T @ self.s + self.t.T @ ds) / 2
        dk += (
            self._gram(u, self.phi, self.difference / self.radius)
            + self._gram(self.phi, u, self.difference / self.radius)
            + self._gram(self.phi, self.phi, dv / self.radius)
        )
        z, zpw = self._matching()
        p, d = self._boundary_value(zpw), self._boundary_slope(z)
        kz, sz = self._products(z)
        dz = self.channel.action(self.phi[-1], u[-1], dm, lapw=self.lapw)
        dzpw = self.channel.action(
            self.phi[-1], u[-1], np.zeros_like(dm), lapw=self.lapw
        )
        # Keep the two distinct differentiated legs before restricting the
        # right band index. Hermitian transposition of rectangular QR would
        # otherwise silently lose one of these terms.
        n = self.retained
        h = (
            np.einsum("min,mip->np", dz.conj(), kz[..., :n], optimize=True)
            + np.einsum("min,mip->np", kz.conj(), dz[..., :n], optimize=True)
            + np.einsum(
                "min,ij,mjp->np", z.conj(), dk, z[..., :n], optimize=True
            )
        )
        o = (
            np.einsum("min,mip->np", dz.conj(), sz[..., :n], optimize=True)
            + np.einsum("min,mip->np", sz.conj(), dz[..., :n], optimize=True)
            + np.einsum(
                "min,ij,mjp->np", z.conj(), ds, z[..., :n], optimize=True
            )
        )
        pd = (
            np.einsum("i,min->mn", u[-1, 0], zpw)
            + np.einsum("i,min->mn", self.phi[-1, 0], dzpw)
        ) / self.r
        dd = (
            np.einsum("i,min->mn", u[-1, 2], z)
            + np.einsum("i,min->mn", self.phi[-1, 2], dz)
        ) / self.r
        h += (
            self.r**2
            / 2
            * (
                pd.conj().T @ d[:, :n]
                + p.conj().T @ dd[:, :n]
                + d.conj().T @ pd[:, :n]
                + dd.conj().T @ p[:, :n]
            )
        )
        dw = (
            u[-1, 2, 0] * self.phi[-1, 0, 1]
            + self.phi[-1, 2, 0] * u[-1, 0, 1]
            - u[-1, 0, 0] * self.phi[-1, 2, 1]
            - self.phi[-1, 0, 0] * u[-1, 2, 1]
        )
        triangle = self._wronskian_triangle(orbit_sum) if self.lapw and dw != 0 else None
        if triangle is not None:
            h += 0.5 * self.r**2 * (2 * self.w - 1) * dw * triangle
        return hermitian_columns(h), hermitian_columns(o)

    def physical_potential_pullback(self, hbar):
        """Contract the explicit potential leg at fixed basis and radial functions.

        Only Z† delta(K) Z varies on this leg, with
        delta(K)_ij = integral phi_i phi_j delta(V) dr (native metric).
        The ordinary radial pullback includes the same quantity as its
        final output, plus primitive and boundary terms for other directions.
        """
        hbar = checked_columns(hbar, self.channel.nband, "physical potential covector")
        if hbar.shape != self.band_shape:
            raise ValueError("physical covector column layout mismatch")
        z = self.z
        zh = right_multiply_columns(z, hbar)
        gk = np.einsum("min,mjn->ij", zh, z.conj(), optimize=True).real
        return self._potential_covector(gk)

    def _potential_covector(self, gk):
        return np.einsum("rci,ij,rcj,rc->r", self.phi, gk, self.phi,
                         self.metric, optimize=True) / self.radius

    def pullback(self, hbar, obar, *, orbit_sum=False):
        hbar = checked_columns(hbar, self.channel.nband, "Hamiltonian covector")
        obar = checked_columns(obar, self.channel.nband, "overlap covector")
        if hbar.shape != self.band_shape or obar.shape != self.band_shape:
            raise ValueError("spherical covector column layout mismatch")
        # Hbar and Obar act only in band space; the small radial matrices
        # k and s act in primitive space. Associate their products so each
        # expensive RR/QR band contraction is performed once. The same zh
        # also supplies the boundary-slope contraction below. This changes
        # no derivative and stores no extra per-k cache.
        z, zpw = self._matching()
        p = self._boundary_value(zpw)
        zh = right_multiply_columns(z, hbar)
        zo = right_multiply_columns(z, obar)
        zbar = 2 * (
            np.einsum("ij,mjn->min", self.k, zh) + np.einsum("ij,mjn->min", self.s, zo)
        )
        gk = np.einsum("min,mjn->ij", zh, z.conj(), optimize=True).real
        gs = (
            np.einsum("min,mjn->ij", zo, z.conj(), optimize=True).real
            + (gk @ self.t.T + self.t @ gk) / 2
        )
        gt = self.s @ gk
        gs[0, 0] = gs[0, 1] = gs[1, 0] = 0.0
        gp = 2 * self.metric[:, :, None] * np.einsum("rcj,ji->rci", self.phi, gs)
        gp += (
            2
            * (self.metric * (self.difference / self.radius)[:, None])[:, :, None]
            * np.einsum("rcj,ji->rci", self.phi, gk)
        )
        difference_bar = self._potential_covector(gk)
        # d = sum_i phi_slope[i]*z[i]/R, hence d@Hbar is already in zh.
        pb = self.r * np.einsum("i,min->mn", self.phi[-1, 2], zh)
        db = self.r**2 * right_multiply_columns(p, hbar)
        gp[-1, 0] += np.einsum("mn,min->i", pb.conj(), zpw).real / self.r
        gp[-1, 2] += np.einsum("mn,min->i", db.conj(), z).real / self.r
        zbar += self.phi[-1, 2][None, :, None] * db[:, None, :] / self.r
        zpwbar = self.phi[-1, 0][None, :, None] * pb[:, None, :] / self.r
        gb, gm = self.channel.pullback(self.phi[-1], zbar, lapw=self.lapw)
        gbpw, _ = self.channel.pullback(self.phi[-1], zpwbar, lapw=self.lapw)
        gp[-1] += gb + gbpw
        triangle = self._wronskian_triangle(orbit_sum) if self.lapw else None
        if triangle is not None:
            wb = 0.5 * self.r**2 * (2 * self.w - 1) * pair_columns(hbar, triangle)
            gp[-1, 2, 0] += wb * self.phi[-1, 0, 1]
            gp[-1, 0, 1] += wb * self.phi[-1, 2, 0]
            gp[-1, 0, 0] -= wb * self.phi[-1, 2, 1]
            gp[-1, 2, 1] -= wb * self.phi[-1, 0, 0]
        return (
            gp,
            gm,
            float(gt[0, 0] + gt[1, 1]),
            np.diag(gt)[2:].copy(),
            difference_bar,
        )
