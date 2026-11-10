"""Prepared spherical H/O in central bands, including native triangle order.

All radial arrays use the discrete ATPARN convention. The symmetric small
primitive kernel can be contracted directly in bands. The only ordered part
is a PW Wronskian term; its geometry/band matrix is prepared once. Thus no
Npw-by-Npw matrix is built inside a forward or reverse response action.
"""

import numpy as np
from scipy.special import eval_legendre
from lapw_density_response import checked, herm
from lapw_response_assembly import NativeBandProjection


def prepare_wronskian_band(ell, rmt, volume, q, phase, coefficients, matching):
    """One site's ordered PW matrix T, with the native upper-triangle copy.

    q is GLOBAL Cartesian k+G; phase is exp(i G.R_atom). Matching contains
    j_l(qR) and q*j'_l(qR). Supply the same matching values as the native
    assembler for exact finite-precision comparisons. Only PW rows appear.
    T=C^dagger F[angular*(dj_i*j_j-j_i*dj_j)] C.
    """
    q = np.asarray(q)
    c = np.asarray(coefficients)
    if q.ndim != 2 or q.shape[1] != 3 or c.ndim != 2 or len(q) != len(c):
        raise ValueError("invalid PW geometry")
    checked(q, q.shape, "reciprocal vectors", real=True)
    phase = checked(phase, (len(q),), "site phase")
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
    angular = angular * phase.conj()[:, None] * phase[None, :]
    radial = b[1, :, None] * b[0, None, :] - b[0, :, None] * b[1, None, :]
    return herm(NativeBandProjection(angular, c).project(radial))


class PreparedSphericalBand:
    """One atom/site/l/k spherical contribution and its exact transpose.

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
        difference=None
    ):
        self.channel = channel
        self.lapw = bool(lapw)
        self.r = channel.rmt
        array = np.asarray(primitive)
        if array.ndim != 3:
            raise ValueError("invalid primitive dimensions")
        nr = len(array)
        m = 2 + channel.nlo
        self.phi = np.array(
            checked(array, (nr, 4, m), "primitive", real=True), dtype=float, copy=True
        )
        self.metric = np.array(
            checked(metric, (nr, 4), "radial metric", real=True), copy=True
        )
        self.radius = np.array(
            checked(radius, (nr,), "radial mesh", real=True), copy=True
        )
        if np.any(self.radius <= 0):
            raise ValueError("nonpositive radial mesh")
        self.lo = np.array(
            checked(lo_coefficients, (m, channel.nlo), "LO factors", real=True),
            copy=True,
        )
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
        self.triangle = np.array(
            checked(
                wronskian_band,
                (channel.nband,) * 2,
                "Wronskian band kernel",
                hermitian=True,
            ),
            complex,
        )
        self.s = self._gram(self.phi, self.phi)
        self.s[0, 0] = 1.0
        self.s[0, 1] = self.s[1, 0] = 0.0
        self.k = (self.s @ self.t + self.t.T @ self.s) / 2 + self._gram(
            self.phi, self.phi, self.difference / self.radius
        )
        self.z = channel.value(self.phi[-1], self.lo, lapw=lapw)
        self.zpw = channel.value(self.phi[-1], np.zeros_like(self.lo), lapw=lapw)
        self.p = np.einsum("i,min->mn", self.phi[-1, 0] / self.r, self.zpw)
        self.d = np.einsum("i,min->mn", self.phi[-1, 2] / self.r, self.z)
        self.w = (
            self.phi[-1, 2, 0] * self.phi[-1, 0, 1]
            - self.phi[-1, 0, 0] * self.phi[-1, 2, 1]
        )
        self.kz = np.einsum("ij,mjn->min", self.k, self.z)
        self.sz = np.einsum("ij,mjn->min", self.s, self.z)

    def _gram(self, a, b, potential=None):
        weight = self.metric if potential is None else self.metric * potential[:, None]
        return np.einsum("rci,rc,rcj->ij", a, weight, b, optimize=True)

    def value(self):
        h = np.einsum("min,mip->np", self.z.conj(), self.kz)
        o = np.einsum("min,mip->np", self.z.conj(), self.sz)
        h += self.r**2 / 2 * (self.p.conj().T @ self.d + self.d.conj().T @ self.p)
        if self.lapw:
            h += 0.5 * self.r**2 * self.w * (self.w - 1) * self.triangle
        return herm(h), herm(o)

    def action(
        self,
        primitive_dot,
        lo_dot,
        energy_dot=0.0,
        lo_energies_dot=None,
        difference_dot=None,
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
        dz = self.channel.action(self.phi[-1], u[-1], dm, lapw=self.lapw)
        dzpw = self.channel.action(
            self.phi[-1], u[-1], np.zeros_like(dm), lapw=self.lapw
        )
        a = np.einsum("min,mip->np", dz.conj(), self.kz)
        b = np.einsum("min,mip->np", dz.conj(), self.sz)
        h = a + a.conj().T + np.einsum("min,ij,mjp->np", self.z.conj(), dk, self.z)
        o = b + b.conj().T + np.einsum("min,ij,mjp->np", self.z.conj(), ds, self.z)
        pd = (
            np.einsum("i,min->mn", u[-1, 0], self.zpw)
            + np.einsum("i,min->mn", self.phi[-1, 0], dzpw)
        ) / self.r
        dd = (
            np.einsum("i,min->mn", u[-1, 2], self.z)
            + np.einsum("i,min->mn", self.phi[-1, 2], dz)
        ) / self.r
        a = pd.conj().T @ self.d + self.p.conj().T @ dd
        h += self.r**2 / 2 * (a + a.conj().T)
        dw = (
            u[-1, 2, 0] * self.phi[-1, 0, 1]
            + self.phi[-1, 2, 0] * u[-1, 0, 1]
            - u[-1, 0, 0] * self.phi[-1, 2, 1]
            - self.phi[-1, 0, 0] * u[-1, 2, 1]
        )
        if self.lapw:
            h += 0.5 * self.r**2 * (2 * self.w - 1) * dw * self.triangle
        return herm(h), herm(o)

    def pullback(self, hbar, obar):
        hbar = checked(
            hbar, self.triangle.shape, "Hamiltonian covector", hermitian=True
        )
        obar = checked(obar, self.triangle.shape, "overlap covector", hermitian=True)
        zbar = 2 * (self.kz @ hbar + self.sz @ obar)
        gk = np.einsum("min,np,mjp->ij", self.z, hbar, self.z.conj()).real
        gs = (
            np.einsum("min,np,mjp->ij", self.z, obar, self.z.conj()).real
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
        difference_bar = (
            np.einsum(
                "rci,ij,rcj,rc->r", self.phi, gk, self.phi, self.metric, optimize=True
            )
            / self.radius
        )
        pb = self.r**2 * self.d @ hbar
        db = self.r**2 * self.p @ hbar
        gp[-1, 0] += np.einsum("mn,min->i", pb.conj(), self.zpw).real / self.r
        gp[-1, 2] += np.einsum("mn,min->i", db.conj(), self.z).real / self.r
        zbar += self.phi[-1, 2][None, :, None] * db[:, None, :] / self.r
        zpwbar = self.phi[-1, 0][None, :, None] * pb[:, None, :] / self.r
        gb, gm = self.channel.pullback(self.phi[-1], zbar, lapw=self.lapw)
        gbpw, _ = self.channel.pullback(self.phi[-1], zpwbar, lapw=self.lapw)
        gp[-1] += gb + gbpw
        if self.lapw:
            wb = 0.5 * self.r**2 * (2 * self.w - 1) * np.vdot(hbar, self.triangle).real
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
