"""Native unmixed output-charge assembly and its exact discrete derivative.

The output of dmft2 and lcore is normalized before mixing. This module follows
MIXER's default CHARGE/SETN/NORMAL operations, not the mixing algorithm and not
uniform rescaling of all density coefficients. The input density must already
be in CLMSUM units and include core. Geometry and target charge are fixed in
ordinary actions; explicit geometry changes have a separate source operation.
"""

from __future__ import annotations
from lapw_lattice import NativeLattice
import numpy as np
from scipy.special import spherical_jn
from lapw_density_response import checked
from lapw_hartree_periodic import native_charge_weights
from stress_response_adjoint import PairedOperator


def native_mixer_inversion(structure):
    """MIXER's exact real-field test: -I with a translation below 1e-5.

    An empty optional in1c file is not a reliable selector. A displaced
    inversion centre is deliberately treated as complex, as native MIXER does.
    """
    rotations = np.asarray(structure.timat)
    translations = np.asarray(structure.tau)
    if (
        rotations.ndim != 3
        or rotations.shape[1:] != (3, 3)
        or translations.shape != (len(rotations), 3)
    ):
        raise ValueError("invalid symmetry metadata")
    if not np.all(np.isfinite(translations)):
        raise ValueError("nonfinite symmetry translation")
    return bool(
        np.any(
            np.all(rotations == -np.eye(3, dtype=int), axis=(1, 2))
            & np.all(abs(translations) < 1e-5, axis=1)
        )
    )


def native_charge_covectors(density, structure, *, inversion):
    """Return (MT,IR) covectors in Density.vector()'s real coordinates.

    MT: MULT * native CHARGE(c00); c00=4*pi*r^2*n_spherical.
    IR: Re sum_G conj(Theta_G)*rho_G, where rho_G is a STAR SUM
    and Theta_G is the integral over the interstitial of exp(-i G.r).
    For a real native calculation SETN discards the imaginary mask part.
    No extra star size, spin weight, or sqrt(4*pi) is allowed here.

    The default mixer CHARGE rule is the same extended Simpson rule used by
    LAPW0. Alternate SIMP/LDMNORM settings require a separate adapter.
    """
    s = structure
    if not isinstance(inversion, (bool, np.bool_)):
        raise ValueError("explicit inversion flag required")
    lattice = NativeLattice.from_structure(s)
    if (
        len(density.grids) != s.nat
        or len(density.mult) != s.nat
        or not np.array_equal(density.mult, s.mult)
    ):
        raise ValueError("density/structure atom multiplicities differ")
    if len(s.pos) != sum(s.mult):
        raise ValueError("missing equivalent-site positions")
    vector = density.vector()
    checked(vector, vector.shape, "central density", real=True)
    mt = np.zeros_like(vector)
    ir = mt.copy()
    offset = 0
    spherical = set()
    for (a, l, m), c in density.rad.items():
        if not 0 <= a < s.nat or len(c) != s.jrj[a] or len(c) != len(density.grids[a]):
            raise ValueError("invalid native radial density layout")
        if l == m == 0:
            if a in spherical:
                raise ValueError("repeated spherical channel")
            spherical.add(a)
            r = np.exp(density.grids[a])
            if not np.allclose([r[0], r[-1]], [s.r0[a], s.rmt[a]], rtol=1e-12, atol=0):
                raise ValueError("density and structure radial grids differ")
            mt[offset : offset + len(c)] = s.mult[a] * native_charge_weights(r)
        offset += len(c)
    if spherical != set(range(s.nat)):
        raise ValueError("missing spherical density")
    reps = np.asarray(list(density.pw))
    if (
        reps.ndim != 2
        or reps.shape[1] != 3
        or not len(reps)
        or not np.issubdtype(reps.dtype, np.integer)
        or np.any(reps[0])
    ):
        raise ValueError("native charge list must start with G=0")
    volume = lattice.volume
    g = np.linalg.norm(lattice.momenta(reps), axis=1)
    mask = np.zeros(len(reps), complex)
    mask[0] = volume
    position = 0
    for a in range(s.nat):
        radius = float(s.rmt[a])
        x = g * radius
        radial = np.divide(
            spherical_jn(1, x), x, out=np.full_like(x, 1 / 3), where=x != 0
        )
        phases = np.exp(
            -2j * np.pi * (np.asarray(s.pos[position : position + s.mult[a]]) @ reps.T)
        ).sum(axis=0)
        mask -= 4 * np.pi * radius**3 * radial * phases
        position += s.mult[a]
    if inversion:
        mask = mask.real.astype(complex)
    ir[offset:] = np.ascontiguousarray(mask).view(float)
    return mt, ir


class NativeChargeNormalization:
    """Derivative N of MIXER's default interstitial-only normalization.

    Let q=(c,p) be the unnormalized total output density, including core;
    M=m.q, I=i.q, and alpha=(Qtarget-M)/I. The nonlinear map is
        F(q)=(c,alpha*p).
    At the central UNNORMALIZED output q0,
        N x = D x - p0 * ell.x,
        D=diag(1_MT,alpha0_PW), ell=(m+alpha0*i)/I0,
    where p0 is q0 with its MT entries set to zero. Hence
        N^T y = D y - ell*(p0.y).
    In general alpha0 != 1. Normalizing q0 before preparing N would give the
    wrong derivative of the actual upstream charge map.

    enabled=False explicitly represents input NORM=NO. Target variation and
    charge-integration geometry variation are external sources, not guessed.
    """

    def __init__(self, density, structure, target_charge, *, inversion, enabled=True):
        if not isinstance(enabled, (bool, np.bool_)):
            raise ValueError("explicit normalization policy required")
        if not np.isfinite(target_charge):
            raise ValueError("nonfinite target charge")
        self.base = np.array(density.vector(), float, copy=True)
        self.mt, self.ir = native_charge_covectors(
            density, structure, inversion=inversion
        )
        self.nrad = sum(map(len, density.rad.values()))
        self.size = len(self.base)
        self.target = float(target_charge)
        self.enabled = bool(enabled)
        self.M = float(self.mt @ self.base)
        self.I = float(self.ir @ self.base)
        if self.enabled and self.I <= 0:
            raise ValueError("nonpositive interstitial charge")
        self.alpha = (self.target - self.M) / self.I if self.enabled else 1.0
        if self.enabled and self.alpha <= 0:
            raise ValueError("target charge leaves no positive interstitial charge")
        self.diagonal = np.ones(self.size)
        self.diagonal[self.nrad :] = self.alpha
        self.p0 = self.base.copy()
        self.p0[: self.nrad] = 0
        self.ell = (
            (self.mt + self.alpha * self.ir) / self.I
            if self.enabled
            else np.zeros(self.size)
        )
        self.central = self.endpoint(self.base)
        self.operator = PairedOperator(
            self.size,
            self.size,
            self.action,
            self.pullback,
            "native unmixed charge normalization",
        )
        for a in [
            self.base,
            self.mt,
            self.ir,
            self.diagonal,
            self.p0,
            self.ell,
            self.central,
        ]:
            a.setflags(write=False)

    def action(self, x):
        x = checked(x, (self.size,), "raw output-density perturbation", real=True)
        return self.diagonal * x - self.p0 * (self.ell @ x)

    def pullback(self, y):
        y = checked(y, (self.size,), "normalized output-density covector", real=True)
        return self.diagonal * y - self.ell * (self.p0 @ y)

    def endpoint(self, x):
        x = checked(x, (self.size,), "raw output density", real=True)
        y = x.copy()
        if self.enabled:
            interstitial = self.ir @ x
            if interstitial <= 0:
                raise ValueError("nonpositive endpoint interstitial charge")
            alpha = (self.target - self.mt @ x) / interstitial
            if alpha <= 0:
                raise ValueError("nonpositive endpoint normalization")
            y[self.nrad :] *= alpha
        return y

    def explicit_source(
        self, *, mt_geometry=0.0, ir_geometry=0.0, target_derivative=0.0
    ):
        """Changing fixed-RMT cell integration/target at fixed coefficients.

        mt_geometry = (partial_nu m).q0 and similarly for IR. These must
        use the same quadrature and Fourier convention as the ordinary map.
        """
        if not np.all(np.isfinite([mt_geometry, ir_geometry, target_derivative])):
            raise ValueError("invalid normalization source")
        if not self.enabled:
            return np.zeros(self.size)
        return (
            self.p0
            * (target_derivative - mt_geometry - self.alpha * ir_geometry)
            / self.I
        )

    def compose(self, raw_charge_map):
        """N [R E B+D+C] K; all upstream physical maps remain explicit."""
        if raw_charge_map.output_size != self.size:
            raise ValueError("output charge layout differs")
        return self.operator @ raw_charge_map
