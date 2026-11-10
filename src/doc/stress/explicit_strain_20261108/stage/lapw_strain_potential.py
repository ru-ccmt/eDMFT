"""Native potential coordinates and explicit fixed-RMT strain ingredients.

The radial file convention is MIXED: spherical VSP stores r*V_00 (physical
spherical average), nonspherical VNS stores V_lm. The interstitial VNS entries
are masked per-member Fourier coefficients, in the file conjugation convention.
Changing geometry changes those coefficients even for a constant potential.

This module supplies strict file-coordinate matching, the analytic moving-mask
source and the hydrostatic MT coordinate-transport term. It does not claim to
generate the entire geometric Hartree/XC/ionic source or a total stress tensor.
Native centered potential samples remain a reference for the unfinished source.
"""

from dataclasses import dataclass
from pathlib import Path
import re
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn
from lapw_density_response import checked


FLOAT = re.compile(r"[-+]?\d+\.\d+[ED][-+]?\d+", re.I)
ATOM = re.compile(r"ATOMNUMBER\s*=\s*(\d+)")
LM = re.compile(r"VLM\(R\) FOR L\s*(-?\d+)\s+M=\s*(-?\d+)")


@dataclass(frozen=True)
class PotentialFile:
    radial: dict
    interstitial: dict


def read_potential_file(path, radial_lengths):
    """Read exactly NPT values per channel and label every Fourier coefficient.

    No positional band/star matching, regex spill into the next channel, silent
    absent harmonics, or extra complex conjugation. Both VSP and VNS use this
    parser; their physical units are assigned by NativePotentialFileLayout.
    """
    lines = Path(path).read_text().splitlines()
    radial, pw = {}, {}
    atom = None
    i = 0
    while i < len(lines):
        line = lines[i]
        am = ATOM.search(line)
        if am:
            atom = int(am[1]) - 1
            if not 0 <= atom < len(radial_lengths):
                raise ValueError("invalid potential atom index")
        lm = LM.search(line)
        if lm:
            if atom is None:
                raise ValueError("potential channel precedes atom header")
            key = (atom, int(lm[1]), int(lm[2]))
            if key in radial:
                raise ValueError("duplicate potential harmonic")
            values = []
            i += 1
            while i < len(lines) and len(values) < radial_lengths[atom]:
                if (
                    ATOM.search(lines[i])
                    or LM.search(lines[i])
                    or "NUMBER OF PW" in lines[i]
                ):
                    break
                values.extend(
                    float(x.replace("D", "E").replace("d", "e"))
                    for x in FLOAT.findall(lines[i])
                )
                i += 1
            a = np.array(values)
            if a.shape != (radial_lengths[atom],) or not np.all(np.isfinite(a)):
                raise ValueError("truncated, oversized or nonfinite potential harmonic")
            a.setflags(write=False)
            radial[key] = a
            continue
        if "NUMBER OF PW" in line:
            if pw:
                raise ValueError("multiple interstitial potential blocks")
            count = int(line.split()[0])
            if count < 1 or i + count >= len(lines):
                raise ValueError("truncated interstitial potential block")
            for row in lines[i + 1 : i + 1 + count]:
                g = tuple(int(row[j : j + 5]) for j in [3, 8, 13])
                z = complex(
                    float(row[18:37].replace("D", "E")),
                    float(row[37:56].replace("D", "E")),
                )
                if g in pw or not np.isfinite(z):
                    raise ValueError("duplicate/nonfinite potential star")
                pw[g] = z
            i += count
        i += 1
    if not radial:
        raise ValueError("no radial potential data")
    return PotentialFile(radial, pw)


class NativePotentialFileLayout:
    """Central, fixed-label selector between prepared fields and VSP/VNS files.

    Layout comes from actual central files. Missing requested labels at an
    endpoint fail. Additional endpoint stars are permitted only if exactly
    negligible: a changing active cutoff invalidates this derivative reference.
    Fixed RMT, radial grids, symmetry/axes and fractional positions must also
    be checked from structures by the caller; potential files omit those data.
    """

    def __init__(self, hartree, vsp, vns, *, complex_interstitial):
        self.lengths = tuple(g.n for g in hartree.integrations)
        self.keys = tuple(hartree.keys)
        self.harmonics = {
            a: tuple(key[1:] for key in vns.radial if key[0] == a)
            for a in range(len(self.lengths))
        }
        self.harmonics = {a: lm for a, lm in self.harmonics.items() if lm}
        self.representatives = tuple(vns.interstitial)
        self.complex = complex_interstitial
        self.selector, self.coordinates = hartree.native_potential_output(
            tuple(range(len(self.lengths))),
            self.harmonics,
            np.asarray(self.representatives, int),
            complex_interstitial=complex_interstitial,
        )
        self.central = self.pack(vsp, vns)
        self.central.setflags(write=False)

    def pack(self, vsp, vns):
        wanted = {(a, 0, 0) for a in range(len(self.lengths))}
        if set(vsp.radial) != wanted:
            raise ValueError("spherical potential atom/harmonic labels differ")
        expected = {(a, *lm) for a, ls in self.harmonics.items() for lm in ls}
        if set(vns.radial) != expected:
            raise ValueError("nonspherical potential labels differ")
        missing = set(self.representatives) - set(vns.interstitial)
        extra = set(vns.interstitial) - set(self.representatives)
        if missing or any(abs(vns.interstitial[g]) > 1e-12 for g in extra):
            raise ValueError("active interstitial potential cutoff differs")
        sphere = np.concatenate([vsp.radial[a, 0, 0] for a in range(len(self.lengths))])
        ns = {
            a: np.array([vns.radial[a, *lm] for lm in ls]).T
            for a, ls in self.harmonics.items()
        }
        ir = np.array([vns.interstitial[g] for g in self.representatives])
        if not self.complex:
            if np.max(abs(ir.imag)) > 1e-12:
                raise ValueError("complex coefficients in real potential layout")
            ir = ir.real
        return self.coordinates.pack(sphere, ns, ir)

    def read(self, vsp_path, vns_path):
        return self.pack(
            read_potential_file(vsp_path, self.lengths),
            read_potential_file(vns_path, self.lengths),
        )

    def covector(self, prepared_field_gradient):
        """Convert a full-field gradient only if its entire support is retained."""
        g = checked(
            prepared_field_gradient,
            (self.selector.input_size,),
            "field gradient",
            real=True,
        )
        result = self.selector.matvec(g)
        residual = g - self.selector.rmatvec(result)
        if np.linalg.norm(residual) > 1e-12 * max(1, np.linalg.norm(g)):
            raise ValueError(
                "functional gradient uses channels absent from potential files"
            )
        return result


def unequal_centered_derivative(minus, center, plus, minus_step, plus_step):
    """Three-point derivative at zero for the ACTUAL signed structure steps."""
    hm, hp = float(minus_step), float(plus_step)
    if not np.isfinite(hm + hp) or not hm < 0 < hp:
        raise ValueError("steps must straddle the central structure")
    center = np.asarray(center)
    if np.iscomplexobj(center) or not np.all(np.isfinite(center)):
        raise ValueError("finite real central potential required")
    minus = checked(minus, center.shape, "minus potential", real=True)
    plus = checked(plus, center.shape, "plus potential", real=True)
    # Difference form avoids summing three large central potentials.
    return (-hp / (hm * (hm - hp))) * (minus - center) + (-hm / (hp * (hp - hm))) * (
        plus - center
    )


def hydrostatic_mt_transport(radial, log_grids):
    """Convert fixed-r potential derivative to material-r derivative.

    Add this term ONLY if changing convention: for F=r*V_00, material F is
    F_lambda(lambda*r)/lambda, giving dF/dlog(r)-F; for nonspherical V_lm,
    material V is V_lambda(lambda*r), giving dV_lm/dlog(r). No density-like
    lambda**3 factor belongs to a potential. The fixed-r B operator must NOT
    receive this addition while its band-geometry source remains unchanged.
    """
    result = {}
    for key, v in radial.items():
        x = log_grids[key[0]]
        v = checked(v, (len(x),), "radial potential", real=True)
        result[key] = CubicSpline(x, v)(x, 1) - (v if key[1:] == (0, 0) else 0)
    return result


def strain_mask_derivative(reanalysis, reciprocal, positions, radii, volume, strain):
    """Analytic derivative of the native fixed-RMT interstitial mask.

    A_eta=(I+eta*strain)A; fractional sites and integer G labels are fixed.
    dG=-strain^T G, dOmega=tr(strain)*Omega. On the FULL native FFT box,
      dM_G=sum_a (4*pi*R_a**3/Omega)*phase_a*
              [tr(strain)*f(x)+direction*x*f'(x)],
      f(x)=j1(x)/x, x=|G|R, direction=G^T strain G/|G|^2.
    For hydrostatic strain=I the bracket reduces to j0(x). Finite RMT surface
    terms survive even for a constant potential. The native G=0 gauge is not
    altered. This term alone is not the full interstitial potential source.
    """
    b = checked(reciprocal, (3, 3), "reciprocal basis", real=True)
    e = checked(strain, (3, 3), "strain direction", real=True)
    positions = np.asarray(positions)
    checked(positions, (len(positions), 3), "fractional sites", real=True)
    radii = checked(radii, (len(positions),), "fixed radii", real=True)
    if (
        not np.allclose(e, e.T, rtol=0, atol=1e-14)
        or np.any(radii <= 0)
        or not np.isfinite(volume)
        or volume <= 0
    ):
        raise ValueError("symmetric strain and positive fixed-RMT geometry required")
    axes = [np.fft.fftfreq(n) * n for n in reanalysis.shape]
    labels = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)
    g = labels @ b.T
    norm = np.linalg.norm(g, axis=-1)
    direction = np.divide(
        np.einsum("...i,ij,...j->...", g, e, g),
        norm**2,
        out=np.zeros_like(norm),
        where=norm != 0,
    )
    result = np.zeros(reanalysis.shape, complex)
    for p, r in zip(positions, radii):
        x = norm * r
        f = np.divide(spherical_jn(1, x), x, out=np.full_like(x, 1 / 3), where=x != 0)
        xfprime = spherical_jn(0, x) - 3 * f
        phase = np.exp(2j * np.pi * (labels @ p))
        if reanalysis.stars.inversion:
            phase = phase.real
        result += (
            4 * np.pi * r**3 / volume * phase * (np.trace(e) * f + direction * xfprime)
        )
    return result


def masked_potential_strain_source(reanalysis, unmasked_star_sums, mask_derivative):
    """d(REAN4)/deta applied to FIXED unmasked potential star sums.

    Returns internal REAN4 per-member coefficients. The existing native file
    selector supplies conjugation once; do not conjugate before that selector.
    The other leg, REAN4*d(unmasked potential), must be added separately.
    """
    z = checked(unmasked_star_sums, (reanalysis.stars.nstar,), "unmasked potential")
    mask = checked(mask_derivative, reanalysis.shape, "mask derivative")
    grid = np.zeros(reanalysis.shape, complex)
    grid[reanalysis.indices] = z[reanalysis.stars.owners] * reanalysis.weight
    physical = np.fft.fftn(grid).real
    transformed = np.fft.ifftn(physical * np.fft.fftn(mask.conj()).real)
    return reanalysis._collect(
        transformed[reanalysis.indices] * reanalysis.weight.conj()
    )
