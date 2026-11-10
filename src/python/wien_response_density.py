"""Read/write and transport WIEN radial-density/Fourier-star expansions.

Radial file coefficients contain r**2*n_lm, not n_lm.  For an affine cell
change lambda, c_lambda(r)=c_0(r/lambda)/lambda.  The reverse pullback of a
physical endpoint density is lambda*c_lambda(lambda*r).  Fourier density
star sums acquire lambda**-3 and lambda**3 respectively.  RMT stays fixed.
The tiny extrapolation across the MT edge must be checked by step halving;
this diagnostic does not assume a zero density at the sphere boundary.
"""

from pathlib import Path
import re
import copy
import numpy as np
from scipy.interpolate import CubicSpline

FLOAT = re.compile(r"[-+]?\d+\.\d+[ED][-+]\d+")


class Density:
    def __init__(self, path, structure):
        self.valence_format = Path(path).suffix.startswith(".clmval")
        self.lines = Path(path).read_text().splitlines()
        self.rad = {}
        self.spans = {}
        self.pw = {}
        self.pw_lines = {}
        self.grids = []
        self.mult = []
        for line in Path(structure).read_text().splitlines():
            m = re.search(r"NPT=\s*(\d+)\s+R0=\s*([\d.]+)\s+RMT=\s*([\d.]+)", line)
            if m:
                self.grids.append(
                    np.linspace(np.log(float(m[2])), np.log(float(m[3])), int(m[1]))
                )
            m = re.search(r"MULT=\s*(\d+)", line)
            if m:
                self.mult.append(int(m[1]))
        i = 0
        atom = -1
        while i < len(self.lines):
            line = self.lines[i]
            m = re.search(r"ATOMNUMBER\s*=?\s*(\d+)", line)
            if m:
                atom = int(m[1]) - 1
            m = re.search(r"CLM\(R\) FOR L\s*(-?\d+)\s+M=\s*(-?\d+)", line)
            if not m:
                i += 1
                continue
            key = (atom, int(m[1]), int(m[2]))
            i += 1
            while not self.lines[i].strip():
                i += 1
            start = i
            values = []
            while len(values) < len(self.grids[atom]):
                values.extend(
                    float(x.replace("D", "E")) for x in FLOAT.findall(self.lines[i])
                )
                i += 1
            if len(values) != len(self.grids[atom]):
                raise ValueError("bad radial length")
            self.rad[key] = np.array(values)
            # MIXER multiplies ONLY the valence l=m=0 channel by sqrt(4*pi)
            # before adding core and writing CLMSUM.  CLMCOR already has the
            # CLMSUM spherical convention.  Higher harmonics are unchanged.
            if self.valence_format and key[1:] == (0, 0):
                self.rad[key] *= np.sqrt(4 * np.pi)
            self.spans[key] = (start, i)
        for j, line in enumerate(self.lines):
            if "NUMBER OF PW" not in line:
                continue
            count = int(line.split()[0])
            for i in range(j + 1, j + 1 + count):
                line = self.lines[i]
                g = tuple(int(line[k : k + 5]) for k in [3, 8, 13])
                self.pw[g] = complex(float(line[18:37]), float(line[37:56]))
                self.pw_lines[g] = i
            break

    def copy(self):
        return copy.deepcopy(self)

    def zeros(self):
        out = self.copy()
        for v in out.rad.values():
            v[:] = 0
        for g in out.pw:
            out.pw[g] = 0j
        return out

    def add(self, other, scale=1):
        for k, v in other.rad.items():
            if k not in self.rad and np.max(abs(v)) > 1e-12:
                raise ValueError(f"missing angular channel {k}")
            if k in self.rad:
                self.rad[k] += scale * v
        for g, v in other.pw.items():
            if g not in self.pw and abs(v) > 1e-12:
                raise ValueError(f"missing Fourier star {g}")
            if g in self.pw:
                self.pw[g] += scale * v
        return self

    def scale(self, t):
        for k in self.rad:
            self.rad[k] *= t
        for k in self.pw:
            self.pw[k] *= t
        return self

    def transport(self, lam, pullback=False):
        out = self.copy()
        shift = np.log(lam) * (1 if pullback else -1)
        factor = lam if pullback else 1 / lam
        for k, v in self.rad.items():
            x = self.grids[k[0]]
            out.rad[k] = factor * CubicSpline(x, v)(x + shift)
        for g, v in self.pw.items():
            out.pw[g] = v * lam ** (3 if pullback else -3)
        return out

    def write(self, path):
        lines = self.lines.copy()
        for k, v in self.rad.items():
            a, b = self.spans[k]
            if self.valence_format and k[1:] == (0, 0):
                v = v / np.sqrt(4 * np.pi)
            rows = [
                "   " + "".join(f"{x:19.12E}" for x in v[i : i + 4])
                for i in range(0, len(v), 4)
            ]
            if len(rows) != b - a:
                raise ValueError("changed radial line layout")
            lines[a:b] = rows
        for g, z in self.pw.items():
            i = self.pw_lines[g]
            lines[i] = lines[i][:18] + f"{z.real:19.12E}{z.imag:19.12E}"
        Path(path).write_text("\n".join(lines) + "\n")

    def vector(self):
        return np.concatenate(
            [v for v in self.rad.values()]
            + [np.array(list(self.pw.values())).view(float)]
        )

    def from_vector(self, v):
        out = self.copy()
        i = 0
        for k in out.rad:
            n = len(out.rad[k])
            out.rad[k] = v[i : i + n].copy()
            i += n
        for g in out.pw:
            out.pw[g] = v[i] + 1j * v[i + 1]
            i += 2
        if i != len(v):
            raise ValueError("wrong density vector dimension")
        return out
