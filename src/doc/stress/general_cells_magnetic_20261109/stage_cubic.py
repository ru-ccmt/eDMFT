from pathlib import Path
p=Path(__file__).parent/'stage'
(p/'lapw_cubic_harmonics.py').write_text('''"""WIEN cubic harmonic coordinate maps (SUMFAC, SROLYL, CUB_XC_BACK).

Stored density channels remain separate real harmonics. Native XC first
contracts each cubic multiplet with c, evaluates the combined harmonic,
and stores the result in its FIRST channel divided by c_first. ENERGY
pairs both fields after contraction: the metric on each group is c c^T.
This metric is not the identity and the output adapter is not C^T.
"""
import numpy as np

# Normalized cubic-harmonic coefficients, using WIEN signed-l labels.
_SQUARED = {
 (0,0):(1,1.), (3,2):(1,1.),
 (4,0):(1,7/12), (4,4):(1,5/12),
 (6,0):(1,1/8), (6,2):(1,11/16), (6,4):(-1,7/8), (6,6):(-1,5/16),
 (7,2):(1,13/24), (7,6):(1,11/24),
 (8,0):(1,33/64), (8,4):(1,7/48), (8,8):(1,65/192),
 (9,2):(1,3/16), (9,4):(1,17/24), (9,6):(-1,13/16), (9,8):(-1,7/24),
 (10,0):(1,65/384), (10,2):(1,247/384), (10,4):(-1,11/32),
 (10,6):(1,19/768), (10,8):(-1,187/384), (10,10):(-1,85/256),
}

def coefficient(l,m):
    try:
        sign, square = _SQUARED[abs(l),m]
    except KeyError as exc:
        raise ValueError("unsupported cubic harmonic") from exc
    return sign*np.sqrt(square)


def cubic_maps(harmonics):
    """Return C (stored -> compressed) and D (compressed -> stored XC)."""
    lm = list(map(tuple, harmonics))
    groups=[]; i=0
    while i<len(lm):
        l,m=lm[i]
        count=1 if (l,m) in [(0,0),(-3,2)] else 2 if l in (4,6,-7,-9) else 3 if l in (8,10) else 0
        if not count or lm[i:i+count] != [(l,m+4*j) for j in range(count)]:
            raise ValueError("invalid or incomplete native cubic multiplet")
        groups.append((i,count)); i+=count
    c=np.zeros((len(groups),len(lm))); d=np.zeros(c.T.shape)
    for row,(i,count) in enumerate(groups):
        for j in range(count): c[row,i+j]=coefficient(*lm[i+j])
        d[i,row]=1/c[row,i]
    return c,d


def boundary_factor(l,m):
    """POTFAC's cubic boundary conversion; radial volume term is unchanged."""
    return 1. if (l,m) in [(0,0),(-3,2)] else .5/coefficient(l,m)**2
''')
f=p/'lapw_hartree_periodic.py';s=f.read_text();s=s.replace('''        if any(i >= 0 for i in s.iatnr):
            raise NotImplementedError(
                "compressed cubic harmonics are not yet supported"
            )
''','''        from lapw_cubic_harmonics import cubic_maps, boundary_factor
        self.cubic_pairing = []
        self.boundary_factors = np.array([
            boundary_factor(l, m) if s.iatnr[a] > 0 else 1.
            for a, l, m in template.rad
        ])
''');needle='''        self.pw_offset = offset
''';s=s.replace(needle,'''        for a in range(s.nat):
            if s.iatnr[a] > 0:
                indices = [i for i, key in enumerate(self.keys) if key[0] == a]
                c, _ = cubic_maps([self.keys[i][1:] for i in indices])
                self.cubic_pairing.append(([self.slices[i] for i in indices], c.T @ c))
'''+needle,1);s=s.replace('''                    factor * spherical_jn(l, norm * s.rmt[atom]) * y[lm_index[l, m]]''','''                    self.boundary_factors[k] * factor
                    * spherical_jn(l, norm * s.rmt[atom]) * y[lm_index[l, m]]''');f.write_text(s)
for name,needle in [('lapw_hydrostatic_potential.py','''                factor
                * (-x'''),('lapw_strain_field.py','''            dB[k, 1:] = factor * (''')]:
 f=p/name;s=f.read_text();assert needle in s
 repl=needle.replace('factor','h.boundary_factors[k] * factor');s=s.replace(needle,repl);f.write_text(s)
f=p/'lapw_xc_response.py';s=f.read_text();start=s.index('class MuffinTinXCResponse:');end=s.index('def cartesian_invariants',start);a=s[:start];b=s[start:end];z=s[end:];b=b.replace('quantity="potential"\n    ):','quantity="potential",\n        cubic=False,\n    ):',1);b=b.replace('''        self.Y = np.array(angular, float)''','''        self.cubic = bool(cubic)
        self.native_angular = np.array(angular, float)
        self.native_angular_derivatives = angular_derivatives
        self.Y = self.native_angular.copy()''',1);b=b.replace('''        yw = self.Y * self.weights[:, None]''','''        projection_y = self.Y
        if self.cubic:
            from lapw_cubic_harmonics import cubic_maps
            if self.harmonics is None:
                raise ValueError("cubic XC requires signed harmonic labels")
            c, d = cubic_maps(self.harmonics)
            # Exported YY has the compressed entries at its front; the tail
            # is workspace and must not be interpreted as ordinary harmonics.
            self.Y = self.native_angular[:, :len(c)] @ c
            projection_y = self.native_angular[:, :len(c)] @ d.T
        yw = projection_y * self.weights[:, None]''',1);b=b.replace('''                np.asarray(a, float) for a in angular_derivatives''','''                (np.asarray(a, float)[:, :len(c)] @ c if self.cubic else np.asarray(a, float))
                for a in angular_derivatives''',1);f.write_text(a+b+z)
f=p/'lapw_potential_core_state.py';s=f.read_text().replace('''                harmonics=lm,
''','''                harmonics=lm,
                cubic=s.iatnr[a] > 0,
''');f.write_text(s)
f=p/'lapw_functional_cancellation.py';s=f.read_text().replace('''        self.weights.setflags(write=False)
''','''        self.cubic_pairing = getattr(layout, "cubic_pairing", [])
        self.weights.setflags(write=False)
''',1);s=s.replace('''        return self.weights * checked(vector, (self.size,), "pairing input", real=True)''','''        result = self.weights * checked(vector, (self.size,), "pairing input", real=True)
        # Native ENERGY contracts each cubic group before pairing. The
        # symmetric C^T C metric is essential for both derivative directions.
        for slices, metric in self.cubic_pairing:
            values = np.column_stack([result[sl] for sl in slices]) @ metric
            for j, sl in enumerate(slices):
                result[sl] = values[:, j]
        return result''',1).replace('''angular_derivatives=source.angular[1:],''','''angular_derivatives=source.native_angular_derivatives,''').replace('''                    source.Y,
''','''                    source.native_angular,
''').replace('''                    harmonics=source.harmonics,
''','''                    harmonics=source.harmonics,
                    cubic=source.cubic,
''');f.write_text(s)
# Malformed cubic multiplets remain a fail-closed input check.
f=p/'tests/test_lapw_hartree_periodic.py';s=f.read_text().replace('s.iatnr = [1]\n        with self.assertRaises(NotImplementedError):','s.iatnr = [1]\n        with self.assertRaises(ValueError):');f.write_text(s)
