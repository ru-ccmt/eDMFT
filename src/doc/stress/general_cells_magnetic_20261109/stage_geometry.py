from pathlib import Path
root=Path(__file__).parent/'stage'
def change(name, pairs):
 p=root/name;s=p.read_text()
 for a,b in pairs:
  if a not in s: raise RuntimeError(name+' missing '+repr(a[:100]))
  s=s.replace(a,b)
 if 'NativeLattice' in s and 'from lapw_lattice import NativeLattice' not in s:
  # Insert after any future import, otherwise immediately before numpy.
  target='import numpy as np';assert target in s;s=s.replace(target,'from lapw_lattice import NativeLattice\n'+target,1)
 p.write_text(s)
change('lapw_charge_response.py',[( '''    if s.lattyp.strip() != "P" or not np.allclose(
        [s.alpha, s.beta, s.gamma], 90, rtol=0, atol=1e-10
    ):
        raise NotImplementedError(
            "native charge adapter requires primitive orthogonal cells"
        )
    cell = checked([s.a, s.b, s.c], (3,), "cell lengths", real=True)
    if np.any(cell <= 0):
        raise ValueError("nonpositive cell")''','    lattice = NativeLattice.from_structure(s)'),('volume = float(np.prod(cell))','volume = lattice.volume'),('g = np.linalg.norm(2 * np.pi * reps / cell, axis=1)','g = np.linalg.norm(lattice.momenta(reps), axis=1)')])
change('wien_response_operators.py',[( '''    if s.lattyp.strip() != "P" or not np.allclose([s.alpha, s.beta, s.gamma], 90.0):
        raise ValueError(
            "charge adapter currently requires a primitive orthogonal cell"
        )''','    lattice = NativeLattice.from_structure(s)'),('    cell = np.array([s.a, s.b, s.c])\n    volume = float(np.prod(cell))','    volume = lattice.volume'),('magnitude = np.linalg.norm(2 * np.pi * g / cell)','magnitude = np.linalg.norm(lattice.momenta(g))')])
change('lapw_native_potential_export.py',[('if not np.isclose(volume, s.a * s.b * s.c, rtol=1e-10):','if not np.isclose(volume, NativeLattice.from_structure(s).volume, rtol=1e-10):')])
change('lapw_hartree_periodic.py',[( '''        if s.lattyp != "P" or not np.allclose(
            [s.alpha, s.beta, s.gamma], 90, atol=1e-10
        ):
            raise NotImplementedError(
                "initial Hartree adapter requires primitive orthogonal cells"
            )''','        lattice = NativeLattice.from_structure(s)\n        self.lattice = lattice'),('reciprocal = np.diag(2 * np.pi / np.array([s.a, s.b, s.c]))','reciprocal = lattice.reciprocal'),('volume = s.a * s.b * s.c','volume = lattice.volume'),('def __init__(self, stars, shape, reciprocal, positions, radii, volume):','def __init__(self, stars, shape, reciprocal, positions, radii, volume, *, lattice=None):'),('        self.mask_fourier = mask','''        self.centering_mask = (lattice.centering_factor(labels) if lattice is not None
                               else np.ones(self.shape))
        mask *= self.centering_mask
        self.mask_fourier = mask'''),('self.stars, fft_shape, reciprocal, s.pos, site_radii, volume','self.stars, fft_shape, reciprocal, s.pos, site_radii, volume, lattice=lattice')])
change('lapw_potential_core_state.py',[('        [s.a, s.b, s.c],\n        functional=m["functional"],','        NativeLattice.from_structure(s).cell,\n        functional=m["functional"],')])
change('lapw_xc_response.py',[( '''        length = checked(cell_lengths, (3,), "orthogonal cell lengths", real=True)
        if np.any(length <= 0):
            raise ValueError("cell lengths must be positive")''','''        cell = np.asarray(cell_lengths)
        if cell.shape == (3,):
            if np.any(cell <= 0):
                raise ValueError("cell lengths must be positive")
            cell = np.diag(cell)
        cell = checked(cell, (3, 3), "Cartesian FFT cell", real=True)
        if abs(np.linalg.det(cell)) < 1e-14:
            raise ValueError("singular FFT cell")'''),('xyz = self.offsets * (length / np.array(self.shape))','xyz = (self.offsets / np.array(self.shape)) @ cell.T')])
change('lapw_hydrostatic_potential.py',[('b = np.diag(2 * np.pi / np.array([s.a, s.b, s.c]))','b = NativeLattice.from_structure(s).reciprocal')])
change('lapw_strain_field.py',[('for rot in s.timat):','for rot in NativeLattice.from_structure(s).cartesian_rotations(s.timat)):'),('b = np.diag(2 * np.pi / np.array([s.a, s.b, s.c]))','b = NativeLattice.from_structure(s).reciprocal')])
# The exact mask derivative must use the SAME centering projector as REAN4.
p=root/'lapw_strain_potential.py';s=p.read_text();start=s.index('def strain_mask_derivative(');tail=s[start:];print(tail[tail.index('    result ='):tail.index('    result =')+1100]);p.write_text(s)
