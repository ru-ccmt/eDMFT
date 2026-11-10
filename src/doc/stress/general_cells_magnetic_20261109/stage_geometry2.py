from pathlib import Path
root=Path(__file__).parent/'stage'
def change(name,pairs):
 p=root/name;s=p.read_text()
 for a,b in pairs:
  assert a in s,(name,a[:100]);s=s.replace(a,b)
 if 'NativeLattice' in s and 'from lapw_lattice import NativeLattice' not in s:s=s.replace('import numpy as np','from lapw_lattice import NativeLattice\nimport numpy as np',1)
 p.write_text(s)
change('lapw_strain_potential.py',[('''            4 * np.pi * r**3 / volume * phase * (np.trace(e) * f + direction * xfprime)
        )
    return result''','''            4 * np.pi * r**3 / volume * phase * (np.trace(e) * f + direction * xfprime)
        )
    return result * getattr(reanalysis, "centering_mask", 1.0)''')])
change('lapw_functional_cancellation.py',[('[state.structure.a, state.structure.b, state.structure.c],','NativeLattice.from_structure(state.structure).cell,')])
change('lapw_strain_charge.py',[( '''    if s.lattyp.strip() != "P" or not np.allclose(
        [s.alpha, s.beta, s.gamma], 90, atol=1e-10, rtol=0
    ):
        raise NotImplementedError("primitive orthogonal charge geometry required")''','    lattice = NativeLattice.from_structure(s)'),('    cell = np.array([s.a, s.b, s.c])\n    q = 2 * np.pi * reps / cell','    q = lattice.momenta(reps)'),('np.trace(e) * np.prod(cell)','np.trace(e) * lattice.volume'),('for rot in s.timat):','for rot in NativeLattice.from_structure(s).cartesian_rotations(s.timat)):'),('s.a * s.b * s.c,','NativeLattice.from_structure(s).volume,')])
change('lapw_geometry_response.py',[('''            delta -= np.rint(delta)
            if max(abs(delta)) < 2e-7:''','''            if NativeLattice.from_structure(s).equivalent_displacement(delta):'''),('axes.append(rot.T)','axes.append(NativeLattice.from_structure(s).cartesian_rotations(rot).T)'),('momenta = [2 * np.pi * g / np.array([s.a, s.b, s.c]) for g in fractional]','momenta = [NativeLattice.from_structure(s).momenta(g) for g in fractional]'),('for rot in s.timat:\n        if np.linalg.norm(e @ rot - rot @ e)','for rot in NativeLattice.from_structure(s).cartesian_rotations(s.timat):\n        if np.linalg.norm(e @ rot - rot @ e)'),('volume = s.a * s.b * s.c','volume = NativeLattice.from_structure(s).volume')])
change('lapw_valence_state.py',[( '''    if s.lattyp.strip() != "P" or not np.allclose(
        [s.alpha, s.beta, s.gamma], 90, atol=1e-9, rtol=0
    ):
        raise NotImplementedError("only primitive orthogonal cells are supported")''','    NativeLattice.from_structure(s)  # Validate both stored and primitive cell maps.'),('''            diff -= np.rint(diff)
            if max(abs(diff)) < 2e-7:''','''            if NativeLattice.from_structure(s).equivalent_displacement(diff):'''),('axes.append(rotation.T)','axes.append(NativeLattice.from_structure(s).cartesian_rotations(rotation).T)'),('q = 2 * np.pi * frac / np.array([s.a, s.b, s.c])','q = NativeLattice.from_structure(s).momenta(frac)'),('volume = s.a * s.b * s.c','volume = NativeLattice.from_structure(s).volume')])
# No publication/preparation guards are removed in this foundation patch.
# End-to-end magnetic/fixed-mu and cubic-harmonic validation remain mandatory.
