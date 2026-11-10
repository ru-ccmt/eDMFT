from pathlib import Path
p=Path(__file__).parent/'stage'
f=p/'lapw_xc_response.py';s=f.read_text();s=s.replace('A quadratic least-squares fit on 27 neighbours reproduces PWXAD4/5.','A quadratic fit on the native lattice-specific neighbours reproduces PWXAD4/5.').replace('27-point variable-coefficient stencil for each forward/reverse call.','variable-coefficient stencil for each forward/reverse call.')
start=s.index('class InterstitialXCResponse:')
a=s[:start];b=s[start:]
b=b.replace('quantity="potential"\n    ):','quantity="potential",\n        lattice_type="P",\n    ):',1)
b=b.replace('''        self.offsets = np.array(
            [(i, j, k) for j in [-1, 0, 1] for i in [-1, 0, 1] for k in [-1, 0, 1]]
        )''','''        # PWXAD4 uses the 120-degree planar metric for H cells. Its
        # neighbour set is independent of the CornUse switch in that case.
        # Cartesian derivatives still use AVEC/FFT-size, never fractional
        # coordinate lengths. Keep native loop ordering for reproducibility.
        extent = 2 if lattice_type.strip() == "H" else 1
        self.offsets = np.array([
            (i, j, k)
            for j in range(-extent, extent + 1)
            for i in range(-extent, extent + 1)
            for k in (-1, 0, 1)
            if lattice_type.strip() != "H" or i*i + j*j - i*j + k*k < 3
        ])''',1).replace('np.ones(27), x, y','np.ones(len(self.offsets)), x, y',1)
f.write_text(a+b)
for name,key in [('lapw_potential_core_state.py','s'),('lapw_functional_cancellation.py','state.structure')]:
 f=p/name;s=f.read_text();needle='NativeLattice.from_structure('+key+').cell,';assert needle in s;s=s.replace(needle,needle+'\n        lattice_type='+key+'.lattyp,');f.write_text(s)
# Replace tests of deliberately lifted restrictions with malformed-input tests.
f=p/'tests/test_lapw_charge_response.py';s=f.read_text().replace('s.alpha = 89.0\n        with self.assertRaises(NotImplementedError):','s.a = 0.0\n        with self.assertRaises(ValueError):');f.write_text(s)
f=p/'tests/test_lapw_hartree_periodic.py';s=f.read_text().replace('s.lattyp = "F"\n        with self.assertRaises(NotImplementedError):','s.lattyp = "INVALID"\n        with self.assertRaises(ValueError):');f.write_text(s)
f=p/'tests/test_lapw_geometry_response.py';s=f.read_text().replace('structure=SimpleNamespace(timat=[rot]),','structure=SimpleNamespace(timat=[rot], lattyp="P", a=4., b=4., c=4.,\n                                      alpha=90., beta=90., gamma=90.),');f.write_text(s)
