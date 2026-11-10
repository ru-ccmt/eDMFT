from pathlib import Path
import sys,json,numpy as np
r=Path(__file__).parent;sys.path.insert(0,str(r/'stage'))
from lapw_valence_state import load_native_valence_state
from lapw_geometry_response import native_band_geometry_source
w=r/'native/MnO/MnO_center'
with load_native_valence_state(w,'MnO',complex_vectors=False,energy_policy='fixed',matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so',diagnostic_gap_cutoff=1e-10,cache_mb=64) as v:
 b,report=native_band_geometry_source(v,w,'MnO',np.eye(3))
 (r/'MnO_hamiltonian_probe.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
