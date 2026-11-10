from pathlib import Path
p=Path(__file__).parent/'stage'
f=p/'lapw_valence_state.py';s=f.read_text().replace('from lapw_lattice import NativeLattice','from lapw_lattice import NativeLattice, native_site_frames');a=s.index('    axes = []\n    locations = []');b=s.index('    geometry = []',a);s=s[:a]+'''    first = np.r_[0, np.cumsum(s.mult)[:-1]]
    site_atom = np.repeat(np.arange(s.nat), s.mult)
    axes, locations = native_site_frames(s)
'''+s[b:];f.write_text(s)
f=p/'lapw_geometry_response.py';s=f.read_text().replace('from lapw_lattice import NativeLattice','from lapw_lattice import NativeLattice, native_site_frames');a=s.index('    axes, locations = [], []');b=s.index('    lo, offset',a);s=s[:a]+'''    axes, locations = native_site_frames(s)
'''+s[b:];f.write_text(s)
