from pathlib import Path
p=Path(__file__).parent/'stage'
f=p/'lapw_total_stress.py';s=f.read_text();a=s.index('    if s.lattyp.strip() != "P"');b=s.index('    if not rotations',a);s=s[:a]+'''    from lapw_lattice import NativeLattice
    rotations = list(NativeLattice.from_structure(s).cartesian_rotations(s.timat))
'''+s[b:];f.write_text(s)
f=p/'lapw_strain_field.py';s=f.read_text().replace('''    if not np.allclose(e, np.diag(np.diag(e)), atol=1e-14, rtol=0):
        raise NotImplementedError("only diagonal cell strains are supplied")''','''    if not np.allclose(e, e.T, atol=1e-14, rtol=0):
        raise ValueError("symmetric cell strain required")''');a=s.index('''        diag = np.diag(e)''');b=s.index('''        point = ''',a);s=s[:a]+'''        direction[..., 1:4] = -fields[..., 1:4] @ e
        pairs = [(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]
        hessian = np.empty(fields.shape[:-1] + (3,3))
        for i,(a,b) in enumerate(pairs):
            hessian[...,a,b] = hessian[...,b,a] = fields[...,4+i]
        derivative = -e @ hessian - hessian @ e
        for i,(a,b) in enumerate(pairs):
            direction[...,4+i] = derivative[...,a,b]
'''+s[b:];f.write_text(s)
