from pathlib import Path
p=Path(__file__).parent/'stage'
f=p/'lapw_magnetic_response.py';s=f.read_text();a=s.index('    if up.bands.input!=down.bands.input:');b=s.index('    complement=CollinearComplement',a);s=s[:a]+'''    for name in ('size','spherical_size','nonspherical_shapes','nstar','complex_interstitial'):
        if getattr(up.bands.input,name)!=getattr(down.bands.input,name):
            raise ValueError("spin potential layout differs")
'''+s[b:];f.write_text(s)
f=p/'lapw_complement_response.py';s=f.read_text();needle='''    nk, nr = model.model.eps.shape
''';assert s.count(needle)==1;s=s.replace(needle,'''    if hasattr(model, "electronic_operator"):
        return model.electronic_operator(fixed_number=fixed_number)
'''+needle);f.write_text(s)
f=p/'lapw_strain_functional.py';s=f.read_text();needle='''    projector, coordinates = projector_functional_gradient(complement)
''';s=s.replace(needle,'''    if hasattr(complement, "functional_gradients"):
        return complement.functional_gradients()
'''+needle);f.write_text(s)
f=p/'lapw_valence_state.py';s=f.read_text();s=s.replace('def _snapshot_inputs(case, source, complex_vectors):','def _snapshot_inputs(case, source, complex_vectors, spin_channel=""):');s=s.replace('(case + ".indmfl")','(case + ".indmfl" + spin_channel)')
# Both public and internal loader receive explicit channel/cutoff options.
s=s.replace('''    basis_cache_mb=128.0,
):''','''    basis_cache_mb=128.0,
    spin_channel="",
    diagnostic_gap_cutoff=None,
):''');s=s.replace('''    if matrix_quadrature not in ("atparn", "wien2k24"):''','''    if spin_channel not in ("", "dn"):
        raise ValueError("spin_channel must be empty or dn")
    response_name = "response" + ("_dn" if spin_channel else "")
    if matrix_quadrature not in ("atparn", "wien2k24"):''',1)
s=s.replace('''            ".indmfl",
''','''            ".indmfl" + spin_channel,
''',1).replace('''            ".clmval",
''','''            ".clmval" + spin_channel,
''',1).replace('''["projectorw.dat", "response.meta"]''','''["projectorw.dat", response_name + ".meta"]''',1).replace('''        "response.k????????.bin",
''','''        response_name + ".k????????.bin",
''',1).replace('''            basis_cache_mb=basis_cache_mb,
''','''            basis_cache_mb=basis_cache_mb,
            spin_channel=spin_channel,
            diagnostic_gap_cutoff=diagnostic_gap_cutoff,
''',1)
s=s.replace('''    mat = _snapshot_inputs(case, central, complex_vectors)
''','''    mat = _snapshot_inputs(case, central, complex_vectors, spin_channel)
''',1).replace('''        central / "response", parallel=parallel, cache_mb=cache_mb
''','''        central / ("response_dn" if spin_channel else "response"),
        parallel=parallel, cache_mb=cache_mb, diagnostic_gap_cutoff=diagnostic_gap_cutoff,
''',1).replace('''        cache_mb=cache_mb,
    )
    del exported''','''        cache_mb=cache_mb,
        diagnostic_gap_cutoff=diagnostic_gap_cutoff,
    )
    del exported''',1).replace('''(case + ".clmval")''','''(case + ".clmval" + spin_channel)''').replace('''inl = Indmfl(str(central / case))''','''inl = Indmfl(str(central / case), extn="indmfl" + spin_channel)''');f.write_text(s)
# Geometry is spin independent, but its band projection is not. Keep each
# branch's KS/orbital gauge; combine only in common spatial coordinates.
f=p/'lapw_geometry_response.py';s=f.read_text();start=s.index('def native_band_geometry_source(');i=s.index('    s = ',start);s=s[:i]+'''    if hasattr(state, "branches"):
        pieces = [native_band_geometry_source(v, directory, case, strain) for v in state.branches]
        return state.bands.output.pack(*(x[0] for x in pieces)), dict(spin_branches=[x[1] for x in pieces])
'''+s[i:];# confirm actual argument below before saving
f.write_text(s)
f=p/'lapw_strain_charge.py';s=f.read_text();needle='''    s = valence.structure
''';assert s.count(needle)==1;s=s.replace(needle,'''    if hasattr(valence, "branches"):
        return .5 * sum(native_spatial_geometry_source(v, directory, case, strain) for v in valence.branches)
'''+needle);f.write_text(s)
f=p/'lapw_energy_response.py';s=f.read_text();needle='''    op = valence.bands.spherical
''';s=s.replace(needle,'''    if hasattr(valence, "branches"):
        pieces = [energy_direction(v, energies) for v in valence.branches]
        return valence.bands.output.pack(*(x[0] for x in pieces)), .5 * sum(x[1] for x in pieces)
'''+needle,1);f.write_text(s)
