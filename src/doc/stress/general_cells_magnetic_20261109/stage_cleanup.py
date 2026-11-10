from pathlib import Path
p=Path(__file__).parent/'stage'
f=p/'lapw_valence_state.py';s=f.read_text();a=s.index('        import os\n        audit = os.environ.get("EDMFT_RESPONSE_RAW_AUDIT")');b=s.index('        difference=',a);s=s[:a]+s[b:];s=s.replace('''    fingerprints = parallel.fingerprints(paths)
''','''    if (directory / "EF.dat").is_file():
        paths["EF.dat"] = directory / "EF.dat"
    fingerprints = parallel.fingerprints(paths)
''',1);f.write_text(s)
f=p/'lapw_magnetic_response.py';s=f.read_text();needle='''    complement=CollinearComplement(up.complement,down.complement)
''';s=s.replace(needle,'''    # Native DMFT uses this precise Ry/eV conversion. Both exports must refer
    # to the common EF.dat, rather than silently retaining independent mus.
    ef_up = getattr(up, "fingerprint_paths", {}).get("EF.dat")
    ef_down = getattr(down, "fingerprint_paths", {}).get("EF.dat")
    if ef_up is None or ef_down is None or up.fingerprints["EF.dat"] != down.fingerprints["EF.dat"]:
        raise ValueError("magnetic response requires the same recorded EF.dat")
    mu = float(ef_up.read_text().split()[0]) / 13.60569253
    if not np.isfinite(mu) or not np.isclose(mu, up.model.meta.mu, rtol=0, atol=1e-10):
        raise ValueError("exported common chemical potential differs from EF.dat")
'''+needle);f.write_text(s)
f=p/'lapw_screened_state.py';s=f.read_text();needle='''        self.metadata = dict(
''';s=s.replace(needle,'''        branches = getattr(valence, "branches", [valence])
        self.metadata = dict(
            dmft_spin_channels=len(branches),
            common_mu_Ry=float(valence.model.meta.mu),
            diagnostic_only=any(getattr(v.model, "diagnostic_only", False) for v in branches),
''',1);f.write_text(s)
changes={
'lapw_hartree_periodic.py':('''The initial adapter deliberately rejects cubic-harmonic compression and
non-primitive/non-orthogonal cells. It uses the default CoulEM/CHARGE rules;''','''Stored G coordinates and primitive volume follow WIEN LATGEN. Cubic site
harmonics use the native POTFAC boundary and ENERGY pairing conventions.
The adapter uses the default CoulEM/CHARGE rules;'''),
'lapw_hydrostatic_potential.py':('Only hydrostatic primitive orthogonal cells are supplied.','Cell geometry follows WIEN LATGEN, including centered/nonorthogonal cells.'),
'lapw_geometry_response.py':('''The adapter supports the prepared nonmagnetic primitive orthogonal snapshots
and symmetry-preserving strains.''','''The adapter supports prepared WIEN cell geometries and symmetry-preserving
strains, with separate electronic branches for DMFT-only magnetism.'''),
'lapw_density_functional.py':('Nonmagnetic LDA/PBE, primitive orthogonal cells, fixed RMT, neutral native','Unpolarized LDA/PBE potential, WIEN cell geometry, fixed RMT, neutral native'),
'lapw_functional_cancellation.py':('''real-coordinate pairing is n^T W v, with W diagonal:''','''real-coordinate pairing is n^T W v. For noncubic site harmonics W is diagonal:'''),
'lapw_total_stress.py':('''Supported: the primitive orthogonal, nonmagnetic/no-SOC, fixed-N projector-5
scope of the native loaders''','''Supported publication: symmetry-enforced diagonal, nonmagnetic/no-SOC,
fixed-N projector-5 scope of the native loaders'''),
}
for name,(old,new) in changes.items():
 f=p/name;s=f.read_text();assert old in s,name;f.write_text(s.replace(old,new))
f=p/'lapw_functional_cancellation.py';s=f.read_text().replace('''    Cubic compressed-harmonic conventions are rejected by the layout builder.''','''    Cubic site multiplets carry the additional C^T C harmonic metric from
    native ENERGY; their radial weights are shared within each multiplet.''');f.write_text(s)
