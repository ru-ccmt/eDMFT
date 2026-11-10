from pathlib import Path
r=Path(__file__).parent/'stage';p=r/'lapw_total_stress_driver.py';s=p.read_text()
s=s.replace('from pathlib import Path','from pathlib import Path\nfrom contextlib import ExitStack',1)
s=s.replace('    save_geometry=False,\n):','    save_geometry=False,\n    ensemble="auto",\n    diagnostic_gap_cutoff=None,\n):',1)
a=s.index('    for suffix in ["inso"');b=s.index('    files = {',a)
s=s[:a]+'''    from lapw_response_policy import detect_policy
    from lapw_magnetic_response import combine_valence
    provenance = snapshot / "preparation.json"
    override = json.loads(provenance.read_text()).get("physical_policy", {}) if provenance.is_file() else {}
    if ensemble not in ("auto", "fixed-n", "fixed-mu"):
        raise ValueError("invalid electronic ensemble")
    if ensemble != "auto":
        override["fixed_number"] = ensemble == "fixed-n"
    policy = detect_policy(snapshot, case, override=override)
'''+s[b:]
a=s.index('    mixer = ');b=s.index('    import re',a);s=s[:a]+s[b:]
a=s.index('    with load_native_valence_state(');b=s.index('        t = time.perf_counter()',a)
s=s[:a]+'''    with ExitStack() as stack:
        branches = [stack.enter_context(load_native_valence_state(
            snapshot, case, complex_vectors=complex_vectors, energy_policy="fixed",
            matrix_quadrature="wien2k24", library=library, parallel=parallel,
            cache_mb=cache_mb, basis_cache_mb=basis_cache_mb, spin_channel=channel,
            diagnostic_gap_cutoff=diagnostic_gap_cutoff,
        )) for channel in policy.channels]
        v = combine_valence(*branches) if policy.magnetic else branches[0]
        timings["prepare_valence_seconds"] = v.preparation_seconds
        state = NativeScreenedState(
            v, pc, core, target_charge=target, fixed_number=policy.fixed_number,
            normalization_enabled=policy.normalization_enabled,
        )
'''+s[b:]
s=s.replace('            fixed_point_tolerance=fixed_point_tolerance,','            fixed_point_tolerance=fixed_point_tolerance,\n            diagnostic=diagnostic_gap_cutoff is not None,',1)
s=s.replace('            case=case,\n            snapshot=', '            case=case,\n            physical_policy=policy.record(),\n            snapshot=',1)
s=s.replace('electronic_cache=v.model.cache_info(),','electronic_cache=[b.model.cache_info() for b in branches],')
s=s.replace('    a = p.parse_args()','''    p.add_argument("--ensemble", choices=["auto", "fixed-n", "fixed-mu"], default="auto",
                   help="auto uses preparation provenance; override for legacy exports")
    p.add_argument("--diagnostic-gap-cutoff", type=float,
                   help="omit singular window pairs ONLY in a non-publishable diagnostic")
    a = p.parse_args()''')
s=s.replace('        save_geometry=a.save_geometry,','        save_geometry=a.save_geometry,\n        ensemble=a.ensemble,\n        diagnostic_gap_cutoff=a.diagnostic_gap_cutoff,')
p.write_text(s)
