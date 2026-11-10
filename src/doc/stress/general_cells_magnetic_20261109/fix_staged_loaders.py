from pathlib import Path
root=Path(__file__).parent;p=root/'stage'
f=root/'probe_native_density.py';s=f.read_text().replace("['struct','in0','clmsum']","['struct','in0','inm','clmsum']");f.write_text(s)
f=p/'lapw_valence_state.py';s=f.read_text().replace('''    result.fingerprints = fingerprints
''','''    result.spin_channel = spin_channel
    result.fingerprints = fingerprints
''',1);f.write_text(s)
f=p/'lapw_geometry_response.py';s=f.read_text().replace('''    def validate():
        for name in [case + ".indmfl", case + ".vns", case + ".struct"]:''','''    channel = getattr(state, "spin_channel", "")
    def validate():
        for name in [case + ".indmfl" + channel, case + ".vns", case + ".struct"]:''').replace('''inl = Indmfl(str(directory / case))''','''inl = Indmfl(str(directory / case), extn="indmfl" + channel)''');f.write_text(s)
