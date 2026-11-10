from pathlib import Path
r=Path(__file__).parent/'stage';p=r/'lapw_response_prepare.py';s=p.read_text()
s=s.replace('from pathlib import Path','from pathlib import Path\nfrom contextlib import chdir\nfrom lapw_response_policy import detect_policy\nfrom dmft2_spin import run_channels, run_fermif',1)
s=s.replace('if set(config) - {"native_mpi", "lapw1_mpi"} != required:', 'if set(config) - {"native_mpi", "lapw1_mpi", "physical_policy", "dmft2dn_definition", "fermif"} != required:')
s=s.replace('    work.mkdir(parents=True, exist_ok=False)','''    policy = detect_policy(reference, case, override=config.get("physical_policy"))
    if policy.magnetic and "dmft2dn_definition" not in config:
        raise ValueError("two-spin preparation requires dmft2dn_definition")
    work.mkdir(parents=True, exist_ok=False)''')
s=s.replace('        "indmf2",','        "indmf2",\n        "indmfldn",\n        "indmf2dn",')
s=s.replace('["EF.dat", "projectorw.dat", "Eorb_imp.dat"]','["EF.dat", "projectorw.dat", "Eorb_imp.dat", "params.dat"]')
a=s.index('    for suffix in ["inso"');b=s.index('    files = ',a);s=s[:a]+s[b:]
s=s.replace('    for label in ["lapw1", "dmft2"]:', '    for label in ["lapw1", "dmft2"] + (["dmft2dn"] if policy.magnetic else []):')
a=s.index('    control = work / (case + ".indmf2")');b=s.index('    times = []',a);s=s[:a]+s[b:]
a=s.index('    run("dmft2_number"');b=s.index('    run("lcore",',a)
s=s[:a]+'''    definitions = [str(work/('dmft2'+channel+'.def')) for channel in policy.channels]
    fermif = Path(config.get("fermif", Path(config["dmft2"]).parent/"fermif")).resolve()
    def common_fermi(bands):
        with (work/"fermif.log").open("w") as stream:
            run_fermif(fermif, bands, stream)
    with chdir(work):
        run_channels(definitions, policy,
            lambda label, definition: run(label+"_number", config["dmft2"], [definition], marker="DMFT2"),
            common_fermi)
        def export(label, definition):
            suffix = "_dn" if label.endswith("dn") else ""
            run(label, config["dmft2"], [definition],
                {"EDMFT_RESPONSE_EXPORT": str(work/("response"+suffix))}, "DMFT2")
        run_channels(definitions, policy, export, common_fermi, export=True)
'''+s[b:]
s=s.replace('    if not (work / "response.meta").is_file() or not list(', '    if any(not (work / ("response"+("_dn" if ch else "")+".meta")).is_file() for ch in policy.channels) or not list(')
s=s.replace('        case=case,\n        parent_fingerprints=', '        case=case,\n        physical_policy=policy.record(),\n        parent_fingerprints=',1)
s=s.replace('DMFT2 first determines the central chemical potential, then its\nfixed-mu export pass records that same state. The analytic derivative later\nenforces fixed N.', 'The configured physical ensemble determines whether DMFT2/fermif first\nupdates the common chemical potential. Fixed-mu exports then record both spin\nbranches. Preparation metadata retains the original response ensemble.')
p.write_text(s)
