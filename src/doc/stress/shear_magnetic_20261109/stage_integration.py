from pathlib import Path
r=Path(__file__).parent/'stage';p=r/'dmft2_step.py';s=p.read_text()
s=s.replace('    "backend",','    "backend",\n    "physical_policy",',1)
s=s.replace('    ".indmf2",','    ".indmf2",\n    ".indmfldn",\n    ".indmf2dn",',1)
a=s.index('    # Backend initially');b=s.index('    control = ',a)
s=s[:a]+'''    for suffix in [".inso", ".clmup", ".clmdn", ".incup", ".incdn"]:
        p = Path(case + suffix)
        if p.is_file() and p.stat().st_size:
            raise ValueError("response requires no SOC and an unpolarized KS/core potential")
'''+s[b:]
s=s.replace('int(control[2]) != 1','int(control[2]) not in (0, 1)',1).replace('response currently requires charge mode c and recomputeEF=1 (fixed N)','response requires charge mode c and recomputeEF=0 or 1',1)
a=s.index('    if s.lattyp.strip() != "P"');b=s.index('    header = ',a)
s=s[:a]+'''    from lapw_lattice import NativeLattice
    from lapw_response_policy import detect_policy
    NativeLattice.from_structure(s)
    policy = detect_policy(".", case, control=defs[2], override=solver.get("physical_policy"))
    if backend == "native" and (policy.magnetic or not policy.fixed_number or s.lattyp.strip() != "P" or not np.allclose([s.alpha,s.beta,s.gamma],90)):
        raise ValueError("legacy native backend requires nonmagnetic fixed-N primitive orthogonal inputs")
    if policy.magnetic:
        if Path(defs[7]).name != case+".indmfl":
            raise ValueError("launch the primary dmft2 definition; response handles both spins")
        companion = Path("dmft2dn.def")
        if not companion.is_file():
            raise ValueError("magnetic response requires prepared dmft2dn.def")
        down = entries(companion.read_text())
        for unit in [5,9,18,20,30,9919]:
            if unit in defs and (unit not in down or Path(defs[unit]).resolve()!=Path(down[unit]).resolve()):
                raise ValueError("spin branches must share native KS inputs at unit "+str(unit))
        dh = Path(down[7]).read_text().splitlines()[0].split()
        if len(dh)<4 or int(dh[2])!=1 or abs(int(dh[3]))!=5:
            raise ValueError("second spin requires normalized projector 5 and fixed bands")
        for name in ["combineud"]+(["fermif"] if policy.ferro and policy.fixed_number else []):
            helper = Path(root)/name
            if not helper.is_file() or not os.access(helper,os.X_OK):
                raise ValueError("missing magnetic helper: "+str(helper))
'''+s[b:]
s=s.replace('from lapw_total_stress import diagonal_invariant_strains\n\n        diagonal_invariant_strains(s)', 'from lapw_total_stress import invariant_strains\n\n        invariant_strains(s)')
s=s.replace('        "backend": backend,\n        "wien_root":','        "backend": backend,\n        "physical_policy": policy.record(),\n        "wien_root":',1)
s=s.replace('p.name in {"EF.dat", "projectorw.dat", "Eorb_imp.dat"}', 'p.name in {"EF.dat", "projectorw.dat", "Eorb_imp.dat", "params.dat"}')
s=s.replace('re.fullmatch(r"sig\\.inp\\d*", p.name)', 're.fullmatch(r"sig\\.inp\\d*(?:dn)?", p.name)')
s=s.replace('        ("lapw1.def", Path("lapw1.def")),\n    ]:', '        ("lapw1.def", Path("lapw1.def")),\n    ] + ([("dmft2dn.def", Path("dmft2dn.def"))] if solver["physical_policy"]["magnetic"] else []):',1)
s=s.replace('        "response_parallel",','        "response_parallel",\n        "lapw_response_policy",\n        "dmft2_spin",\n        "lapw_magnetic_response",\n        "lapw_lattice",',1)
s=s.replace('            "case": case,\n            "reference": str(snapshot),','            "case": case,\n            "physical_policy": solver["physical_policy"],\n            "fermif": str(Path(root)/"fermif"),\n            "reference": str(snapshot),',1)
s=s.replace('        atomic_json(work / "prepare.json", cfg)','        if solver["physical_policy"]["magnetic"]:\n            cfg["dmft2dn_definition"] = str(snapshot/"dmft2dn.def")\n        atomic_json(work / "prepare.json", cfg)',1)
s=s.replace('or report.get("full_tensor_within_symmetry") is not True\n        ):','or report.get("full_tensor_within_symmetry") is not True\n            or report.get("diagnostic_only", False)\n        ):',1)
s=s.replace('k not in {"lapw1.def", "dmft2.def"}','k not in {"lapw1.def", "dmft2.def", "dmft2dn.def"}')
s=s.replace('        with execution_definition(definition) as native:\n            run_core(root, mpi, native, output)', '''        magnetic = due and Path(case+".indmfldn").is_file() and Path(case+".indmfldn").stat().st_size > 0
        if magnetic:
            from lapw_response_policy import detect_policy
            from dmft2_spin import run_channels, run_fermif, average_charge
            policy = detect_policy(".",case,control=req["files"][2],
                                   override=req["config"].get("solver",{}).get("physical_policy"))
            run_channels([definition,"dmft2dn.def"],policy,
                lambda label,df: run_core(root,mpi,df,output),
                lambda bands: run_fermif(Path(root)/"fermif",bands,output))
            average_charge(case,root,output)
            status["physical_policy"] = policy.record()
        else:
            with execution_definition(definition) as native:
                run_core(root, mpi, native, output)''')
p.write_text(s)
p=r/'x_dmft.py';s=p.read_text();needle='''        if not options.def_only:
            #checkSigmas(cixs)
            runExternal('dmft'+idmf, dmfe.ROOT, dmfe.MPI2, options.nom, options.ntail, 'dmft2', options.m_ext, response_config=options.response_config, no_response=options.no_response)'''
replace='''        if not options.def_only:
            # A requested total response owns the two-spin measurement. Prepare
            # and split the companion here, where ordinary input generation is
            # already authorized. --execute-only deliberately never does this.
            if (not options.m_ext and mode == 'c' and
                os.path.isfile(w2k.case+'.indmfldn') and
                dmft2_step.request(w2k.case, 'dmft2.def', options.response_config, options.no_response)['enabled']):
                down = indmffile.Indmfl(w2k.case, 'indmfldn')
                down.read()
                PrepareDefinitionFile_dmft2('2dn', mode, w2k.case, list(down.siginds.keys()), updn, dnup, so, para, w2k.SCRATCH, cmplx, 'dn')
                PrepareInFile('2dn', mode, w2k.case, updn, 0, 0, options.mixEF, options.WL, Temperature, options.Q_ETOT)
                subprocess.check_call([sys.executable, os.path.join(dmfe.ROOT,'ssplit.py'), '-n', str(options.nom), '-t', str(options.ntail), '-l', 'dn'])
            runExternal('dmft'+idmf, dmfe.ROOT, dmfe.MPI2, options.nom, options.ntail, 'dmft2', options.m_ext, response_config=options.response_config, no_response=options.no_response)'''
assert needle in s;s=s.replace(needle,replace);p.write_text(s)
p=r/'run_dmft.py';s=p.read_text().replace('    The current backend publishes potential tangents, not a corrected tensor.','    The executor detects both DMFT channels, establishes their common mu,\n    averages charge once and computes one shared complete stress tensor.');p.write_text(s)
