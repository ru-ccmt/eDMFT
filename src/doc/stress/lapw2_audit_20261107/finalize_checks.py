from pathlib import Path
import json,re
from run_cases import ROOT,CASES
folders=[]
for case,source,_ in CASES:
    folder=ROOT/'runs'/(case+'_'+source.name)
    assert (folder/'audit/relaxation.json').exists()
    scfm=folder/'audit/corrected_total'/(case+'.scfm')
    text=scfm.read_text()
    dis=float(re.findall(r'^:DIS.*$',text,re.M)[-1].split()[-1])
    assert dis<=2e-7,(folder,dis)
    assert len(re.findall(r'^:STRESS_GPa\d+:.*total',text,re.M))==3
    for program in ['lapw0','lapw1','lapw2','lcore','mixer']:
        error=folder/(program+'.error')
        assert not error.exists() or not error.read_text().strip(),error
    for ext in ['scf0','scf1','scf2','scfc','scfm']:
        p=folder/(case+'.'+ext)
        assert not re.search(r'\b(?:nan|infinity|inf)\b',p.read_text(),re.I),p
    folders.append(str(folder))
(ROOT/'relaxed_runs.json').write_text(json.dumps(folders,indent=2)+'\n')
print('Six completed, converged DFT chains: no error records or nonfinite SCF values.')
