from pathlib import Path
import sys,json,numpy as np
root=Path(__file__).parent;sys.path.insert(0,str(root/'stage'))
from frozen_sigma_tangent import CachedFrozenSigmaResponse
rows=[]
for suffix in ['', '_dn']:
 m,_=CachedFrozenSigmaResponse.from_export(root/'native/MnTe/dmft'/('response'+suffix),diagnostic_gap_cutoff=1e-10,cache_mb=0,precontract=False)
 h=np.zeros((len(m.eps),m.eps.shape[1],m.eps.shape[1]),complex)
 rho=m.density(h,t=0,u=m.uc)
 rows.append(dict(spin=suffix or 'up',physical_valence_N=float(.5*np.trace(rho,axis1=1,axis2=2).real.sum()),mu_Ry=m.meta.mu))
result=dict(scope='central frozen-Sigma native exports; branch factor 2 divided out',rows=rows,spin_imbalance=rows[0]['physical_valence_N']-rows[1]['physical_valence_N'])
(root/'spin_compensation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
