from pathlib import Path
import sys,json,numpy as np
root=Path(__file__).parent;sys.path.insert(0,str(root/'stage'))
from lapw_charge_response import NativeChargeNormalization,native_mixer_inversion
from lapw_lattice import native_structure
from wien_response_density import Density
from wstruct import Struct
rows=[]
for item in json.loads((root/'native_jobs.json').read_text()):
 w=Path(item['directory']);c=item['case'];s=Struct();s.ReadStruct(str(w/(c+'.struct')));s=native_structure(s)
 inp=Density(w/(c+'.clmsum'),w/(c+'.struct'));raw=Density(w/(c+'.clmval'),w/(c+'.struct'))
 if item['material']=='MnTe':raw.add(Density(w/(c+'.clmvaldn'),w/(c+'.struct'))).scale(.5)
 raw.add(Density(w/(c+'.clmcor'),w/(c+'.struct')))
 normalization=NativeChargeNormalization(raw,s,float(np.dot(s.mult,s.Znuc)),inversion=native_mixer_inversion(s),enabled=((w/(c+'.inm')).read_text().split()[2].upper()=='YES'))
 delta=normalization.central-inp.vector()
 rows.append(dict(material=item['material'],directory=item['name'],residual_norm=float(np.linalg.norm(delta)),maximum=float(abs(delta).max()),alpha=float(normalization.alpha),normalization_enabled=normalization.enabled))
print(json.dumps(rows,indent=2));(root/'fixedpoint_probe.json').write_text(json.dumps(rows,indent=2)+'\n')
