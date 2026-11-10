from pathlib import Path
import sys,json,re,numpy as np
root=Path(__file__).parent;sys.path.insert(0,str(root/'stage'))
from frozen_sigma_tangent import CachedFrozenSigmaResponse
from frozen_sigma_klein import last_tag_value
from wstruct import Struct
from lapw_lattice import NativeLattice
out={}
for mat,dirs,c,number in [('MnO',['MnO_center','MnO_m_0.001','MnO_p_0.001'],'MnO',21),('MnTe',['dmft','dmft_m_0.001','dmft_p_0.001'],'0_800_MnTe',62)]:
 vals=[]
 for d in dirs:
  w=root/'native'/mat/d;s=Struct();s.ReadStruct(str(w/(c+'.struct')))
  lines=(w/(c+'.inc')).read_text().splitlines();i=0;occup=[]
  for a in range(s.nat):
   n=int(lines[i].split()[0]);i+=1
   for j in range(n):occup.append(float(lines[i+j].replace(',',' ').split()[2])*s.mult[a])
   i+=n
  energies=[float(x) for x in re.findall(r'^:[^:]+:.*?([-+]?\d+\.\d+)\s+Ry', (w/(c+'.scfc')).read_text(),re.M)]
  assert len(occup)==len(energies)
  core=float(np.dot(occup,energies));den=last_tag_value(w/(c+'.scf0'),':DEN')
  ls=[]
  for suffix in (['','_dn'] if mat=='MnTe' else ['']):
   m,_=CachedFrozenSigmaResponse.from_export(w/('response'+suffix),diagnostic_gap_cutoff=1e-10,cache_mb=0,precontract=False)
   dh=np.zeros((len(m.eps),m.eps.shape[1],m.eps.shape[1]),complex)
   ls.append(m.lattice_scalar(dh,t=0,u=m.uc)+m.meta.mu*number)
  value=dict(directory=d,core_ry=core,den_ry=den,lattice_ry=float(np.mean(ls)),K_ry=core+den+float(np.mean(ls)),mu_ry=m.meta.mu,volume=NativeLattice.from_structure(s).volume)
  vals.append(value);print(mat,json.dumps(value),flush=True)
 fd=(vals[2]['K_ry']-vals[1]['K_ry'])/.002*14710.507641078/(3*vals[0]['volume'])
 out[mat]=dict(scope='One native update on saved endpoint densities, NOT newly reconverged endpoints. Independent reconstruction C_DEN+core+L+mu*Ntarget.',values=vals,mean_numeric_stress_GPa=fd)
 (root/'native_scalar_probe.json').write_text(json.dumps(out,indent=2)+'\n')
 print(mat,'mean numeric',fd,flush=True)
