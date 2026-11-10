from pathlib import Path
import re,json,numpy as np,hashlib
r=Path(__file__).resolve().parent

def t(p,tag):
 s=p.read_text();return np.array([list(map(float,re.findall('^:'+tag+f'{i:03d}:'+r'\s*(.*)',s,re.M)[-1].split())) for i in [1,2,3]])
b=r/'SrVO3_DMFT/baseline';d=t(b/'SrVO3.scf0','STR_DENS');core=t(b/'SrVO3.scfc','STR_CORE');newd=t(r/'SrVO3_DMFT/stage/SrVO3.scf0','STR_DENS');out={}
for n in [43,100,200,400,680]:
 folder=r/f'SrVO3_bands/results_{n}';val=t(folder/'SrVO3.scf2','STR_VAL');m=(d-core-val)*14710.5164057;newm=(newd-core-val)*14710.5164057;s=(folder/'SrVO3.scf2').read_text()
 out[n]={'diagonal_GPa':m.diagonal().tolist(),'mean_GPa':float(np.trace(m)/3),'with_lapw0_fixes_GPa':float(np.trace(newm)/3),'EF_eV':float((folder/'EF.dat').read_text()),'SUM_Ry':float(re.findall(r'^:SUM.*=\s*([-+0-9.]+)',s,re.M)[-1]),'VSUM_Ry':float(re.findall(r'^:VSUM.*=\s*([-+0-9.]+)',s,re.M)[-1]),'SOlapm_sha256':hashlib.sha256((folder/'SOlapm.dat').read_bytes()).hexdigest()}
 print(n,out[n])
(r/'band_results.json').write_text(json.dumps(out,indent=2)+'\n')
