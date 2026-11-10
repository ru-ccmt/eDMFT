from pathlib import Path
import re,json,struct,hashlib
import numpy as np
root=Path('/private/tmp/numeric_stress_audit_20261107');out={};case='SrVO3'
def tag(s,t):return float(re.findall('^'+re.escape(t)+r'\s*(.*)$',s,re.M)[-1].split('=')[-1].split()[0])
for key,name in [('minus','V_strain_m001'),('center','V_plus_force'),('plus','V_strain_p001')]:
 d=Path('/Users/haule/projects/materials/SrVO3')/name;s=(d/(case+'.scf')).read_text();s=re.split(r'(?=^:ITE\d+:)',s,flags=re.M)[-1]
 mult=[int(v) for v in re.findall(r'MULT=\s*(\d+)',(d/(case+'.struct')).read_text())];ls=(d/(case+'.inc')).read_text().splitlines();occ=[];i=0
 for m in mult:
  n=int(ls[i].split()[0]);i+=1
  for j in range(n):occ.append(float(ls[i].split('#')[0].split('(')[0].replace(',',' ').split()[2])*m);i+=1
 corelines=[l for l in s.splitlines() if re.match(r'^:[1-7][SPDFGHKL][A-Z ]*\d{3}:',l) and l.endswith('Ry')]
 corevalues=[float(l.split()[-2]) for l in corelines];assert len(corevalues)==len(occ),(key,len(corevalues),len(occ))
 core=sum(a*b for a,b in zip(corevalues,occ));density=tag(s,':DEN  :');vs=tag(s,':VSUM :');K=tag(s,':ENE  :')-tag(s,':SUM  :')+vs-tag(s,':EORB')
 out[key]={'core_eigenvalue_sum_Ry':core,'density_energy_Ry':density,'VSUM_Ry':vs,'independent_K_Ry':density+core+vs,'mixer_based_K_Ry':K,'difference_Ry':density+core+vs-K,'core_electron_count':sum(occ)}
 print(key,out[key])
fac=14710.507641078/(.006*382.70904)
out['independent_mean_GPa']=(out['plus']['independent_K_Ry']-out['minus']['independent_K_Ry'])*fac
out['components_mean_GPa']={k:(out['plus'][k]-out['minus'][k])*fac for k in ['core_eigenvalue_sum_Ry','density_energy_Ry','VSUM_Ry']}
print('independent FD',out['independent_mean_GPa'],'pieces',out['components_mean_GPa'])

def vectors(d):
 nat=4;nlo_values=set(map(int,re.findall(r'number of local orbitals, nlo \(hamilt\)\s+(\d+)',(d/(case+'.output1')).read_text())));assert len(nlo_values)==1;nlo=nlo_values.pop();a=[]
 files=sorted(p for p in d.glob(case+'.vector_*') if re.fullmatch(re.escape(case)+r'\.vector_\d+',p.name))
 for p in files:
  with p.open('rb') as f:
   def rec():
    raw=f.read(4)
    if not raw:return None
    n=struct.unpack('<i',raw)[0];data=f.read(n);assert len(data)==n and f.read(4)==raw;return data
   for i in range(2*nat):rec()
   while (h:=rec()) is not None:
    xx,yy,zz,name,n,nb,w=struct.unpack('<3d10s2id',h);g=rec();assert len(g)==12*n;pw=np.frombuffer(g,dtype='<i4').reshape(n,3)[:n-nlo];digest=hashlib.sha256(np.array(sorted(map(tuple,pw.tolist())),dtype='<i4').tobytes()).hexdigest()
    for j in range(nb):rec();rec()
    a.append(dict(k=[xx,yy,zz],weight=w,basis=n,bands=nb,pw_hash=digest))
 a.sort(key=lambda r:r['k']);return a
v={name:vectors(root/name) for name in ['minus','plus']};out['fresh_vector_consistency']={k:[a[k] for a in v['minus']]==[a[k] for a in v['plus']] for k in v['minus'][0]};out['fresh_vectors']=v
print('actual G sets',out['fresh_vector_consistency'],'nk',len(v['minus']))
(root/'independent_checks.json').write_text(json.dumps(out,indent=2)+'\n')
