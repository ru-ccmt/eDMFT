from pathlib import Path
import sys,re,json
import numpy as np
from scipy.interpolate import CubicSpline
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from wstruct import Struct
from scipy.special import spherical_jn
r=Path(__file__).parent/'SrVO3_affine';out={}
for name in ['center','p001','p0005','p00025']:
 p=r/name;ss=Struct();ss.ReadStruct(str(p/'SrVO3.struct'));a=(p/'SrVO3.clmsum').read_text().splitlines();i=0;atom=-1;mt=0.
 while i<len(a):
  if 'ATOMNUMBER' in a[i]:atom=int(a[i].split('=')[1])-1
  if not re.search(r'CLM\(R\) FOR L\s+0\s+M=\s*0',a[i]):i+=1;continue
  i+=1
  while not a[i].strip():i+=1
  v=[]
  while len(v)<ss.jrj[atom]:v.extend(map(float,re.findall(r'[-+]?\d+\.\d+E[-+]\d+',a[i])));i+=1
  x=np.linspace(np.log(ss.r0[atom]),np.log(ss.rmt[atom]),ss.jrj[atom]);mt+=ss.mult[atom]*CubicSpline(x,np.array(v)*np.exp(x)).integrate(x[0],x[-1])
 j=next(i for i,x in enumerate(a) if 'NUMBER OF PW' in x);n=int(a[j].split()[0]);ir=0j;irc=0j;vol=ss.a*ss.b*ss.c
 for line in a[j+1:j+1+n]:
  g=np.array([int(line[j:j+5]) for j in [3,8,13]]);z=float(line[18:37])+1j*float(line[37:56]);q=np.sqrt(np.sum((2*np.pi*g/np.array([ss.a,ss.b,ss.c]))**2));th=complex(int(np.all(g==0)));posindex=0
  for at in range(ss.nat):
   x=q*ss.rmt[at];f=1 if x==0 else 3*spherical_jn(1,x)/x
   for mu in range(ss.mult[at]):
    th-=4*np.pi*ss.rmt[at]**3/(3*vol)*f*np.exp(-2j*np.pi*np.dot(ss.pos[posindex],g));posindex+=1
  ir+=z*th.conjugate()*vol;irc+=z*th*vol
 out[name]={'MT_charge':mt,'interstitial_charge':ir.real,'total':mt+ir.real,'alternate_conjugation_total':mt+irc.real}
print(json.dumps(out,indent=2));(r/'charge_checks.json').write_text(json.dumps(out,indent=2))
