import sys
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from wstruct import Struct
from wien_response_density import Density
from setup import ROOT,SOURCES

def integrate(d,struct):
 s=Struct();s.ReadStruct(str(struct));vol=s.a*s.b*s.c;cell=np.array([s.a,s.b,s.c]);mt=0
 for (at,l,m),v in d.rad.items():
  if l==m==0:x=d.grids[at];mt+=d.mult[at]*CubicSpline(x,v*np.exp(x)).integrate(x[0],x[-1])
 ir=0j
 for g,z in d.pw.items():
  g=np.array(g);k=np.linalg.norm(2*np.pi*g/cell);th=complex(np.all(g==0));ipos=0
  for at in range(s.nat):
   x=k*s.rmt[at];form=1 if x==0 else 3*spherical_jn(1,x)/x
   for j in range(s.mult[at]):th-=4*np.pi*s.rmt[at]**3/(3*vol)*form*np.exp(-2j*np.pi*np.dot(s.pos[ipos],g));ipos+=1
  ir+=z*th.conjugate()*vol
 return mt,ir.real,mt+ir.real
if __name__=='__main__':
 for case in ['SrVO3','FeSe']:
  for name in ['central','source_p_0.0002500','source_m_0.0002500']:
   p=ROOT/case/name
   for suf in ['clmsum','clmval','clmcor','clmout','qout']:
    if not (p/(case+'.'+suf)).exists():continue
    d=Density(p/(case+'.'+suf),p/(case+'.struct'));ss=SOURCES[case]/(case+'.struct') if suf=='qout' else p/(case+'.struct');print(case,name,suf,integrate(d,ss))
