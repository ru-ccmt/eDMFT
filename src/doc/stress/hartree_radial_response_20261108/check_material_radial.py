from pathlib import Path
import sys,ctypes,json,time
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_hartree_radial import NativeCoulombIntegral,RadialHartreeDirichlet
from wien_response_density import Density
ref=ctypes.CDLL(str(ROOT/'native/libhartree_reference.dylib'));ptr=ctypes.POINTER(ctypes.c_double)
report={}
for case,label in [('FeSe','FeSe_b509_one'),('SrVO3','SrVO3_b703_one')]:
 src=Path('/private/tmp/edmft_coupled_valence_response_20261108')/label/'central'
 template=Density(src/(case+'.clmval'),src/(case+'.struct'))
 dn=np.load(Path('/private/tmp/edmft_valence_charge_checks_20261108')/(case+'_joint_fixedN.npz'))['density']
 density=template.from_vector(dn);integrations={a:NativeCoulombIntegral(np.exp(g)) for a,g in enumerate(template.grids)}
 rows=[];actions=[]
 for (a,ll,m),c in density.rad.items():
  l=abs(ll);op=RadialHartreeDirichlet(integrations[a],l);r=integrations[a].r;n=len(r);dx=integrations[a].dx;bdy=.02
  out=np.empty(n)
  ref.radial_dirichlet_ref_(ctypes.byref(ctypes.c_int(n)),r.ctypes.data_as(ptr),ctypes.byref(ctypes.c_double(dx)),ctypes.byref(ctypes.c_int(l)),c.ctypes.data_as(ptr),ctypes.byref(ctypes.c_double(bdy)),out.ctypes.data_as(ptr))
  actual=op.action(c,bdy);error=np.linalg.norm(actual-out)/max(np.linalg.norm(out),1e-30)
  rows.append(dict(atom=a,ell=ll,m=m,nrad=n,native_relative=float(error),transpose=op.operator.check_transpose(samples=1,tolerance=1e-11)))
  actions.append((op,c,bdy))
 def timing(reverse):
  times=[]
  for _ in range(7):
   t=time.perf_counter()
   for op,c,b in actions:
    if reverse:op.pullback(c)
    else:op.action(c,b)
   times.append(time.perf_counter()-t)
  return float(np.median(times))
 report[case]=dict(channels=rows,whole_MT_forward_seconds=timing(False),whole_MT_reverse_seconds=timing(True))
 print(case,'harmonics',len(rows),'max native error',max(v['native_relative'] for v in rows),'transpose',max(v['transpose'] for v in rows),'seconds',report[case]['whole_MT_forward_seconds'],report[case]['whole_MT_reverse_seconds'],flush=True)
(ROOT/'material_radial_hartree.json').write_text(json.dumps(report,indent=2)+'\n')
