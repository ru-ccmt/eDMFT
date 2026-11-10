from pathlib import Path
import subprocess,ctypes,json
import numpy as np
root=Path('/private/tmp/stress_remaining_audit_20261107')
subprocess.run(['gfortran','-dynamiclib','-ffree-form','-O2',str(root/'SRC_lapw0/gga_metric_tensor.f'),'-o',str(root/'metric.dylib')],check=True)
lib=ctypes.CDLL(str(root/'metric.dylib'));fn=lib.gga_metric_tensor_
p=ctypes.POINTER(ctypes.c_double)
def kernel(js,gu,gd,su,sd,sud):
 q=np.zeros((3,3),order='F');gu=np.asarray(gu,dtype=float);gd=np.asarray(gd,dtype=float);j=ctypes.c_int(js);s=[ctypes.c_double(x) for x in [su,sd,sud]]
 fn(ctypes.byref(j),gu.ctypes.data_as(p),gd.ctypes.data_as(p),*[ctypes.byref(x) for x in s],q.ctypes.data_as(p));return q
rng=np.random.default_rng(7011);maxerr=0;maxrot=0;maxspin=0
for it in range(60):
 gu=rng.normal(size=3);gd=rng.normal(size=3);a,b,c=rng.normal(size=3)
 q=kernel(2,gu,gd,a/2,b/2,c/2)
 for ia in range(3):
  for ib in range(ia,3):
   e=np.zeros((3,3));e[ia,ib]=e[ib,ia]=.5 if ia!=ib else 1
   h=2e-5
   def energy(t):
    u=np.linalg.solve(np.eye(3)+t*e,gu);d=np.linalg.solve(np.eye(3)+t*e,gd)
    return a*u@u+b*d@d+c*u@d
   fd=(-energy(2*h)+8*energy(h)-8*energy(-h)+energy(-2*h))/(12*h)
   maxerr=max(maxerr,abs(fd+q[ia,ib]))
 rot=np.linalg.qr(rng.normal(size=(3,3)))[0]
 maxrot=max(maxrot,float(np.max(abs(kernel(2,rot@gu,rot@gd,a/2,b/2,c/2)-rot@q@rot.T))))
 # An unpolarized functional a*|grad n|^2 has spin derivatives a,a,2a.
 maxspin=max(maxspin,float(np.max(abs(kernel(2,gu,gd,a/2,a/2,a)-kernel(1,gu,gd,a,a,0)))))
 # Exchange of the spin labels must leave the tensor invariant.
 assert np.max(abs(q-kernel(2,gd,gu,b/2,a/2,c/2)))<1e-13
assert maxerr<1e-8 and maxrot<1e-12 and maxspin<1e-12
res={'six_strain_FD_max_error':maxerr,'rotation_covariance_max_error':maxrot,'spin_representation_max_error':maxspin,'cases':60,'status':'PASS'}
(root/'metric_test_results.json').write_text(json.dumps(res,indent=2)+'\n');print(res)
# Independent Cartesian polynomial oracle for the compiled lcore fixture.
lines=(root/'core_fixture/output_fixed.txt').read_text().splitlines();old=(root/'core_fixture/output.txt').read_text().splitlines();error=0;before={}
for l in old:
 if l.startswith('FIXTURE'):
  v=l.split();before[int(v[1])]=np.array(v[4:],float).reshape(3,3)
for l in lines:
 if not l.startswith('FIXTURE'):continue
 v=l.split();mode=int(v[1]);obs=np.array(v[4:],float).reshape(3,3);Q=np.zeros((3,3))
 if mode==0:Q=np.eye(3)
 if mode==1:Q=np.diag([-1,-1,2])*np.sqrt(5/(16*np.pi))
 if mode in [2,4]:
  a=0 if mode==2 else 1;Q[a,2]=Q[2,a]=np.sqrt(15/(4*np.pi))/2
 if mode==3:Q=np.diag([1,-1,0])*np.sqrt(15/(16*np.pi))
 if mode==5:Q[0,1]=Q[1,0]=np.sqrt(15/(4*np.pi))/2
 ref=2*Q/15;error=max(error,float(np.max(abs(obs-ref))))
 print('core fixture',mode,'old error',float(np.max(abs(before[mode]-ref))),'fixed error',float(np.max(abs(obs-ref))))
assert error<2e-8
(root/'core_fixture/results.json').write_text(json.dumps({'fixed_max_error':error,'status':'PASS','old_missing_yz_xy':True},indent=2)+'\n')
