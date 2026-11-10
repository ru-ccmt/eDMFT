from pathlib import Path
import re,ctypes as ct,subprocess,json
import numpy as np
from scipy.integrate import simpson
root=Path(__file__).resolve().parent
src=Path('/Users/haule/wien2k_24/SRC_lapw2/kin_str_integrals.F')
body=re.search(r'(?ims)^subroutine g1_radial_derivative.*?^end subroutine g1_radial_derivative',src.read_text()).group()
(root/'radial_probe.f90').write_text(body+'\n')
subprocess.run(['/opt/homebrew/bin/gfortran','-O0','-g','-fcheck=all','-ffree-form','-shared','-fPIC','radial_probe.f90','/Users/haule/wien2k_24/SRC_lapw2/spline.f','/Users/haule/wien2k_24/SRC_lapw2/dergl.f','-o','radial_probe.dylib'],cwd=root,check=True)
lib=ct.CDLL(str(root/'radial_probe.dylib'))
def deriv(r,f,parity=None):
 n=np.array([len(r)],np.int32);out=np.zeros(len(r));f=np.array(f,float);r=np.array(r,float)
 args=[n,r,f,out]
 if parity is not None:args+=[np.array([parity],np.int32),np.zeros(1)]
 fun=lib.g1_radial_derivative_ if parity is None else lib.dergl_
 fun(*[a.ctypes.data_as(ct.c_void_p) for a in args]);return out
out={}
r=np.geomspace(1e-5,2.,801)
for power in [0,1,2,3]:
 f=r**power;df=power*r**(power-1) if power else np.zeros_like(r)
 d=deriv(r,f);ref2=power*(power-1)*r**(power-2) if power>1 else np.zeros_like(r)
 dd=deriv(r,d)
 # Norms are weighted as the actual G1 integrals, not singular relative
 # errors of a vanishing derivative near the origin.
 e1=float(simpson(r*abs(d-df),x=r));e2=float(simpson(r*r*abs(dd-ref2),x=r))
 out['power_'+str(power)]={'first_weighted_abs_error':e1,'second_weighted_abs_error':e2}
 assert e1<1e-8 and e2<1e-7
out['constant_old_odd_max_abs']=float(max(abs(deriv(r,np.ones_like(r),0))))
out['constant_positive_max_abs']=float(max(abs(deriv(r,np.ones_like(r)))))
# Cusp-containing l=0 and a noninteger RELA-like exponent.
for name,p,z in [('cusp',0.,2.7),('fractional',.83,3.1)]:
 res=[]
 for n in [401,801,1601]:
  r=np.geomspace(1e-5,2.,n);f=r**p*np.exp(-z*r)
  d=f*(p/r-z);dd=f*((p/r-z)**2-p/r**2)
  got=deriv(r,f);got2=deriv(r,got)
  e1=simpson(r*abs(got-d),x=r);e2=simpson(r*r*abs(got2-dd),x=r)
  res.append({'n':n,'first_error':float(e1),'second_error':float(e2)})
 assert res[-1]['second_error']<res[0]['second_error']/10
 out[name]=res
(root/'radial_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
