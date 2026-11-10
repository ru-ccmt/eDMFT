from material_augmentation import *
import ctypes
lib=ctypes.CDLL(str(R/'bessel_native.so'));fn=lib.native_bessel
ptr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');out=np.ctypeslib.ndpointer(dtype=np.float64,flags='F_CONTIGUOUS')
fn.argtypes=[ctypes.c_int,ctypes.c_int,ptr,out,out];fn.restype=None
for case in ['FeSe','SrVO3']:
 mat=Material(case);ev,c,labels=mat.cv[0];npw=len(c)-mat.nlo_total
 with (mat.source/'native.k00000001.bin').open('rb') as f:f.seek(44);kf=np.fromfile(f,'f8',3)
 qn=np.linalg.norm(2*np.pi*(labels[:,:npw].T+kf)/np.array([mat.s.a,mat.s.b,mat.s.c]),axis=1)
 for ia in range(mat.s.nat):
  x=np.ascontiguousarray(qn*mat.s.rmt[ia]);j=np.empty((5,npw),order='F');dj=j.copy(order='F');fn(npw,4,x,j,dj)
  j0=np.array([spherical_jn(l,x) for l in range(5)]);d0=np.array([spherical_jn(l,x,derivative=True) for l in range(5)])
  print(case,ia,'jmax',np.max(abs(j-j0)),'djmax',np.max(abs(dj-d0)))
