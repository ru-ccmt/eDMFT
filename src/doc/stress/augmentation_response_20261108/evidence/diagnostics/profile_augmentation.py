from material_augmentation import *
import time
from lapw_augmentation_response import harmonics
reports=[]
for case,ia in [('FeSe',0),('SrVO3',1)]:
 mat=Material(case);s=mat.s;l=2;latom=sum(s.mult[:ia]);ev,c,labels=mat.cv[0];npw=len(c)-mat.nlo_total
 n=s.jrj[ia];r=s.r0[ia]*np.exp(np.log(s.rmt[ia]/s.r0[ia])/(n-1)*np.arange(n));dv=.07*r*np.sin(1.3*r)
 for ll in range(5):mat.radial(ia,ll,dv)
 mat.coefficients(ia,latom,0);channel=mat.channels[ia,latom,0,l];ctx=mat.context[ia,l]
 with (mat.source/'native.k00000001.bin').open('rb') as f:f.seek(44);kf=np.fromfile(f,'f8',3)
 q=2*np.pi*(labels.T+kf)/np.array([s.a,s.b,s.c]);qn=np.linalg.norm(q,axis=1)
 qloc=(s.rotloc[ia]@mat.rotations[latom]@q.T).T
 geometry=4*np.pi/np.sqrt(s.a*s.b*s.c)*1j**l*harmonics(l,qloc).conj()*np.exp(2j*np.pi*(labels.T@s.pos[latom]))
 rmt=s.rmt[ia];j=spherical_jn(l,qn[:npw]*rmt);dj=qn[:npw]*spherical_jn(l,qn[:npw]*rmt,derivative=True)
 du=ctx['dprimitive'][-1];u=ctx['primitive'][-1];nr=u.shape[1]
 def direct():
  dm=np.zeros((nr,len(c)))
  if ctx['lapw']:
   dm[0,:npw]=rmt*(dj*du[0,1]-j*du[2,1]);dm[1,:npw]=rmt*(j*du[2,0]-dj*du[0,0])
  else:dm[0,:npw]=-rmt*j*du[0,0]/u[0,0]**2
  for slot in range(ctx['lo'].shape[1]):dm[:,npw+mat.lo_indices[ia,l,slot]]=rmt*rmt*ctx['dlo'][:,slot,None]
  return ((geometry[:,None,:]*dm[None,:,:]).reshape((2*l+1)*nr,len(c))@c).reshape(2*l+1,nr,len(ev))
 def cached():return channel.action(u,du,ctx['dlo'],lapw=ctx['lapw'])
 error=np.max(abs(direct()-cached()));np.testing.assert_allclose(direct(),cached(),atol=1e-13,rtol=1e-11)
 def elapsed(fn):
  values=[]
  for _ in range(7):
   t=time.perf_counter()
   for repeat in range(200):fn()
   values.append((time.perf_counter()-t)/200)
  return float(np.median(values))
 td,tc=elapsed(direct),elapsed(cached)
 row=dict(material=case,npw=npw,nbands=len(ev),direct_seconds=td,cached_seconds=tc,speedup=td/tc,maximum_difference=float(error),cache_bytes=channel.fj.nbytes+channel.dfj.nbytes+channel.local.nbytes)
 reports.append(row);print(row,flush=True)
(R/'profile.json').write_text(json.dumps(reports,indent=2)+'\n')
