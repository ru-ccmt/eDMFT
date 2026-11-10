from pathlib import Path
import sys,json,time
import numpy as np
from scipy.special import spherical_jn
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'stage'),'/Users/haule/projects/eDMFT_new/src/python','/private/tmp/edmft_augmentation_response_20261108']
from material_augmentation import Material
# material_augmentation changes sys.path for its standalone diagnostic; keep
# the staged ABI and spherical assembly first for this validation.
sys.path.insert(0,str(ROOT/'stage'))
from lapw_radial_native import NativeRadialContext
from lapw_spherical_response import PreparedSphericalBand,prepare_wronskian_band
from lapw_augmentation_response import PreparedAugmentationChannel
from lapw_density_response import herm
case=sys.argv[1];mat=Material(case);s=mat.s;ev,c,labels=mat.cv[0]
npw=len(c)-mat.nlo_total;nb=c.shape[1]
with (mat.source/'native.k00000001.bin').open('rb') as f:f.seek(44);kf=np.fromfile(f,'f8',3)
q=2*np.pi*(labels.T+kf)/np.array([s.a,s.b,s.c]);qn=np.linalg.norm(q,axis=1)
ltop=int((mat.source/(case+('.in1c' if mat.complex else '.in1'))).read_text().splitlines()[1].split()[1])
h=np.zeros((nb,nb),complex);o=h.copy();latom=0
rng=np.random.default_rng(156);hb=herm(rng.normal(size=(nb,nb))+1j*rng.normal(size=(nb,nb)))
ob=herm(rng.normal(size=(nb,nb))+1j*rng.normal(size=(nb,nb)))
reverse=0.;timings=[]
for a in range(s.nat):
 for l in range(ltop+1):
  n=s.jrj[a];r0=s.r0[a];rmt=s.rmt[a];dx=np.log(rmt/r0)/(n-1)
  r=r0*np.exp(dx*np.arange(n));dv=.07*r*np.sin(1.3*r)
  e=mat.header[a][l];lapw=bool(e<150);e=e if lapw else e-200
  elo=mat.lo[a][l,:mat.ilo[a][l]] if l<4 else np.array([])
  with NativeRadialContext(mat.v[a][0,0],r0=r0,dx=dx,rmt=rmt,energy=e,lo_energies=elo,
             ell=l,charge=s.Znuc[a],relativistic=bool('RELA' in s.mode),lapw=lapw,library=ROOT/'liblapw_response.so') as ctx:
   groups=[npw+mat.lo_indices[a,l,j] for j in range(len(elo))]
   b=np.array([spherical_jn(l,qn[:npw]*rmt),qn[:npw]*spherical_jn(l,qn[:npw]*rmt,derivative=True)])
   dp,dm=ctx.action(dv);gp=np.zeros_like(dp);gm=np.zeros_like(dm);vb=np.zeros(n)
   for site in range(latom,latom+s.mult[a]):
    phase=np.exp(2j*np.pi*(labels.T@s.pos[site]))
    ch=PreparedAugmentationChannel(l,rmt,s.a*s.b*s.c,q,qn,phase,c,npw,groups,matching_values=b)
    triangle=prepare_wronskian_band(l,rmt,s.a*s.b*s.c,q[:npw],phase[:npw],c[:npw],b)
    op=PreparedSphericalBand(ch,ctx.primitive,ctx.lo_coefficients,ctx.overlap_metric,
              ctx.radius,e,elo,triangle,lapw=lapw)
    start=time.perf_counter();dh,do=op.action(dp,dm,difference_dot=-dv);timings.append(time.perf_counter()-start)
    h+=dh;o+=do
    pbar,mbar,eb,elb,db=op.pullback(hb,ob);gp+=pbar;gm+=mbar;vb-=db
   vb+=ctx.pullback(gp,gm)[0];reverse+=vb@dv
 latom+=s.mult[a];print(case,'atom',a,'done',flush=True)
old=Path('/private/tmp/edmft_augmentation_response_20261108')/case/'moving_probe'
reference=np.load(old/'complete_matrices.npz');ns=np.load(old/'central_response.npz')['h_nonspherical']
for a,b in [(h,reference['H']-ns),(o,reference['O'])]:np.testing.assert_allclose(a,b,rtol=1e-7,atol=1e-10)
forward=np.vdot(hb,h).real+np.vdot(ob,o).real
report=dict(case=case,H_relative=float(np.linalg.norm(h-(reference['H']-ns))/np.linalg.norm(h)),
  O_relative=float(np.linalg.norm(o-reference['O'])/np.linalg.norm(o)),transpose=float(abs(forward-reverse)/max(1.,abs(forward),abs(reverse))),
  forward=forward,reverse=reverse,action_seconds=float(sum(timings)),channels=len(timings))
(ROOT/(case+'_spherical.json')).write_text(json.dumps(report,indent=2)+'\n');print(report)
