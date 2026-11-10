from pathlib import Path
import sys,re,subprocess,json
import numpy as np
from scipy.io import FortranFile
from scipy.special import sph_harm_y,spherical_jn
from numpy.polynomial.legendre import leggauss
R=Path(__file__).resolve().parent
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
sys.path.insert(0,str(R/'stage'))
sys.path.insert(0,'/private/tmp/edmft_radial_all_l_20261108')
from wstruct import Struct
from lapw_augmentation_response import PreparedAugmentationChannel,RadialAngularLayout,AugmentationContractions
from test_wien_matrices import vectors
from check_oneparticle_stress_fd import read_basis_matrix
FLOAT=re.compile(r'[-+]?\d+\.\d+[DE][-+]\d+')

def read_vns(path):
 blocks=re.split(r'ATOMNUMBER\s*=\s*\d+',path.read_text())[1:];out=[]
 for b in blocks:
  channels=re.split(r'VLM\(R\) FOR L\s*(-?\d+)\s+M=\s*(\d+)',b)
  out.append({(int(channels[i]),int(channels[i+1])):np.array([float(x.replace('D','E')) for x in FLOAT.findall(channels[i+2])]) for i in range(1,len(channels),3)})
 return out

def angular_harmonics(q,lmax):
 lengths=np.linalg.norm(q,axis=-1);z=np.divide(q[:,2],lengths,out=np.ones_like(lengths),where=lengths>0)
 theta=np.arccos(np.clip(z,-1,1));phi=np.arctan2(q[:,1],q[:,0]);out=[]
 for l in range(lmax+1):
  for m in range(-l,l+1):out.append(sph_harm_y(l,m,theta,phi))
 return np.array(out)

def gaunt_real(lmax,L,M):
 # Product quadrature exactly integrates these finite harmonic products.
 x,w=leggauss(24);ph=2*np.pi*np.arange(49)/49;th=np.arccos(x)
 q=np.array([np.sin(th)[:,None]*np.cos(ph),np.sin(th)[:,None]*np.sin(ph),np.broadcast_to(x[:,None],(24,49))]).reshape(3,-1).T
 y=angular_harmonics(q,max(lmax,abs(L)));yl=y[abs(L)*(abs(L)+1)+M]
 if M:
  minus=y[abs(L)*(abs(L)+1)-M]
  potential=(((-1)**M)*yl+minus)/np.sqrt(2) if L>=0 else -1j*(((-1)**M)*yl-minus)/np.sqrt(2)
 else:potential=yl
 y=y[:(lmax+1)**2]
 return (y.conj()*(np.repeat(w,49)*2*np.pi/49*potential))@y.T

class Material:
 def __init__(self,case,source=None):
  self.case=case;self.complex=case=='SrVO3'
  self.source=source or Path('/private/tmp/edmft_native_operator_20261108/native_probe')/case/'center'
  self.work=R/case;self.work.mkdir(exist_ok=True)
  self.s=Struct();self.s.ReadStruct(str(self.source/(case+'.struct')));s=self.s
  self.v=read_vns(self.source/(case+'.vsp'));self.vns=read_vns(self.source/(case+'.vns'))
  for fields in [self.v,self.vns]:
   for ia,channels in enumerate(fields):
    for key in channels:channels[key]=channels[key][:s.jrj[ia]]
  self.header=[];self.lo=[]
  with FortranFile(self.source/(case+'.vector')) as f:
   for ia in range(s.nat):
    self.header.append(f.read_record('f8'));self.lo.append(f.read_record('f8').reshape(-1,4).T)
  self.context={};self.primitive={};self.rotations=[];self.channels={}
  start=0
  for ia in range(s.nat):
   for pos in s.pos[start:start+s.mult[ia]]:
    for rotation,tau in zip(s.timat,s.tau):
     diff=rotation@s.pos[start]+tau-pos;diff-=np.rint(diff)
     if np.max(abs(diff))<2e-7:
      self.rotations.append(rotation.T);break
    else:raise ValueError('missing rotation')
   start+=s.mult[ia]
  self.ilo=[np.sum(e<995,axis=1) for e in self.lo]
  self.lo_indices={};index=0
  for ia in range(s.nat):
   for l in range(4):
    for j in range(self.ilo[ia][l]):
     count=s.mult[ia]*(2*l+1);self.lo_indices[ia,l,j]=np.arange(index,index+count);index+=count
  self.nlo_total=index
  self.cv=list(vectors(self.source/(case+'.vector'),s.nat,self.complex))

 def radial(self,ia,l,dv=None,de=0):
  s=self.s;n=s.jrj[ia];r0=s.r0[ia];rmt=s.rmt[ia];dx=np.log(rmt/r0)/(n-1);r=r0*np.exp(dx*np.arange(n))
  e=self.header[ia][l];lapw=e<150;e=e if lapw else e-200
  elo=self.lo[ia][l,:self.ilo[ia][l]] if l<4 else np.array([])
  inp=self.work/f'radial_a{ia}_l{l}.input';out=inp.with_suffix('.bin')
  v=self.v[ia][0,0];dv=np.zeros(n) if dv is None else dv
  with inp.open('w') as f:
   f.write(f'{n} {len(elo)} {int("RELA" in s.mode)} {l} {int(lapw)}\n')
   f.write(' '.join(f'{x:.17e}' for x in [r0,dx,rmt,e,s.Znuc[ia],de])+'\n')
   if len(elo):f.write(' '.join(map(str,elo))+'\n'+' '.join(map(str,np.zeros(len(elo))))+'\n')
   for arr in [v,dv]:f.write(' '.join(f'{x:.17e}' for x in arr)+'\n')
  subprocess.run([str(R/'export_radial'),str(inp),str(out)],check=True)
  data=np.fromfile(out,'f8');off=0
  def get(shape):
   nonlocal off
   size=np.prod(shape);a=data[off:off+size].reshape(shape,order='F');off+=size;return a
  primitive=get((n,4,2+len(elo)));du=get(primitive.shape);weight=get((n,4));wp=get((n,4));m=get((2+len(elo),len(elo)));dm=get(m.shape)
  assert off==len(data)
  ctx=dict(primitive=primitive,dprimitive=du,weight=weight,wp=wp,lo=m,dlo=dm,r=r,lapw=lapw,energy=e,lo_energy=elo)
  self.context[ia,l]=ctx;return ctx

 def coefficients(self,ia,latom,ik,lmax=4):
  s=self.s;ev,c,labels=self.cv[ik];n=len(labels.T);npw=n-self.nlo_total
  p=self.source/f'native.k{ik+1:08d}.bin'
  with p.open('rb') as f:f.seek(44);kf=np.fromfile(f,'f8',3)
  q=2*np.pi*(labels.T+kf)/np.array([s.a,s.b,s.c]);qn=np.linalg.norm(q,axis=1)
  qloc=(s.rotloc[ia]@self.rotations[latom]@q.T).T
  y=angular_harmonics(qloc,lmax)
  phase=np.exp(2j*np.pi*(labels.T@s.pos[latom]));volume=s.a*s.b*s.c;rmt=s.rmt[ia]
  zs=[];dzs=[];ls=[];ms=[];rs=[]
  for l in range(lmax+1):
   ctx=self.context[ia,l];u=ctx['primitive'];du=ctx['dprimitive'];nr=u.shape[-1]
   mat=np.zeros((nr,n));dmat=mat.copy();j=spherical_jn(l,qn[:npw]*rmt);dj=qn[:npw]*spherical_jn(l,qn[:npw]*rmt,derivative=True)
   p,d=u[-1,1-1,0]/rmt,u[-1,3-1,0]/rmt;pe,dpe=u[-1,0,1]/rmt,u[-1,2,1]/rmt
   dp,dd=du[-1,0,0]/rmt,du[-1,2,0]/rmt;dpe0,ddpe=du[-1,0,1]/rmt,du[-1,2,1]/rmt
   if ctx['lapw']:
    mat[0,:npw]=rmt*rmt*(dj*pe-j*dpe);mat[1,:npw]=rmt*rmt*(j*d-dj*p)
    dmat[0,:npw]=rmt*rmt*(dj*dpe0-j*ddpe);dmat[1,:npw]=rmt*rmt*(j*dd-dj*dp)
   else:mat[0,:npw]=j/p;dmat[0,:npw]=-j*dp/p**2
   for jlo in range(ctx['lo'].shape[1]):
    indices=npw+self.lo_indices[ia,l,jlo]
    mat[:,indices]=rmt*rmt*ctx['lo'][:,jlo,None];dmat[:,indices]=rmt*rmt*ctx['dlo'][:,jlo,None]
   groups=[npw+self.lo_indices[ia,l,jlo] for jlo in range(ctx['lo'].shape[1])]
   channel=PreparedAugmentationChannel(l,rmt,volume,qloc,qn,phase,c,npw,groups)
   self.channels[ia,latom,ik,l]=channel
   cached=channel.value(u[-1],ctx['lo'],lapw=ctx['lapw'])
   cached_dot=channel.action(u[-1],du[-1],ctx['dlo'],lapw=ctx['lapw'])
   for m in range(-l,l+1):
    geometry=4*np.pi/np.sqrt(volume)*1j**l*y[l*(l+1)+m].conj()*phase
    direct=(mat*geometry)@c;direct_dot=(dmat*geometry)@c
    np.testing.assert_allclose(cached[m+l],direct,rtol=1e-11,atol=2e-12)
    np.testing.assert_allclose(cached_dot[m+l],direct_dot,rtol=1e-11,atol=2e-12)
    zs.extend(cached[m+l]);dzs.extend(cached_dot[m+l])
    ls.extend([l]*nr);ms.extend([m]*nr);rs.extend(range(nr))
  return np.array(zs),np.array(dzs),np.array(ls),np.array(ms),np.array(rs)

 def kernel(self,ia,channels,indices,transpose=False):
  ls,ms,rs=indices;ng=len(ls);s=self.s;n=s.jrj[ia]
  p=np.column_stack([self.context[ia,l]['primitive'][:,0,r] for l,r in zip(ls,rs)])
  small=np.column_stack([self.context[ia,l]['primitive'][:,1,r] for l,r in zip(ls,rs)])
  w=self.context[ia,0]['weight'];out=np.zeros((ng,ng),complex)
  for (L,M),v in channels.items():
   g=gaunt_real(max(ls),L,M);angular=ls*(ls+1)+ms;g=g[np.ix_(angular,angular)]
   if transpose:g=g.T
   integral=p.T@((w[:,0]*v)[:,None]*p)+small.T@((w[:,1]*v)[:,None]*small)
   out+=g*integral
  return out

 def small_kernels(self,ia,channels,projector=None):
  counts={l:self.context[ia,l]['primitive'].shape[-1] for l in range(5)}
  layout=RadialAngularLayout(counts,channels.keys());pairs={};dpairs={};metrics={};dmetrics={}
  v=np.array(list(channels.values())).T
  for l,nl in counts.items():
   ctx=self.context[ia,l];u=ctx['primitive'];du=ctx['dprimitive'];w=ctx['weight']
   metrics[l]=sum(u[:,c].T@(w[:,c,None]*u[:,c]) for c in range(2))
   dmetrics[l]=sum(du[:,c].T@(w[:,c,None]*u[:,c])+u[:,c].T@(w[:,c,None]*du[:,c]) for c in range(2))
   metrics[l][0,0]=1;metrics[l][0,1]=metrics[l][1,0]=0
   dmetrics[l][0,0]=dmetrics[l][0,1]=dmetrics[l][1,0]=0
   for lp,nlp in counts.items():
    up=self.context[ia,lp]['primitive'];dup=self.context[ia,lp]['dprimitive'];values=[];dvalues=[]
    for vq in v.T:
     values.append(sum(u[:,c].T@((w[:,c]*vq)[:,None]*up[:,c]) for c in range(2)))
     dvalues.append(sum(du[:,c].T@((w[:,c]*vq)[:,None]*up[:,c])+u[:,c].T@((w[:,c]*vq)[:,None]*dup[:,c]) for c in range(2)))
    pairs[l,lp]=np.array(values).transpose(1,2,0);dpairs[l,lp]=np.array(dvalues).transpose(1,2,0)
  p=np.zeros((layout.size,5),complex);dp=p.copy()
  if projector is not None:
   ctx=self.context[ia,2];u=ctx['primitive'];du=ctx['dprimitive'];w=ctx['wp']
   q=sum(u[:,c].T@(w[:,c]*projector[:,c]) for c in range(2));dq=sum(du[:,c].T@(w[:,c]*projector[:,c]) for c in range(2))
   p=layout.projector(2,q,np.eye(5));dp=layout.projector(2,dq,np.eye(5))
  return layout,layout.potential(pairs),layout.overlap(metrics),p,layout.potential(dpairs),layout.overlap(dmetrics),dp
