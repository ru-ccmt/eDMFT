"""Recorded nonmagnetic orthogonal-cell material fixture for native checks."""
from pathlib import Path
import sys
import numpy as np
from scipy.special import spherical_jn
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'stage'),'/Users/haule/projects/eDMFT_new/src/python','/private/tmp/edmft_augmentation_response_20261108']
from material_augmentation import Material
sys.path.insert(0,str(ROOT/'stage'))
from lapw_radial_native import NativeRadialContext
from lapw_augmentation_response import PreparedAugmentationChannel
from lapw_spherical_response import PreparedSphericalBand,prepare_wronskian_band
from lapw_nonspherical_response import PreparedNonspherical
from lapw_projector_response import PreparedProjector5
from lapw_spherical_potential_response import SphericalPotentialBands
from lapw_density_response import MuffinTinDensitySum,SpatialValenceDensity,InterstitialDensity,FourierStars
from lapw_basis_density_response import BasisDensityResponse
from frozen_sigma_tangent import CachedFrozenSigmaResponse
from indmffile import Indmfl
from wien_response_density import Density


def build(case,source=None,keep_bands=None):
 central=source or ROOT/case/'central';mat=Material(case,central);s=mat.s
 model,records=CachedFrozenSigmaResponse.from_export(central/'response')
 exported_bands=model.eps.shape[1]
 # The native dmft2 reader has a 100-Ry ceiling. In the complete-basis
 # diagnostic two empty SrVO3 states lie above it. Their eigenvectors and
 # energies are present in case.vector; assemble their raw overlaps below.
 # Provisional zero overlaps outside W do not change the central density.
 extend=keep_bands is None and len(mat.cv[0][0])>exported_bands
 if extend:
  n=len(mat.cv[0][0]);raw_native=np.zeros((len(records),n,records[0].maxdim,records[0].norbitals),complex)
  normalized=np.zeros_like(raw_native)
  for orb,ic in enumerate(model.meta.cix):
   for j in range(model.meta.nindo[orb]):raw_native[:,:exported_bands,j,orb]=model.raw[:,:,model.meta.isx[j,orb],ic]
  for k,r in enumerate(records):
   assert len(mat.cv[k][0])==n and np.min(mat.cv[k][0][exported_bands:])>model.meta.mu+100*model.meta.temperature
   normalized[k,:exported_bands]=r.projector
   r.nstored=n;r.nemax=r.nemin+n-1;r.coefficients=mat.cv[k][1][:,:,None];r.energies=mat.cv[k][0]
  model=CachedFrozenSigmaResponse(model.meta,np.array([cv[0] for cv in mat.cv]),raw_native,normalized,model.kw,model.W)
 if keep_bands is not None:
  from copy import copy
  n=int(keep_bands)
  if n<=max(model.W) or n>model.eps.shape[1]:raise ValueError('invalid retained-band control')
  normalized=np.array([r.projector[:n] for r in records]);raw_native=np.zeros_like(normalized)
  for orb,ic in enumerate(model.meta.cix):
   for j in range(model.meta.nindo[orb]):raw_native[:,:,j,orb]=model.raw[:,:n,model.meta.isx[j,orb],ic]
  model=CachedFrozenSigmaResponse(model.meta,model.eps[:,:n],raw_native,normalized,model.kw,model.W)
  records=[copy(r) for r in records]
  for r in records:
   r.nstored=n;r.nemax=r.nemin+n-1;r.coefficients=r.coefficients[:,:n];r.projector=r.projector[:n];r.energies=r.energies[:n]
  mat.cv=[(ev[:n],c[:,:n],labels) for ev,c,labels in mat.cv]
 rho=model.central_density();nk,nb=model.eps.shape
 reference=Density(central/(case+'.clmval'),central/(case+'.struct'))
 inl=Indmfl(str(central/case));inl.read();ctx={};sphere=[];nonsphere=[];project=[]
 ltop,lns=map(int,(central/(case+('.in1c' if mat.complex else '.in1'))).read_text().splitlines()[1].split()[1:3])
 axes=[];locations=[];first=np.r_[0,np.cumsum(s.mult)[:-1]];site_atom=np.repeat(np.arange(s.nat),s.mult)
 for site,a in enumerate(site_atom):
  for rotation,tau in zip(s.timat,s.tau):
   diff=rotation@s.pos[first[a]]+tau-s.pos[site];diff-=np.rint(diff)
   if max(abs(diff))<2e-7:break
  else:raise RuntimeError('missing site rotation')
  axes.append(rotation.T);locations.append(rotation@s.pos[first[a]]+tau)
 geometry=[]
 for k,(ev,c,labels) in enumerate(mat.cv):
  np.testing.assert_allclose(c,records[k].coefficients[:,:,0],atol=1e-13,rtol=0)
  with (central/f'native.k{k+1:08d}.bin').open('rb') as f:f.seek(44);kf=np.fromfile(f,'f8',3)
  frac=labels.T+kf;q=2*np.pi*frac/np.array([s.a,s.b,s.c]);qn=np.linalg.norm(q,axis=1)
  geometry.append((c,labels,frac,q,qn,len(c)-mat.nlo_total))
 density_atoms={};density_matching={};local={}
 volume=s.a*s.b*s.c
 for a in range(s.nat):
  n=s.jrj[a];ctx[a]={}
  for l in range(ltop+1):
   e=mat.header[a][l];lapw=bool(e<150);e=e if lapw else e-200
   elo=mat.lo[a][l,:mat.ilo[a][l]] if l<4 else np.array([])
   ctx[a][l]=NativeRadialContext(mat.v[a][0,0],r0=s.r0[a],dx=np.log(s.rmt[a]/s.r0[a])/(n-1),rmt=s.rmt[a],energy=e,
       lo_energies=elo,ell=l,charge=s.Znuc[a],relativistic=bool('RELA' in s.mode),lapw=lapw,library=ROOT/'liblapw_response.so')
  dens=[];density_matching[a]=[];ns=[];ks=[]
  for site in range(first[a],first[a]+s.mult[a]):
   for k,(c,labels,frac,q,qn,npw) in enumerate(geometry):
    phase=np.exp(2j*np.pi*(labels.T@s.pos[site]));loc=(s.rotloc[a]@axes[site]@q.T).T
    channels={};z={}
    for l,rc in ctx[a].items():
     groups=[npw+mat.lo_indices[a,l,j] for j in range(rc.nlo)]
     matching=np.array([spherical_jn(l,qn[:npw]*rc.rmt),qn[:npw]*spherical_jn(l,qn[:npw]*rc.rmt,derivative=True)])
     ch=PreparedAugmentationChannel(l,rc.rmt,volume,q,qn,phase,c,npw,groups,matching_values=matching)
     triangle=prepare_wronskian_band(l,rc.rmt,volume,q[:npw],phase[:npw],c[:npw],matching)
     sphere.append((a,l,k,PreparedSphericalBand(ch,rc.primitive,rc.lo_coefficients,rc.overlap_metric,rc.radius,rc.energy,rc.lo_energies,triangle,lapw=rc.lapw)))
     if l<=8:
      ch=PreparedAugmentationChannel(l,rc.rmt,volume,loc,qn,phase,c,npw,groups,matching_values=matching)
      channels[l]=ch;z[l]=ch.value(rc.primitive[-1],rc.lo_coefficients,lapw=rc.lapw).reshape(-1,nb)
    dz=np.concatenate([z[l] for l in range(9)]);dens.append((k,dz,1/s.mult[a]));density_matching[a].append(channels)
    ns.append((np.concatenate([z[l] for l in range(lns+1)]),{l:channels[l] for l in range(lns+1)},
                {l:ctx[a][l].lo_coefficients for l in range(lns+1)},{l:ctx[a][l].lapw for l in range(lns+1)}));ks.append(k)
  hs=[key[1:] for key in reference.rad if key[0]==a]
  small=1/137.0359895**2 if 'RELA' in s.mode else 4e-22
  density_atoms[a]=MuffinTinDensitySum(dens,rho,{l:ctx[a][l].primitive[:,:2] for l in range(9)},hs,small)
  nonsphere.append((a,PreparedNonspherical({l:ctx[a][l].primitive for l in range(lns+1)},ctx[a][0].overlap_metric,
                      list(mat.vns[a]),np.array(list(mat.vns[a].values())).T,ns),ks))
  print(case,'prepared atom',a,flush=True)
 radial_projector={}
 with (central/'projectorw.dat').open() as f:
  for _ in range(int(f.readline().split()[1])):
   nr,nr0,a,l=map(int,f.readline().split()[1:]);rows=np.array([list(map(float,f.readline().split())) for _ in range(nr)])
   if nr!=nr0:raise NotImplementedError('fixture needs explicit projector-tail construction')
   radial_projector[a-1,l]=rows[:,1:3]
 raw=np.zeros_like(model.raw)
 for ic,orbs in inl.cix.items():
  if len(orbs)!=1:raise NotImplementedError('fixture is single-site block only')
  site,l,_=orbs[0];site-=1;a=int(site_atom[site]);rc=ctx[a][l]
  locrot,xyz,shift,rmt2=inl.atoms[site+1]
  if shift is not None:raise NotImplementedError('fixture has no shifted orbital')
  xyz=np.array(xyz) if len(xyz) else s.rotloc[a]
  if xyz.shape!=(3,3):raise NotImplementedError('fixture needs explicit three-axis frame')
  for k,(c,labels,frac,q,qn,npw) in enumerate(geometry):
   phase=np.exp(2j*np.pi*(frac@locations[site]));qloc=(xyz@axes[site]@q.T).T
   groups=[npw+mat.lo_indices[a,l,j] for j in range(rc.nlo)]
   ch=PreparedAugmentationChannel(l,rc.rmt,volume,qloc,qn,phase,c,npw,groups)
   op=PreparedProjector5(ch,rc.primitive,rc.lo_coefficients,rc.projector_metric,radial_projector[a,l],inl.cftrans[ic].T,lapw=rc.lapw)
   cols=slice(0,2*l+1);raw[k,:,cols,ic-1]=op.value();project.append((a,l,k,int(ic)-1,cols,op))
 np.testing.assert_allclose(raw[:,:exported_bands],model.raw[:,:exported_bands],atol=1e-12,rtol=1e-10)
 if extend:
  normalized_local=model.apply_x(raw,model.X);raw_native=np.zeros_like(raw);normalized_native=np.zeros_like(raw)
  for orb,ic in enumerate(model.meta.cix):
   for j in range(model.meta.nindo[orb]):
    raw_native[:,:,j,orb]=raw[:,:,model.meta.isx[j,orb],ic]
    normalized_native[:,:,j,orb]=normalized_local[:,:,model.meta.isx[j,orb],ic]
  model=CachedFrozenSigmaResponse(model.meta,model.eps,raw_native,normalized_native,model.kw,model.W)
  np.testing.assert_allclose(model.central_density(),rho,atol=1e-13,rtol=1e-13)
 reps=np.array(list(reference.pw));allg=[r.labels[:,:r.nbasis-r.nlo].T for r in records]
 extent=np.max([np.ptp(g,axis=0) for g in allg],axis=0)
 star_extent=np.max(abs(np.einsum('sij,qj->sqi',s.timat.transpose(0,2,1),reps)),axis=(0,1))
 shape=tuple(2*np.maximum(extent,star_extent)+1)
 stars=FourierStars(shape,reps,s.timat.transpose(0,2,1),s.tau,complex_phases=mat.complex)
 ir=[InterstitialDensity(r.coefficients[:r.nbasis-r.nlo,:,0],g,volume,stars) for r,g in zip(records,allg)]
 spatial=SpatialValenceDensity(reference,density_atoms,ir)
 direct=BasisDensityResponse(spatial,{a:{l:ctx[a][l] for l in range(9)} for a in ctx},density_matching,energy_policy='fixed')
 bands=SphericalPotentialBands(ctx,sphere,nonsphere,project,rho.shape,raw.shape,energy_policy='fixed')
 return dict(model=model,records=records,mat=mat,contexts=ctx,bands=bands,spatial=spatial,direct=direct,reference=reference)
