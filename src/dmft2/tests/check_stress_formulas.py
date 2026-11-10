"""Regression checks using production Fortran and independent functional oracles."""
import ctypes as ct
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.special import sph_harm_y

import os
HERE=Path(os.environ['EDMFT_STRESS_TEST_DIR'])
lib=ct.CDLL(str(HERE/'audit_routines.dylib'))
lib.init_audit_fact_()
F=lambda x: np.asfortranarray(x)
I=lambda x: np.asfortranarray(x,dtype=np.int32)
Z=lambda shape: np.zeros(shape,complex,order='F')
def call(name,*args):
    getattr(lib,name+'_')(*[a.ctypes.data_as(ct.c_void_p) for a in args])
def herm(a):return (a+a.conj().T)/2
def invsqrt(a):
    d,v=eigh(a);return (v/np.sqrt(d))@v.conj().T
rng=np.random.default_rng(61037)
def rand(*shape):return rng.normal(size=shape)+1j*rng.normal(size=shape)
report={}

# Verify matrix-square-root derivative against independently diagonalized endpoints.
a=rand(3,3);s=a.conj().T@a+np.eye(3);ds=herm(rand(3,3));x=invsqrt(s)
dxrel=Z((3,3,1,1));si=I(np.ones((3,3,1)))
call('computenormalizationderivative',dxrel,F(ds[:,:,None,None]),F(x[:,:,None]),
     si,I([3]),I([3]),I([1]),I([0]),I([1]))
h=1e-5
dxfd=(invsqrt(s+h*ds)-invsqrt(s-h*ds))/(2*h)
dx=x@dxrel[:,:,0,0]
report['normalization_sylvester']={'absolute_error':float(np.max(abs(dx-dxfd))),
 'hermiticity_error':float(np.max(abs(dx-dx.conj().T)))}

# Adjoint identity with an independently computed Frechet derivative.
m=rand(3,3);qbar=Z((3,3,1))
call('buildnormalizationadjoint',qbar,F(m[:,:,None]),F(x[:,:,None]),si,I([3]),I([3]),I([1]))
lhs=2*np.real(np.sum(m@np.eye(3)*np.linalg.solve(x,dxfd)))
rhs=2*np.real(np.sum(qbar[:,:,0]*ds))
report['normalization_adjoint']={'direct_fd':float(lhs),'adjoint':float(rhs),'error':float(abs(lhs-rhs))}

# Exact finite-dimensional fixed-X projector example, no degeneracies.
nb,nw,nc=6,3,2
energies=np.array([-1.4,-.9,-.25,.4,1.1,1.8])
u=rand(nb,nc)*.3;sig=herm(rand(nc,nc))*.35
dH=np.zeros((nb,nb),complex);dH[nw:,:nw]=rand(nb-nw,nw)*.3;dH[:nw,nw:]=dH[nw:,:nw].conj().T
eff=np.diag(energies[:nw])+u[:nw]@sig@u[:nw].conj().T
ev,v=eigh(eff);temp=.23;chem=-.5;ff=1/(1+np.exp((ev-chem)/temp))
rho=(v*ff)@v.conj().T;rhoe=eff@rho;c=sig@u[:nw].conj().T@rho
rot=dH[nw:,:nw]/(energies[:nw][None,:]-energies[nw:,None])
du=rot.conj().T@u[nw:]
direct=2*np.real(np.trace(c@du))

def window_code(u,rho,rhoe,dh):
    rr=Z((nb,nb));er=Z((nb,nb));grad=Z((nw,nc,1))
    call('buildwindowresponsematrices',rr,er,grad,F(u[:,:,None]),F(rho),F(rhoe),F(energies),
         *[I([q]) for q in [1,nb,1,nw,nb,nw,nc,1,nc,1]],
         I([1]),I([nc]),I([nc]),I(np.arange(1,nc+1)[:,None]),I(np.ones((nc,nc,1))),F((sig@u[:nw].conj().T@rho).T[:,:,None]),I([1]))
    return float(np.trace(rr@dh).real),float(np.max(abs(grad[:,:,0].T-(sig@u[:nw].conj().T@rho)))),rr,er
code,cerr,rr,er=window_code(u,rho,rhoe,dH)

def endpoint(t):
    en,vv=eigh(np.diag(energies)+t*dH)
    # Set each near-identity KS state to positive overlap with the reference.
    vv=vv*np.exp(-1j*np.angle(np.diag(vv)))[None,:]
    pp=vv[:,:nw].conj().T@u
    en2=eigh(np.diag(en[:nw])+pp@sig@pp.conj().T,eigvals_only=True)
    return -temp*np.logaddexp(0,-(en2-chem)/temp).sum(),pp
fd=(endpoint(1e-5)[0]-endpoint(-1e-5)[0])/2e-5
dufd=(endpoint(1e-5)[1]-endpoint(-1e-5)[1])/2e-5
phase=np.exp(1j*np.array([.23,-.8,1.2,.6,-1.5,.9]));q=np.diag(phase);qw=q[:nw,:nw]
codephase,_,_,_=window_code(q.conj().T@u,qw.conj().T@rho@qw,qw.conj().T@rhoe@qw,q.conj().T@dH@q)
corrected_rr=rr.conj()
report['fixed_X_window']={'functional_fd':float(fd),'correct_projector_derivative':float(direct),
 'production_matrix_pair':code,'production_after_band_rephase':codephase,
 'legacy_transposed_matrix_pair':float(np.trace(corrected_rr@dH).real),
 'projector_derivative_fd_error':float(np.max(abs(du-dufd))),
 'C_reconstruction_error':cerr}

# The same outside-window orientation in the global-normalization branch.
qb=herm(rand(nc,nc));rrn=Z((nb,nb));ern=Z((nb,nb));gn=Z((nw,nc,1));wg=.37
call('buildnormalizationwindowresponsematrices',rrn,ern,gn,F(u[:,:,None]),
     F(qb.conj()[:,:,None]),F([wg]),F(energies),
     *[I([q]) for q in [1,nb,1,nw,nb,nw,nc,1,nc,1]],
     I([1]),I([nc]),I([nc]),I(np.arange(1,nc+1)[:,None]))
ds=wg*(du.conj().T@u[:nw]+u[:nw].conj().T@du)
normal_direct=2*np.real(np.trace(qb@ds))
normal_fd=2*np.real(np.trace(qb@(wg*(endpoint(1e-5)[1].conj().T@endpoint(1e-5)[1]-endpoint(-1e-5)[1].conj().T@endpoint(-1e-5)[1])/2e-5)))
report['normalization_window']={'direct_derivative':float(normal_direct),'finite_difference':float(normal_fd),
 'gradient_contraction':float(2*np.real(np.sum(gn[:,:,0]*du))),
 'production_matrix_pair':float(np.trace(rrn@dH).real),
 'legacy_transposed_matrix_pair':float(np.trace(rrn.conj()@dH).real)}

# Direct angular quadrature, independent of the Gaunt summation in the code.
cth,ww=np.polynomial.legendre.leggauss(32);ph=np.arange(64)*2*np.pi/64
cth,ph=np.meshgrid(cth,ph,indexing='ij');ww=np.repeat(ww,64)*2*np.pi/64
vec=np.stack([np.sqrt(1-cth**2)*np.cos(ph),np.sqrt(1-cth**2)*np.sin(ph),cth],axis=-1).reshape((-1,3))
def ylm(l,m,vec):
    theta=np.arccos(np.clip(vec[:,2]/np.linalg.norm(vec,axis=1),-1,1))
    phi=np.arctan2(vec[:,1],vec[:,0]);return sph_harm_y(l,m,theta,phi)
angular_rows=[]
quadrature_rows=[]
for lb,mb,lk,mk,lp,mp in [(1,1,1,0,2,1),(1,1,1,-1,2,0),(2,1,2,-1,2,0),
                           (2,0,2,0,2,0),(4,0,4,0,2,0),(8,0,8,0,2,0),(8,1,8,1,2,0),(0,0,0,0,2,0)]:
    rad=Z((3,3));ang=Z((3,3))
    call('angular_probe',*[I([q]) for q in [lb,mb,lk,mk,lp,mp]],rad,ang)
    radref=Z((3,3));angref=Z((3,3))
    yb=ylm(lb,mb,vec);yk=ylm(lk,mk,vec);yp=ylm(lp,mp,vec)
    for aa in range(3):
        for bb in range(3):
            gen=np.zeros((3,3));gen[aa,bb]+=.5;gen[bb,aa]+=.5
            eps=1e-6
            dy=(ylm(lb,mb,vec@(np.eye(3)+eps*gen))-ylm(lb,mb,vec@(np.eye(3)-eps*gen)))/(2*eps)
            radref[aa,bb]=np.sum(ww*yb.conj()*vec[:,aa]*vec[:,bb]*yp*yk)
            angref[aa,bb]=np.sum(ww*dy.conj()*yp*yk)
    radquad=Z((3,3));angquad=Z((3,3))
    call('quadrature_probe',*[I([q]) for q in [lb,mb,lk,mk,lp,mp]],radquad,angquad)
    quadrature_rows.append({'channel':[lb,mb,lk,mk,lp,mp],'radial_max_error':float(np.max(abs(radquad-radref))),'angular_max_error':float(np.max(abs(angquad-angref)))})
    angular_rows.append({'channel':[lb,mb,lk,mk,lp,mp],
      'radial_max_error':float(np.max(abs(rad-radref))),
      'angular_max_error':float(np.max(abs(ang-angref))),
      'radial_xy':[str(rad[0,1]),str(radref[0,1])],
      'radial_yz':[str(rad[1,2]),str(radref[1,2])],
      'angular_trace':[str(np.trace(ang)),str(np.trace(angref))]})
report['nonspherical_geometry']=angular_rows
report['complete_G2_angular_integrals']=quadrature_rows
# General nonsingular basis change: stationarity cancels all coefficient terms.
t=rand(nw,nw)*.1
hprime=t.conj().T@np.diag(energies[:nw])+np.diag(energies[:nw])@t
oprime=t.conj().T+t
uprime=t.conj().T@u[:nw]
oneparticle=np.trace(rho@hprime-rhoe@oprime).real
projector=2*np.trace(c@uprime).real
report['basis_connection_cancellation']={'oneparticle':float(oneparticle),'projector':float(projector),'sum':float(oneparticle+projector)}
# Normalization in a nonorthogonal band frame: the metric is indispensable.
praw=u[:nw]
ds_naive=uprime.conj().T@praw+praw.conj().T@uprime
ds_metric=ds_naive-praw.conj().T@oprime@praw
report['normalization_metric_identity']={'naive_derivative_norm':float(np.linalg.norm(ds_naive)),'metric_complete_derivative_norm':float(np.linalg.norm(ds_metric))}

# Generalized Hermitian eigenproblem: validate both H' and O' blocks.
dO=np.zeros((nb,nb),complex);dO[nw:,:nw]=rand(nb-nw,nw)*.1;dO[:nw,nw:]=dO[nw:,:nw].conj().T
rg=(dH[nw:,:nw]-dO[nw:,:nw]*energies[:nw][None,:])/(energies[:nw][None,:]-energies[nw:,None])
dug=rg.conj().T@u[nw:]
def generalized_endpoint(t):
    en,vv=eigh(np.diag(energies)+t*dH,np.eye(nb)+t*dO)
    pp=vv[:,:nw].conj().T@u
    vals=eigh(np.diag(en[:nw])+pp@sig@pp.conj().T,eigvals_only=True)
    return -temp*np.logaddexp(0,-(vals-chem)/temp).sum()
report['generalized_window_H_and_O']={'functional_fd':float((generalized_endpoint(1e-5)-generalized_endpoint(-1e-5))/2e-5),'direct_projector':float(2*np.trace(c@dug).real),'production_H_minus_O':float(np.trace(rr@dH-er@dO).real),'legacy_transposed_H_minus_O':float(np.trace(rr.conj()@dH-er.conj()@dO).real)}
# Default normalization response must include the within-window metric,
# not just the complementary-band rotation.  Differentiate the selected
# generalized-eigenvector span at independent endpoints, including O'_WW.
metric_change=herm(np.arange(nb*nb).reshape(nb,nb)*.0007+1j*dH*.21)
metric_change[:nw,:nw]+=np.diag([.07,-.03,.05])
def normalization_endpoint(t):
    _,vv=eigh(np.diag(energies)+t*dH,np.eye(nb)+t*metric_change)
    pw=vv[:,:nw].conj().T@u
    return 2*wg*np.trace(qb@pw.conj().T@pw).real
norm_metric_fd=(normalization_endpoint(1e-5)-normalization_endpoint(-1e-5))/2e-5
norm_metric_actual=np.trace(rrn@dH-ern@metric_change).real
report['normalization_window_metric']={
    'functional_fd':float(norm_metric_fd),
    'production_H_minus_O':float(norm_metric_actual),
    'absolute_error':float(abs(norm_metric_actual-norm_metric_fd))}
assert abs(norm_metric_actual-norm_metric_fd)<1e-9,report['normalization_window_metric']
print(json.dumps(report,indent=2))
(HERE/'results.json').write_text(json.dumps(report,indent=2)+'\n')

# These are independent physical/algebraic oracles, not snapshots of the code.
assert report['normalization_sylvester']['absolute_error']<1e-9
assert report['normalization_adjoint']['error']<1e-9
w=report['fixed_X_window']
assert abs(w['production_matrix_pair']-w['functional_fd'])<1e-9,w
assert abs(w['production_after_band_rephase']-w['functional_fd'])<1e-9,w
assert w['C_reconstruction_error']<1e-10,w
w=report['normalization_window']
assert abs(w['production_matrix_pair']-w['finite_difference'])<1e-9,w
w=report['generalized_window_H_and_O']
assert abs(w['production_H_minus_O']-w['functional_fd'])<1e-9,w
for table in ['nonspherical_geometry','complete_G2_angular_integrals']:
 for row in report[table]:
  assert row['radial_max_error']<1e-9,row
  assert row['angular_max_error']<1e-8,row
assert abs(report['basis_connection_cancellation']['sum'])<1e-12
assert report['normalization_metric_identity']['metric_complete_derivative_norm']<1e-12
print('PASS: complex window, generalized H/O, phase covariance, normalization, and G2 angular identities')

# The production direct-C accumulators must differentiate the log determinant,
# including the ordered advanced/negative-frequency factor for complex Sigma.
uc=rand(4,2)*.2;su=herm(rand(2,2))*.2-1j*np.eye(2)*.11
h0=np.diag([-1.2,-.3,.4,1.1]);z=.15+.31j;uuprime=rand(4,2)*.1
g=np.linalg.inv(z*np.eye(4)-h0-uc@su@uc.conj().T)
for routine,kind in [('accumulatedirectprojectorpair','matsubara'),('accumulatedirectprojectorspectral','real_axis')]:
 cp=Z((2,4,1));work=[Z((2,4)) for _ in range(4)]
 call(routine,cp,F(g),F(uc[:,:,None]),F(su[:,:,None]),I([2]),I([4]),I([2]),I([1]),*work)
 def functional(t):
  v=uc+t*uuprime
  ld=np.log(np.linalg.det(z*np.eye(4)-h0-v@su@v.conj().T))
  return -2*ld.real if kind=='matsubara' else ld.imag/np.pi
 fd=(functional(1e-5)-functional(-1e-5))/2e-5
 predicted=2*np.trace(cp[:,:,0]@uuprime).real
 assert abs(fd-predicted)<1e-9,(kind,fd,predicted)
 print('PASS: direct C',kind,'log-determinant derivative',fd,predicted)
