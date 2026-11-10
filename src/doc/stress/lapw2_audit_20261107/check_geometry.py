"""Extract the actual lapw2 G2 loops; compare with independent sphere quadrature."""
from pathlib import Path
import ctypes as ct, subprocess, json, re, os
import numpy as np
from scipy.special import sph_harm_y

root=Path(__file__).resolve().parent
src=Path(os.environ.get('WIEN_LAPW2_SOURCE',str(root/'SRC_lapw2')))
if not src.exists():src=Path.home()/'wien2k_24/SRC_lapw2'
text=(src/'val_corr_stress.F').read_text()
rad=text.split('             radial_geom = zeroc',1)[1].split('             ! The radial matrix',1)[0]
ang=text.split('             angular_geom = zeroc',1)[1].split('             DO num = nemin,nemax',1)[0]
guard=text.split('       duplicate_radial_channel = .FALSE.',1)[1].split('       IF (duplicate_radial_channel) CYCLE',1)[0]
header='''subroutine probe(l,mbra,lp,mket,lpot,mpot,radial_geom,angular_geom)
implicit none
integer,intent(in)::l,mbra,lp,mket,lpot,mpot
complex*16,intent(out)::radial_geom(3,3),angular_geom(3,3)
integer,parameter::lmax2n=8
complex*16,parameter::zeroc=(0d0,0d0)
real*8,parameter::hf=0.5d0
integer::ih,sp,tp,t,u,v1,v,alpha,beta,lint,nu
integer::lv(1),lpv(1),llv(1),mmv(1)
complex*16::c_am(3,-1:1),ca_st_lm(3,-1:1,-1:1),radial_factor,angular_factor
complex*16::vpuu(1),vpdd(1),vpud(1),vpdu(1),wpuu(1),wpdd(1),wpud(1),wpdu(1)
logical::duplicate_radial_channel
real*8::gaunt1,gaunt_pot,gaunt_couple,channel_diff,channel_scale
real*8,external::gaunt
integer,external::notri
ih=1;lv=lp;lpv=l;llv=lpot;mmv=mpot
vpuu=1;vpdd=1;vpud=1;vpdu=1;wpuu=1;wpdd=1;wpud=1;wpdu=1
call c_alpha_m(c_am)
call c_factor(l,mbra,ca_st_lm)
radial_geom=zeroc;angular_geom=zeroc
'''
code=header+guard+rad+ang+'end subroutine\n'
(root/'geometry_probe.f90').write_text(code)
subprocess.run(['/opt/homebrew/bin/gfortran','-O0','-g','-fcheck=all',
    '-ffree-form','-ffree-line-length-none','-shared','-fPIC',str(root/'geometry_probe.f90'),
    *[str(src/f) for f in ['c_alpha_m.f','gaunt.f','t3j.f','t3j0.f','notri.f','ph.f']],
    '-o',str(root/'geometry_probe.dylib')],check=True)
lib=ct.CDLL(str(root/'geometry_probe.dylib'))
x,w=np.polynomial.legendre.leggauss(40)
x,phi=np.meshgrid(x,np.arange(96)*2*np.pi/96,indexing='ij')
w=np.repeat(w,96)*2*np.pi/96
xyz=np.stack([np.sqrt(1-x*x)*np.cos(phi),np.sqrt(1-x*x)*np.sin(phi),x],axis=-1).reshape(-1,3)
def y(l,m,r):
    return sph_harm_y(l,m,np.arccos(np.clip(r[:,2]/np.linalg.norm(r,axis=1),-1,1)),np.arctan2(r[:,1],r[:,0]))
rows=[]
for ch in [(1,1,1,-1,2,0),(1,1,1,0,2,1),(2,1,2,-1,2,0),(4,0,4,0,2,0),(8,0,8,0,2,0),(8,1,8,1,2,0)]:
    args=[np.array([n],np.int32) for n in ch]
    rad=np.zeros((3,3),complex,order='F');ang=rad.copy(order='F')
    lib.probe_(*[a.ctypes.data_as(ct.c_void_p) for a in args+[rad,ang]])
    lb,mb,lk,mk,lp,mp=ch
    bra=y(lb,mb,xyz);ket=y(lk,mk,xyz);pot=y(lp,mp,xyz)
    rr=np.zeros((3,3),complex);aa=rr.copy()
    for a in range(3):
        for b in range(3):
            strain=np.zeros((3,3));strain[a,b]+=.5;strain[b,a]+=.5
            step=1e-6
            dy=(y(lb,mb,xyz@(np.eye(3)+step*strain))-y(lb,mb,xyz@(np.eye(3)-step*strain)))/(2*step)
            rr[a,b]=np.sum(w*bra.conj()*xyz[:,a]*xyz[:,b]*pot*ket)
            aa[a,b]=np.sum(w*dy.conj()*pot*ket)
    rows.append(dict(channel=ch,radial_error=float(np.max(abs(rad-rr))),
                     angular_error=float(np.max(abs(ang-aa))),angular_trace=float(abs(np.trace(ang)))))
    assert rows[-1]['radial_error']<1e-9,rows[-1]
    assert rows[-1]['angular_error']<1e-8,rows[-1]

# Independently integrate two regular free radial channels.  Check the exact
# indexing expression extracted from production, not a retyped correction.
q,wt=np.polynomial.legendre.leggauss(200);r=(q+1)*1.7/2;wt=wt*1.7/2
k=np.array([1.1,2.3]);f=np.sin(k[:,None]*r)/r
df=k[:,None]*np.cos(k[:,None]*r)/r-np.sin(k[:,None]*r)/r**2
table=np.einsum('ir,jr,r->ij',f,df,wt*r**3)
coef=np.array([.7+.2j,-.4+.1j]);energy=k*k;ev=.31
expected=2*np.real(np.einsum('i,j,ir,jr,j,r->',coef.conj(),coef,df,f,energy-ev,wt*r**3))
expr=re.search(r'\(elo\(lp,jlop\)-e\(num\)\)\*jr3upu_22\(([^\n]*)\)',text).group(1)
indices=expr.split(',')
assert len(indices)==4
actual=0j
for i in range(2):
    for j in range(2):
        mapping={'jlo':i,'jlop':j,'l':0,'lp':0}
        ij=[mapping[s.strip()] for s in indices]
        actual+=coef[i].conjugate()*coef[j]*(energy[j]-ev)*table[ij[0],ij[1]]
error=abs(2*actual.real-expected)
assert error<1e-12,error
report=dict(angular=rows,lo_radial=dict(reference=float(expected),production=float(2*actual.real),error=float(error)),bounds_check='passed at ih=1 with one-based arrays')
(root/'geometry_results.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
