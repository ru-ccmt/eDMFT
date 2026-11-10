from pathlib import Path
import numpy as np,re,json,ctypes,subprocess
r=Path(__file__).resolve().parent
# Cartesian angular quadrature, independent of the Fortran Gaunt contractions.
z,w=np.polynomial.legendre.leggauss(24);phi=np.arange(64)*2*np.pi/64
n=np.stack([np.sqrt(1-z[:,None]**2)*np.cos(phi),np.sqrt(1-z[:,None]**2)*np.sin(phi),np.broadcast_to(z[:,None],(24,64))],axis=-1).reshape(-1,3);weights=np.repeat(w,64)*2*np.pi/64
Q=[np.zeros((3,3)),np.diag([-1,-1,2])*np.sqrt(5/(16*np.pi))]
for mode in [2,3,4,5]:
 q=np.zeros((3,3))
 if mode in [2,4]:a=0 if mode==2 else 1;q[a,2]=q[2,a]=np.sqrt(15/(4*np.pi))/2
 if mode==3:q=np.diag([1,-1,0])*np.sqrt(15/(16*np.pi))
 if mode==5:q[0,1]=q[1,0]=np.sqrt(15/(4*np.pi))/2
 Q.append(q)
res={}
for version in ['before','fixed']:
 chunks=(r/f'elect_fixture/output_{version}.txt').read_text().split('CASE')[1:];mx=0
 for ch in chunks:
  side,mode=map(int,ch.splitlines()[0].split());ys=np.einsum('ni,ij,nj->n',n,Q[mode],n);fac=1/3 if mode==0 else 1/(7 if side==3 else 5);field=np.ones(len(n)) if mode==0 else ys**(2 if side==3 else 1)
  ref=-.5*fac*np.einsum('n,na,nb->ab',weights*field,n,n)
  rows=re.findall(r'^:DBG_STR_DENS:\s*(.*)',ch,re.M);obs=np.zeros((3,3))
  for row in rows:
   t=row.split();obs[int(t[0])-1,int(t[1])-1]=float(t[-1])
  err=float(np.max(abs(obs-ref)));mx=max(mx,err)
  if err>1e-6:print('ELECT',version,side,mode,err)
 res[version+'_max_error']=mx
assert res['fixed_max_error']<2e-8,res
res['status']='PASS';(r/'elect_fixture/results.json').write_text(json.dumps(res,indent=2)+'\n');print(res)
# Differentiate the native PBEX/PBEC energy, not a reimplementation of Q.
d=r/'pbe_fixture';d.mkdir(exist_ok=True);src=r/'SRC_lapw0'
(d/'stubs.f90').write_text('module xcparam\nreal(8)::xc1=-1000,xc2=-1000,xc3=-1000\ninteger::icount=0\nend module\n')
code=(src/'vxclm2.F').read_text();(d/'vsigma.f90').write_text(code[code.index('      subroutine pbe_vsigma_builtin'):])
(d/'energy.f90').write_text('''subroutine test_energy(nu,nd,gu,gd,energy,q)
 implicit none
 real(8)::nu,nd,gu(3),gd(3),energy,q(3,3),u,v,g,exu,exd,vxu,vxd,ec,vcu,vcd,su,sd,sud
 logical::cu,cd
 u=sqrt(sum(gu*gu));v=sqrt(sum(gd*gd));g=sqrt(sum((gu+gd)**2))
 call pbex(nu,u,0d0,0d0,nd,v,0d0,0d0,exu,exd,vxu,vxd,1d-18,cu,cd)
 call pbec(nu,u,0d0,nd,v,0d0,g,0d0,ec,vcu,vcd,cu,cd)
 energy=2d0*(nu*exu+nd*exd+(nu+nd)*ec)
 call pbe_vsigma_builtin(nu,nd,u,v,g,.2195149727645171d0,.804d0,.06672455060314922d0,1d0,2,cu,cd,su,sd,sud)
 call gga_metric_tensor(2,gu,gd,su,sd,sud,q)
end subroutine
''')
subprocess.run(['/opt/homebrew/bin/gfortran','-dynamiclib','-O2','-ffree-form','-ffree-line-length-none','-I.','stubs.f90','vsigma.f90','energy.f90']+[str(src/(f+'.f')) for f in ['pbex','pbec','gcor','xcpot','xcener','gga_metric_tensor']]+['-o','pbe.dylib'],cwd=d,check=True)
lib=ctypes.CDLL(str(d/'pbe.dylib'));ptr=ctypes.POINTER(ctypes.c_double)
def calc(nu,nd,gu,gd):
 vals=[ctypes.c_double(x) for x in [nu,nd]];gu=np.asarray(gu,float);gd=np.asarray(gd,float);e=ctypes.c_double();q=np.zeros((3,3),order='F')
 lib.test_energy_(*[ctypes.byref(v) for v in vals],gu.ctypes.data_as(ptr),gd.ctypes.data_as(ptr),ctypes.byref(e),q.ctypes.data_as(ptr));return e.value,q
rng=np.random.default_rng(317);err=0
for k in range(40):
 nu,nd=10**rng.uniform(-2,1,2);gu=rng.normal(size=3)*nu;gd=rng.normal(size=3)*nd;e,q=calc(nu,nd,gu,gd)
 for a in range(3):
  for b in range(a,3):
   direction=np.zeros((3,3));direction[a,b]=direction[b,a]=1 if a==b else .5
   def en(t):
    f=np.eye(3)+t*direction;return calc(nu,nd,np.linalg.solve(f,gu),np.linalg.solve(f,gd))[0]
   h=1e-3;fd=(en(-2*h)-8*en(-h)+8*en(h)-en(2*h))/(12*h);err=max(err,abs(fd+q[a,b]))
assert err<1e-8,err
res={'native_PBE_metric_FD_max_error':err,'samples':40,'status':'PASS'};(d/'results.json').write_text(json.dumps(res,indent=2)+'\n');print(res)
