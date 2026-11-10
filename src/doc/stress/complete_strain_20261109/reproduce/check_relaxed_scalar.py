"""Independent nonlinear lattice scalar at relaxed +/- halfstep endpoints."""
from pathlib import Path
import os,sys,shutil,json,subprocess,re,time
import numpy as np
from scipy.linalg import eigvalsh,lu_factor
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from frozen_sigma_response import FrozenSigmaResponse,herm
from frozen_sigma_tangent import CachedFrozenSigmaResponse

def logdet(a):
 return np.array([np.log(abs(np.diag(lu_factor(x,check_finite=False)[0]))).sum() for x in a])
def scalar(m,mu):
 T=m.meta.temperature;total=0.
 for k,u in enumerate(m.uc):
  eps=m.eps[k];ev=eigvalsh(herm(np.diag(eps-mu)+u@m.sinf@u.conj().T));value=-T*np.logaddexp(0,-ev/T).sum()
  for j in range(0,len(m.meta.omega),64):
   sl=slice(j,j+64);z=mu+1j*m.meta.omega[sl]
   ov=np.einsum('ni,wn,nj->wij',u.conj(),1/(z[:,None]-eps),u,optimize=True)
   value+=2*T*np.dot(m.meta.weights[sl],logdet(np.eye(u.shape[1])-m.sinf@ov)-logdet(np.eye(u.shape[1])-m.sigma[sl]@ov))
  total+=2*m.kw[k]*value
 return float(total)
env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','WIEN_','DFTKS_'))};env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
center,_=FrozenSigmaResponse.from_export(root/'FeSe/refined/response');report={};t=time.perf_counter()
for label in ['minus','plus']:
 src=root/'FeSe/halfstep'/label;out=root/'FeSe/scalar_endpoints'/label;out.mkdir(parents=True,exist_ok=True)
 if not (out/'response.meta').exists():
  for p in src.iterdir():
   if p.is_file() and (p.name.startswith('FeSe.') or p.name in ['EF.dat','Eorb_imp.dat','projectorw.dat','dmft2.def'] or p.name.startswith('sig.inp')):shutil.copy2(p,out/p.name)
  p=out/'FeSe.indmf2';lines=p.read_text().splitlines();i=next(i for i,l in enumerate(lines) if l.strip() and not l.lstrip().startswith('#'));w=lines[i].split();w[1]=(src/'EF.dat').read_text().split()[0];w[2]='0';lines[i]=' '.join(w);p.write_text('\n'.join(lines)+'\n')
  with (out/'export.log').open('w') as f:subprocess.run(['/Users/haule/dbin/dmft2_response','dmft2.def'],cwd=out,env=env|{'EDMFT_RESPONSE_EXPORT':str(out/'response')},stdout=f,stderr=subprocess.STDOUT,check=True)
  if 'DMFT2 END' not in (out/'export.log').read_text():raise ValueError('failed endpoint export')
 m,rec=CachedFrozenSigmaResponse.from_export(out/'response')
 for name in ['temperature','omega','weights','sigma_inf','sigma','dimensions','sym_real','sym_imag']:np.testing.assert_allclose(getattr(m.meta,name),getattr(center.meta,name),atol=1e-13,rtol=0)
 rho=m.central_density();N=float(np.trace(rho,axis1=1,axis2=2).real.sum())
 L=scalar(m,m.meta.mu);vsum=float(re.findall(r'^:VSUM.*?=\s*(\S+)',(out/'FeSe.scf2').read_text(),re.M)[-1]);hist=json.loads((src/'history.json').read_text())[-1]
 report[label]=dict(mu=m.meta.mu,N=N,grand_lattice_Ry=L,legendre_N=(vsum-L)/m.meta.mu,VSUM_Ry=vsum,previous_VSUM_Ry=hist['VSUM_Ry'],K_Ry=hist['K_Ry'])
 print(label,report[label],flush=True)
hm=json.loads((root/'FeSe/halfstep/minus/prepared.json').read_text())['actual'];hp=json.loads((root/'FeSe/halfstep/plus/prepared.json').read_text())['actual'];d=hp-hm
report['derivatives']=dict(dmu=(report['plus']['mu']-report['minus']['mu'])/d,grand_lattice=(report['plus']['grand_lattice_Ry']-report['minus']['grand_lattice_Ry'])/d,native_VSUM=(report['plus']['VSUM_Ry']-report['minus']['VSUM_Ry'])/d)
report['seconds']=time.perf_counter()-t;(root/'FeSe/scalar_endpoints/REPORT.json').write_text(json.dumps(report,indent=2)+'\n');print(report['derivatives'],flush=True)
