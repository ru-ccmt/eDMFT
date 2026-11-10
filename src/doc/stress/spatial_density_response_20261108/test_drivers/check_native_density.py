from pathlib import Path
import sys,os,subprocess,shutil,json,time
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path[:0]=['/private/tmp/edmft_augmentation_response_20261108','/Users/haule/projects/eDMFT_new/src/python']
import material_augmentation as material
material.R=ROOT
sys.path.insert(0,str(ROOT/'stage'))
from lapw_density_response import MuffinTinDensity,InterstitialDensity,FourierStars,herm
from frozen_sigma_response import FrozenSigmaResponse
from wien_response_density import Density

if not (ROOT/'export_radial').exists():shutil.copy2('/private/tmp/edmft_augmentation_response_20261108/export_radial',ROOT/'export_radial')
case=sys.argv[1];central=ROOT/case/'central'
model,records=FrozenSigmaResponse.from_export(central/'response')
mat=material.Material(case,central);s=mat.s
rng=np.random.default_rng(821);rho=[]
for k in range(len(records)):
    r=herm(rng.normal(size=(model.eps.shape[1],)*2)+1j*rng.normal(size=(model.eps.shape[1],)*2))*.01
    rho.append(r)
rho=np.array(rho)
model.write_density(central/'signed_rho',rho,records)
dest=ROOT/case/'signed'
if not dest.exists():
    dest.mkdir()
    for p in central.iterdir():
        if p.is_file() and (p.name.startswith(case+'.') or p.name in ['EF.dat','projectorw.dat','SOlapm.dat','Eorb_imp.dat','dmft2.def'] or p.name.startswith('sig.inp')):
            shutil.copy2(p,dest/p.name)
    env=os.environ.copy()
    for key in list(env):
        if key.startswith(('EDMFT_','DFTKS_','WIEN_LAPW1')):del env[key]
    env.update(EDMFT_RESPONSE_DENSITY_PREFIX=str(central/'signed_rho'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    with (dest/'run.log').open('w') as log:p=subprocess.run(['/Users/haule/dbin/dmft2_response','dmft2.def'],cwd=dest,env=env,stdout=log,stderr=subprocess.STDOUT)
    p.check_returncode()
    if 'DMFT2 END' not in (dest/'run.log').read_text():raise RuntimeError('density reference failed')
reference=Density(dest/(case+'.clmval'),dest/(case+'.struct'))
out=reference.zeros()
for ia in range(s.nat):
    for l in range(9):mat.radial(ia,l)
print(case,'radials prepared',flush=True)
small=1/137.0359895**2 if 'RELA' in s.mode else 4e-22
latom=0
for ia in range(s.nat):
    hs=[key[1:] for key in reference.rad if key[0]==ia]
    phi={l:mat.context[ia,l]['primitive'][:,:2] for l in range(9)}
    for mu in range(s.mult[ia]):
        for k,r in enumerate(rho):
            # The exact central eigenvector gauge has to agree with signed input.
            np.testing.assert_allclose(mat.cv[k][1],records[k].coefficients[:,:,0],atol=1e-13,rtol=0)
            z,*_=mat.coefficients(ia,latom,k,lmax=8)
            adapter=MuffinTinDensity(z,r,phi,hs,small)
            values=adapter.value()/s.mult[ia]
            for h,v in zip(hs,values):out.rad[(ia,*h)]+=v
        latom+=1
    print(case,'atom',ia,'assembled',flush=True)

reps=np.array(list(reference.pw));allg=[r.labels[:,:r.nbasis-r.nlo].T for r in records]
extent=np.max([np.ptp(g,axis=0) for g in allg],axis=0)
star_extent=np.max(abs(np.einsum('sij,qj->sqi',s.timat.transpose(0,2,1),reps)),axis=(0,1))
shape=tuple(2*np.maximum(extent,star_extent)+1)
stars=FourierStars(shape,reps,s.timat.transpose(0,2,1),s.tau,complex_phases=mat.complex)
ir=np.zeros(len(reps),complex)
for k,r in enumerate(records):
    pw=r.nbasis-r.nlo
    adapter=InterstitialDensity(r.coefficients[:pw,:,0],allg[k],s.a*s.b*s.c,stars)
    ir+=adapter.value(rho[k])
for g,v in zip(reference.pw,ir):out.pw[g]=v

rows=[]
for key,a in reference.rad.items():
    b=out.rad[key]
    rows.append(dict(channel=key,relative=float(np.linalg.norm(b-a)/max(1e-30,np.linalg.norm(a))),maximum=float(np.max(abs(b-a))),scale=float(np.max(abs(a)))))
a=np.array(list(reference.pw.values()));b=np.array(list(out.pw.values()))
report=dict(material=case,mt_channels=rows,interstitial=dict(relative=float(np.linalg.norm(b-a)/np.linalg.norm(a)),maximum=float(np.max(abs(b-a))),shape=list(map(int,shape))))
(ROOT/case/'native_density.json').write_text(json.dumps(report,indent=2)+'\n')
out.write(ROOT/case/(case+'.assembled.clmval'))
print(json.dumps(dict(material=case,mt_maximum=max(x['maximum'] for x in rows),mt_max_relative=max(x['relative'] for x in rows),interstitial=report['interstitial']),indent=2),flush=True)
