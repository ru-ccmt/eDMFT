from pathlib import Path
import sys,os,subprocess,shutil,json,re,time
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path[:0]=['/private/tmp/edmft_augmentation_response_20261108','/Users/haule/projects/eDMFT_new/src/python']
import material_augmentation as material
material.R=ROOT
sys.path.insert(0,str(ROOT/'stage'))
from lapw_density_response import MuffinTinDensity,InterstitialDensity,FourierStars,herm
from frozen_sigma_response import FrozenSigmaResponse
from wien_response_density import Density,FLOAT

case=sys.argv[1];central=ROOT/case/'central'
model,records=FrozenSigmaResponse.from_export(central/'response')
mat=material.Material(case,central);s=mat.s
rng=np.random.default_rng(821);rho=[]
for k in range(len(records)):
    rho.append(herm(rng.normal(size=(model.eps.shape[1],)*2)+1j*rng.normal(size=(model.eps.shape[1],)*2))*.01)
rho=np.array(rho)
drho=herm(rng.normal(size=rho[0].shape)+1j*rng.normal(size=rho[0].shape))*.01
drho=np.array([(k+1)*drho for k in range(len(records))])
reference=Density(ROOT/case/'signed'/(case+'.clmval'),central/(case+'.struct'))
joint=reference.zeros();basis=reference.zeros();fixed=reference.zeros();y=reference.from_vector(rng.normal(size=len(reference.vector())))
dv={};gp={};gm={};rhobar=np.zeros_like(rho)
for ia in range(s.nat):
    n=s.jrj[ia];r=s.r0[ia]*np.exp(np.log(s.rmt[ia]/s.r0[ia])/(n-1)*np.arange(n))
    dv[ia]=.07*r*np.sin(1.3*r)
    for l in range(9):
        ctx=mat.radial(ia,l,dv[ia]);gp[ia,l]=np.zeros_like(ctx['primitive']);gm[ia,l]=np.zeros_like(ctx['lo'])
small=1/137.0359895**2 if 'RELA' in s.mode else 4e-22
latom=0
for ia in range(s.nat):
    hs=[key[1:] for key in reference.rad if key[0]==ia]
    phi={l:mat.context[ia,l]['primitive'][:,:2] for l in range(9)}
    dp={l:mat.context[ia,l]['dprimitive'][:,:2] for l in range(9)}
    ym=np.array([y.rad[(ia,*h)] for h in hs])/s.mult[ia]
    for mu in range(s.mult[ia]):
        for k,rho0 in enumerate(rho):
            z,dz,*_=mat.coefficients(ia,latom,k,lmax=8)
            adapter=MuffinTinDensity(z,rho0,phi,hs,small)
            a=adapter.action(drho[k]);b=adapter.action(np.zeros_like(rho0),dz,dp)
            for h,av,bv in zip(hs,a,b):
                fixed.rad[(ia,*h)]+=av/s.mult[ia]
                basis.rad[(ia,*h)]+=bv/s.mult[ia]
                joint.rad[(ia,*h)]+=(av+bv)/s.mult[ia]
            rb,zb,pb=adapter.pullback(ym);rhobar[k]+=rb
            for l in range(9):
                ctx=mat.context[ia,l];gp[ia,l][:,:2]+=pb[l]
                zbl=zb[adapter.layout.slices[l]].reshape(2*l+1,phi[l].shape[2],len(rho0))
                channel=mat.channels[ia,latom,k,l]
                bb,mb=channel.pullback(ctx['primitive'][-1],zbl,lapw=ctx['lapw'])
                gp[ia,l][-1]+=bb;gm[ia,l]+=mb
        latom+=1
    print(case,'moving atom',ia,'assembled',flush=True)
reps=np.array(list(reference.pw));allg=[r.labels[:,:r.nbasis-r.nlo].T for r in records]
extent=np.max([np.ptp(g,axis=0) for g in allg],axis=0)
star_extent=np.max(abs(np.einsum('sij,qj->sqi',s.timat.transpose(0,2,1),reps)),axis=(0,1))
shape=tuple(2*np.maximum(extent,star_extent)+1)
stars=FourierStars(shape,reps,s.timat.transpose(0,2,1),s.tau,complex_phases=mat.complex)
ir=np.zeros(len(reps),complex)
for k,r in enumerate(records):
    pw=r.nbasis-r.nlo;adapter=InterstitialDensity(r.coefficients[:pw,:,0],allg[k],s.a*s.b*s.c,stars)
    ir+=adapter.value(drho[k]);rhobar[k]+=adapter.pullback(np.array(list(y.pw.values())))
for g,v in zip(reference.pw,ir):joint.pw[g]=v;fixed.pw[g]=v
reverse=float(np.vdot(rhobar,drho).real)
for ia in range(s.nat):
    for l in range(9):
        inp=ROOT/case/f'radial_a{ia}_l{l}.input';cov=inp.with_suffix('.covector');out=inp.with_suffix('.pullback')
        with cov.open('wb') as f:
            gp[ia,l].ravel(order='F').tofile(f);gm[ia,l].ravel(order='F').tofile(f)
        subprocess.run([str(ROOT/'radial_pullback'),str(inp),str(cov),str(out)],check=True)
        vb=np.fromfile(out,'f8')[:s.jrj[ia]];reverse+=float(vb@dv[ia])
forward=float(y.vector()@joint.vector())
report=dict(material=case,transpose=dict(forward=forward,reverse=reverse,
    relative=abs(forward-reverse)/max(1.,abs(forward),abs(reverse))),
    direct_basis_norm=float(np.linalg.norm(basis.vector())),
    density_matrix_norm=float(np.linalg.norm(fixed.vector())),steps=[])
work=ROOT/case/'moving';work.mkdir(exist_ok=True)
env=os.environ.copy()
for key in list(env):
    if key.startswith(('EDMFT_','DFTKS_','WIEN_LAPW1')):del env[key]
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
for h in [.002,.001,.0005]:
    ends=[]
    for sign in [-1,1]:
        d=work/(('m_' if sign<0 else 'p_')+f'{h:.0e}')
        if not d.exists():
            d.mkdir()
            for p in central.iterdir():
                if p.is_file() and (p.name.startswith(case+'.') or p.name in ['EF.dat','projectorw.dat','SOlapm.dat','Eorb_imp.dat','dmft2.def'] or p.name.startswith('sig.inp')):
                    shutil.copy2(p,d/p.name)
            p=d/(case+'.vsp');text=p.read_text();patches=[]
            for block in re.finditer(r'ATOMNUMBER\s*=\s*(\d+)',text):
                ia=int(block[1])-1;end=text.find('ATOMNUMBER',block.end());end=len(text) if end<0 else end
                channel=re.search(r'VLM\(R\) FOR L\s*0\s+M=\s*0',text[block.end():end]);offset=block.end()+channel.end()
                for i,match in enumerate(list(FLOAT.finditer(text[offset:end]))[:s.jrj[ia]]):
                    value=float(match.group().replace('D','E'))+sign*h*dv[ia][i]
                    patches.append((offset+match.start(),offset+match.end(),value))
            for a,b,value in reversed(patches):text=text[:a]+f'{value:{b-a}.12E}'+text[b:]
            p.write_text(text)
            model.write_density(d/'rho',rho+sign*h*drho,records)
            with (d/'run.log').open('w') as log:
                run=subprocess.run(['/Users/haule/dbin/dmft2_response','dmft2.def'],cwd=d,
                    env=env|{'EDMFT_RESPONSE_DENSITY_PREFIX':str(d/'rho')},stdout=log,stderr=subprocess.STDOUT)
            run.check_returncode()
            if 'DMFT2 END' not in (d/'run.log').read_text():raise RuntimeError('native endpoint failed')
        ends.append(Density(d/(case+'.clmval'),d/(case+'.struct')))
    fd=(ends[1].vector()-ends[0].vector())/(2*h)
    row=dict(step=h,relative=float(np.linalg.norm(fd-joint.vector())/np.linalg.norm(fd)),maximum=float(np.max(abs(fd-joint.vector()))),
             frozen_basis_relative=float(np.linalg.norm(fd-fixed.vector())/np.linalg.norm(fd)))
    report['steps'].append(row);print(case,row,flush=True)
    (work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
joint.write(work/(case+'.tangent.clmval'));basis.write(work/(case+'.basis.clmval'))
print(case,'transpose',report['transpose'],flush=True)
