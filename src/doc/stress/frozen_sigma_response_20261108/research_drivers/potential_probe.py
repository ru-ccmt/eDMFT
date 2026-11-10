"""Measure a local-potential operator in the fixed central LAPW basis.

This diagnostic uses the existing WIEN_LAPW1_BASIS_VSP switch.  LAPW1 is
only a matrix-element evaluator: radial augmentation and cell stay central.
No charge loop is run.  Projection uses the dual of the COMPLETE central
LAPW eigensystem, so omitted endpoint eigenvectors affect the central
derivative only at higher order.  The central metric and H are checked.
"""
import os,sys,re,shutil,time,json
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from scipy.special import spherical_jn
from setup import ROOT,SOURCES,ENV
from strain_source import command,WIEN
from wien_response_density import Density
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
sys.path.insert(0,'/private/tmp/projector_cross_material_20261107')
from analyze_cases import vectors
from check_interstitial_overlap_stress import read_atoms
from wstruct import Struct

def align_rows(v,labels,d):
    npw=d.nbasis-d.nlo
    lookup={tuple(x):i for i,x in enumerate(labels[:,:npw].T)}
    wanted=[tuple(x) for x in d.labels[:,:npw].T]
    if len(v)!=d.nbasis or set(lookup)!=set(wanted):raise ValueError('plane-wave basis changed')
    order=[lookup[x] for x in wanted]+list(range(npw,d.nbasis))
    return v[order]

class PotentialProbe:
    def __init__(self,case,dumps):
        self.case=case;self.src=SOURCES[case];self.root=ROOT/case;self.dumps=dumps
        self.structure=self.src/(case+'.struct');self.ss=Struct();self.ss.ReadStruct(str(self.structure))
        self.cell=np.array([self.ss.a,self.ss.b,self.ss.c]);self.vol=float(np.prod(self.cell));self.atoms=read_atoms(self.structure)
        self.rho=Density(self.src/(case+'.clmsum'),self.structure)
        full=Path('/private/tmp/projector_actual_stress_20261107/complete_center/SrVO3.vector_x') if case=='SrVO3' else Path('/private/tmp/projector_cross_material_20261107/FeSe/complete_center/FeSe.vector')
        cache=self.root/'dual.npz'
        if cache.exists():
            a=np.load(cache);self.duals=[a[f'd{k}'] for k in range(len(dumps))]
        else:
            self.duals=[];err=0.
            for d,(e,v,labels,rad) in zip(dumps,vectors(full,self.ss.nat)):
                if len(e)!=len(v):raise ValueError('incomplete central dual basis')
                v=align_rows(v,labels,d)
                inv=solve(v,np.eye(len(v)),check_finite=False)
                if case=='SrVO3':
                    gauge=np.load('/private/tmp/projector_actual_stress_20261107/dftks_wien_gauge.npy')[len(self.duals)]
                    err=max(err,float(np.max(abs(gauge@gauge.conj().T-np.eye(d.nstored)))))
                    self.duals.append(gauge@inv[:d.nstored])
                else:
                    t=inv@d.coefficients[:,:,0];err=max(err,float(np.max(abs(t.conj().T@t-np.eye(d.nstored)))))
                    self.duals.append(t.conj().T@inv)
            if err>1e-8:raise ValueError(('central metric',err))
            np.savez_compressed(cache,**{f'd{k}':v for k,v in enumerate(self.duals)})
        self.base=self.prepare('potential_base');self.rho.write(self.base/(case+'.clmsum'))
        self.run_lapw0(self.base);self.h0=self.measure(self.base)
        self.h0error=float(np.max(abs(self.h0-np.array([np.diag(d.energies) for d in dumps]))))
        text=(self.base/(case+'.output0')).read_text();self.shape=tuple(map(int,re.search(r'global FFT-parameters:\s*(\d+)\s+(\d+)\s+(\d+)',text).groups()))
        self.xyz=np.array(np.meshgrid(*[np.fft.fftfreq(n)*n for n in self.shape],indexing='ij'))
        self.qnorm=np.sqrt(np.sum((self.xyz*(2*np.pi/self.cell)[:,None,None,None])**2,axis=0))
        self.theta=self.make_theta()
        errs={}
        for trans in [False,True]:
            v,mask=self.stars(self.base/(case+'.vtotal'),trans);w,_=self.stars(self.base/(case+'.vns'),trans)
            errs[trans]=float(np.max(abs((self.product(v,self.theta)-w.conj())[mask])))
        self.transpose=min(errs,key=errs.get)
        if errs[self.transpose]>1e-8:raise ValueError(('Fourier reconstruction',errs))
        self.v0,self.mask=self.stars(self.base/(case+'.vtotal'))
        self.w0,_=self.stars(self.base/(case+'.vns'))

    def prepare(self,name):
        d=self.root/name;d.mkdir(parents=True,exist_ok=True)
        for suf in ['struct','in0','in1','in1c','inc','inm','klist','kgen','rotlm','clmcor']:
            p=self.src/(self.case+'.'+suf)
            if p.exists():shutil.copy2(p,d/p.name)
        p=d/(self.case+'.in0');a=p.read_text().splitlines();a[0]='TOT '+a[0][4:];a[1]=a[1].replace('NR2V','R2V ');p.write_text('\n'.join(a)+'\n')
        df=Path('/private/tmp/projector_actual_stress_20261107/complete_center/lapw1.def') if self.case=='SrVO3' else self.src/'lapw1.def'
        (d/'lapw1.def').write_text(df.read_text().replace('_x',''))
        return d
    def run_lapw0(self,d):return command(d,'lapw0',[str(WIEN/'x'),'lapw0','-f',self.case])
    def measure(self,potential_dir):
        d=self.prepare(potential_dir.name+'_matrix')
        for suf in ['vsp','vns']:shutil.copy2(potential_dir/(self.case+'.'+suf),d/(self.case+'.'+suf))
        command(d,'lapw1',[str(WIEN/('lapw1c' if self.case=='SrVO3' else 'lapw1')),'lapw1.def'],{'WIEN_LAPW1_BASIS_VSP':str(self.src/(self.case+'.vsp'))})
        h=[]
        for dump,dual,(e,v,labels,rad) in zip(self.dumps,self.duals,vectors(d/(self.case+'.vector'),self.ss.nat)):
            v=align_rows(v,labels,dump)
            t=dual@v;h.append((t*e)@t.conj().T)
        if len(h)!=len(self.dumps):raise ValueError('k count changed')
        return np.array(h)
    def stars(self,p,transposed=None):
        if transposed is None:transposed=self.transpose
        lines=p.read_text().splitlines();j=next(i for i,x in enumerate(lines) if 'TOTAL POTENTIAL IN INTERSTITIAL' in x)+2;n=int(lines[j].split()[0]);a=np.zeros(self.shape,complex);mask=np.zeros(self.shape,bool)
        for line in lines[j+1:j+1+n]:
            g=np.array([int(line[i:i+5]) for i in [3,8,13]]);z=float(line[18:37])+1j*float(line[37:56]);star={}
            for rot,tau in zip(self.ss.timat,self.ss.tau):
                gg=tuple(((rot.T if transposed else rot)@g)%self.shape);star.setdefault(gg,[]).append(np.exp(2j*np.pi*np.dot(tau,g)))
            for gg,ph in star.items():
                phase=np.mean(ph);a[gg]=z*(phase if p.suffix=='.vtotal' else phase.conjugate())/(len(star) if p.suffix=='.vtotal' else 1);mask[gg]=True
        return a,mask
    def make_theta(self):
        t=np.zeros(self.shape,complex)
        for atom in self.atoms:
            x=self.qnorm*atom.rmt;f=np.ones_like(x);np.divide(3*spherical_jn(1,x),x,out=f,where=x!=0)
            for p in atom.positions:t-=4*np.pi*atom.rmt**3/(3*self.vol)*f*np.exp(-2j*np.pi*np.einsum('a,aijk->ijk',p,self.xyz))
        t[0,0,0]+=1;return t
    def product(self,a,b):return np.fft.fftn(np.fft.ifftn(a).real*np.fft.ifftn(b).real)*np.prod(self.shape)
    def fourier_matrix(self,op):
        result=[]
        for d in self.dumps:
            npw=d.nbasis-d.nlo;labels=d.labels[:,:npw];diff=(labels[:,:,None]-labels[:,None,:])%np.array(self.shape)[:,None,None];a=d.coefficients[:npw,:,0]
            result.append(a.conj().T@op[tuple(diff)]@a)
        return np.array(result)
    def affine_operator(self,h):
        tag=f'{h:.7f}';plus=self.root/('source_p_'+tag);minus=self.root/('source_m_'+tag)
        hp=json.loads((plus/'source.json').read_text())['actual_step'];hm=json.loads((minus/'source.json').read_text())['actual_step']
        raw=(self.measure(plus)-self.measure(minus))/(hp-hm)
        vp,_=self.stars(plus/(self.case+'.vtotal'));vm,_=self.stars(minus/(self.case+'.vtotal'))
        wp,_=self.stars(plus/(self.case+'.vns'));wm,_=self.stars(minus/(self.case+'.vns'))
        true=self.product((vp-vm)/(hp-hm),self.theta);true[~self.mask]=0
        interstitial_delta=self.fourier_matrix(true-((wp-wm)/(hp-hm)).conj())
        if self.case=='SrVO3':adv=np.load('/private/tmp/projector_actual_stress_20261107/advection_half_potential_derivative.npy')
        else:adv=np.load('/private/tmp/projector_cross_material_20261107/FeSe/matrices.npz')['advection']
        return raw+interstitial_delta+adv
    def apply_kernel(self,density,name,step=1e-5,core=False):
        d=self.prepare(name);self.rho.copy().add(density,step).write(d/(self.case+'.clmsum'))
        self.run_lapw0(d);mat=(self.measure(d)-self.h0)/step
        core_response=None
        if core:
            command(d,'lcore',[str(WIEN/'x'),'lcore','-f',self.case])
            if not (self.base/(self.case+'.clmcor_reference')).exists():
                command(self.base,'lcore',[str(WIEN/'x'),'lcore','-f',self.case])
                shutil.copy2(self.base/(self.case+'.clmcor'),self.base/(self.case+'.clmcor_reference'))
            now=Density(d/(self.case+'.clmcor'),self.structure);old=Density(self.base/(self.case+'.clmcor_reference'),self.structure)
            core_response=self.rho.zeros().add(now).add(old,-1).scale(1/step)
        return mat,core_response
