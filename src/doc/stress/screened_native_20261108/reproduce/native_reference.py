"""Isolated native hydrostatic sources and independent charge-map actions.

All writes stay below this checkpoint. No user calculation is run or edited.
Fixed G labels, stored bands, radial energies and Sigma are retained. Endpoint
mu is obtained from an independent nonlinear trace, never analytic dmu.
"""
from pathlib import Path
import hashlib, json, os, re, shutil, subprocess, sys, time
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent
sys.path[:0] = ['/Users/haule/projects/eDMFT_new/src/python',
                '/private/tmp/edmft_full_k_response_20261108']
from frozen_sigma_response import FrozenSigmaResponse
from number_reference import NumberReference
from wien_response_density import Density
from wien_frozen_sigma_response import vector_signature
from lapw_charge_response import NativeChargeNormalization, native_mixer_inversion
from wstruct import Struct


class NativeReference:
    def __init__(self, case, mesh='full_bz'):
        self.case = case
        self.mesh = mesh
        old = Path('/private/tmp/edmft_full_k_response_20261108')
        self.source = (old / mesh if mesh == 'full_bz' else old) / case / 'charge/0'
        self.xcroot = Path('/private/tmp/edmft_xc_core_response_20261108')
        self.work = ROOT / mesh / case
        self.work.mkdir(parents=True, exist_ok=True)
        self.structure_path = self.source / (case + '.struct')
        self.structure = Struct(); self.structure.ReadStruct(str(self.structure_path))
        self.rho = Density(self.source/(case+'.clmsum'), self.structure_path)
        self.cell = np.array([self.structure.a, self.structure.b, self.structure.c])
        self.base, _ = FrozenSigmaResponse.from_export(self.source/'response')
        self.target_band = NumberReference(self.base)(self.base.meta.mu)
        self.conversion = float((self.source/'EF.dat').read_text().split()[0])/self.base.meta.mu
        self.window = tuple(map(int,(self.source/(case+'.indmfl')).read_text().splitlines()[0].split()[:2]))
        self.signature = vector_signature(self.source/(case+'.vector'), self.structure.nat, self.window)
        self.inversion = native_mixer_inversion(self.structure)
        core = Density(self.xcroot/case/'hxc_core/0'/(case+'.clmcor'), self.structure_path)
        raw = Density(self.source/(case+'.clmval'), self.structure_path).add(core)
        from lapw_charge_response import native_charge_covectors
        m,i=native_charge_covectors(self.rho,self.structure,inversion=self.inversion)
        self.target = float(round((m+i)@self.rho.vector()))
        self.normalization = NativeChargeNormalization(raw,self.structure,self.target,inversion=self.inversion)
        self.env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))}
        self.env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
        self.timings=[]

    def command(self, folder, label, exe, extra=None):
        start=time.monotonic()
        with (folder/(label+'.log')).open('w') as f:
            subprocess.run([str(exe), label.split('_')[0]+'.def'],cwd=folder,
                           env=self.env|(extra or {}),stdout=f,stderr=subprocess.STDOUT,check=True)
        text=(folder/(label+'.log')).read_text()
        if any(x in text for x in ['Fortran runtime error','Error termination','ERROR STOP']):
            raise RuntimeError(str(folder/(label+'.log')))
        native=label.split('_')[0]
        if native in ['lapw1','dmft2'] and native.upper()+' END' not in text:
            raise RuntimeError(str(folder/(label+'.log')))
        for p in folder.glob(native+'.error*'):
            if p.read_text().strip(): raise RuntimeError(str(p))
        self.timings.append(dict(folder=str(folder),program=label,seconds=time.monotonic()-start))

    def evaluate(self, label, density=None, scale=1.0):
        if density is None: density=self.rho
        folder=self.work/label
        request=dict(density_sha256=hashlib.sha256(density.vector().tobytes()).hexdigest(),scale=scale,
                     target=self.target,target_band=self.target_band)
        if (folder/'done.json').exists():
            old=json.loads((folder/'done.json').read_text())
            if old['request'] != request: raise ValueError('cached endpoint request differs')
            saved=np.load(folder/'normalized.npy')
            if saved.shape != self.rho.vector().shape:
                layout=Density(folder/(self.case+'.clmval'),folder/(self.case+'.struct'))
                saved=self.rho.zeros().add(layout.from_vector(saved)).vector()
                np.save(folder/'normalized.npy',saved)
            return saved,old['actual_scale']
        folder.mkdir(exist_ok=True)
        # Copy only inputs, not prior output, SOlapm, vectors or matrix exports.
        suffixes=['struct','in1','in1c','in2','in2c','indmfl','indmfi','indmf2','klist','kgen','rotlm']
        names=[self.case+'.'+s for s in suffixes]+['EF.dat','Eorb_imp.dat','projectorw.dat','lapw1.def','dmft2.def']
        names += [p.name for p in self.source.glob('sig.inp*')]
        for name in names:
            p=self.source/name
            if p.exists(): shutil.copy2(p,folder/name)
        for name in [self.case+'.in0',self.case+'.inm','lapw0.def']:
            shutil.copy2(self.xcroot/self.case/'xc/0'/name,folder/name)
        for name in [self.case+'.inc','lcore.def']:
            shutil.copy2(self.xcroot/self.case/'hxc_core/0'/name,folder/name)
        lines=self.structure_path.read_text().splitlines()
        lines[3]=''.join(f'{a*scale:10.6f}' for a in self.cell)+lines[3][30:]
        (folder/(self.case+'.struct')).write_text('\n'.join(lines)+'\n')
        actual=np.array([float(lines[3][i:i+10]) for i in [0,10,20]])/self.cell
        if np.ptp(actual)>2e-7: raise ValueError('rounded cell no longer hydrostatic')
        scale=float(actual.mean())
        density.transport(scale).write(folder/(self.case+'.clmsum'))
        self.command(folder,'lapw0',self.xcroot/'native/SRC_lapw0/lapw0')
        if label=='central': (folder/'lapw1.fixgvec.write').touch()
        else:
            shutil.copy2(self.work/'central/lapw1.fixgvec',folder/'lapw1.fixgvec')
            (folder/'lapw1.fixgvec.read').touch()
        exe=Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if self.case=='SrVO3' else 'lapw1')
        self.command(folder,'lapw1',exe)
        sig=vector_signature(folder/(self.case+'.vector'),self.structure.nat,self.window)
        if sig != self.signature: raise ValueError('native G/band/k signature differs')
        self.command(folder,'dmft2','/Users/haule/dbin/dmft2_response',{'EDMFT_RESPONSE_EXPORT':str(folder/'response')})
        endpoint,_=FrozenSigmaResponse.from_export(folder/'response')
        oracle=NumberReference(endpoint)
        mu=brentq(lambda x:oracle(x)-self.target_band,self.base.meta.mu-.2,self.base.meta.mu+.2,xtol=5e-15,rtol=1e-14)
        p=folder/(self.case+'.indmf2'); lines=p.read_text().splitlines()
        j=next(j for j,l in enumerate(lines) if l.strip() and not l.startswith('#'))
        fields=lines[j].split();fields[1]=f'{mu*self.conversion:.16g}';fields[2]='0';lines[j]=' '.join(fields)
        p.write_text('\n'.join(lines)+'\n');(folder/'EF.dat').write_text(f'{mu*self.conversion:.16g}\n')
        self.command(folder,'dmft2_fixed_n','/Users/haule/dbin/dmft2_response')
        self.command(folder,'lcore',self.xcroot/'native/SRC_lcore/lcore')
        s=Struct();s.ReadStruct(str(folder/(self.case+'.struct')))
        raw=Density(folder/(self.case+'.clmval'),folder/(self.case+'.struct'))
        raw.add(Density(folder/(self.case+'.clmcor'),folder/(self.case+'.struct')))
        physical=NativeChargeNormalization(raw,s,self.target,inversion=self.inversion).central
        # Native output may append zero high-G stars at a geometric cutoff.
        # Align to the central labels, rejecting any nonzero discarded mode.
        pulled=self.rho.zeros().add(raw.from_vector(physical)).transport(scale,True)
        # The native reference defines material-coordinate density using a
        # second normalization after pullback; N^2=N at the central output.
        out=NativeChargeNormalization(pulled,self.structure,self.target,inversion=self.inversion).central
        np.save(folder/'normalized.npy',out)
        record=dict(request=request,actual_scale=scale,mu_Ry=mu,number_residual=oracle(mu)-self.target_band)
        (folder/'done.json').write_text(json.dumps(record,indent=2)+'\n')
        print(self.case,label,'native completed',flush=True)
        return out,scale

    def sources(self):
        base,_=self.evaluate('central')
        central_error=float(np.linalg.norm(base-self.normalization.central)/np.linalg.norm(base))
        if central_error>2e-9: raise ValueError(f'central native output differs: {central_error}')
        report=dict(central_relative_error=central_error,target_charge=self.target,target_band_number=self.target_band,
                    scope='hydrostatic one-step source at the archived central state; not final stress',steps={})
        for h in [.00025,.000125]:
            p,lp=self.evaluate('source_p_'+str(h),scale=1+h)
            m,lm=self.evaluate('source_m_'+str(h),scale=1-h)
            b=(p-m)/(lp-lm);np.save(self.work/('source_'+str(h)+'.npy'),b)
            report['steps'][h]=dict(norm=float(np.linalg.norm(b)),actual_difference=lp-lm,
                                  charge=float((self.normalization.mt+self.normalization.ir)@b))
        a=np.load(self.work/'source_0.00025.npy');b=np.load(self.work/'source_0.000125.npy')
        report['step_relative_difference']=float(np.linalg.norm(a-b)/np.linalg.norm(b))
        report['timings']=self.timings
        (self.work/'source_report.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__': NativeReference(sys.argv[1],sys.argv[2] if len(sys.argv)>2 else 'full_bz').sources()
