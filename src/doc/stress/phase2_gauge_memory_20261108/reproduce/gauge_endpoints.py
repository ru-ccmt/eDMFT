"""Native constant-gauge endpoints; isolated copies, no source/material edits."""
from pathlib import Path
import sys, os, re, shutil, subprocess, json
import numpy as np
from scipy.io import FortranFile
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_native_snapshot import read_basis_matrix
from wien_response_density import Density
case=sys.argv[1];central=Path('/private/tmp/edmft_full_k_response_20261108')/case/'charge/0'
num=re.compile(r'[-+]?\d+\.\d+[DE][-+]\d+')
with load_native_valence_state(central,case,complex_vectors=case=='SrVO3',energy_policy='fixed',matrix_quadrature='wien2k24',library='/private/tmp/edmft_phase2_gauge_memory_20261108/build/liblapw_response.so') as f:
    state=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
    select,_=state.hartree.native_potential_output(tuple(f.bands.spherical.slices),
        {a:op.layout.harmonics for a,(op,ks) in f.bands.ns.items()},f.bands.ir[0].representatives,complex_interstitial=case=='SrVO3')
    field=np.zeros(state.hartree.size);pw=np.zeros(state.hartree.stars.nstar,complex)
    pw[np.flatnonzero(np.all(state.hartree.stars.representatives==0,axis=1))[0]]=1.
    field[state.hartree.pw_offset:]=state.hartree.reanalysis.action(pw).view(float)
    _,_,theta=f.bands.input.unpack(select.matvec(field))
    nr=f.model.eps.shape[1];gauge=[np.zeros((len(b.retained),nr),complex) for b in f.bases]
    radial={a:{l:c.action(c.radius,1.,np.ones(c.nlo)) for l,c in cs.items()} for a,cs in f.contexts.items()}
    for a,l,k,op in f.bands.spherical.spherical:
        c=f.contexts[a][l]
        dh,do=op.action(*radial[a][l],energy_dot=1.,lo_energies_dot=np.ones(c.nlo));gauge[k]+=dh
    for a,op,ks in f.bands.spherical.nonspherical:
        for k,dh in zip(ks,op.action({l:radial[a][l][0] for l in op.phi},{l:radial[a][l][1] for l in op.phi},retained=nr)):gauge[k]+=dh
    for k,ir in enumerate(f.bands.ir):
        a,b=ir.action(theta,retained=nr);gauge[k]+=np.concatenate((a,b))
    analytic=np.concatenate([a.ravel() for a in gauge]);ident=np.concatenate([np.vstack((np.eye(nr),np.zeros((b.solver.nq,nr)))).ravel() for b in f.bases])
    conversion=float((central/'EF.dat').read_text().split()[0])/f.model.meta.mu
    env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))}
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    exe=Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if case=='SrVO3' else 'lapw1')
    ends={};report=dict(case=case,analytic_Hprime_minus_identity=float(np.linalg.norm(analytic-ident)),steps={})
    for h in [.125,-.125]:
        folder=root/case/str(h);folder.mkdir(parents=True,exist_ok=True)
        if not (folder/'done').exists():
            for p in central.iterdir():
                if p.is_file() and (p.name.startswith(case+'.') or p.name in ['EF.dat','projectorw.dat','Eorb_imp.dat','dmft2.def','lapw1.def'] or p.name.startswith('sig.inp')):
                    if p.suffix in ['.vector','.clmval','.output1','.output2'] or p.name.startswith(case+'.scf'):continue
                    shutil.copy2(p,folder/p.name)
            # VSP contains r*V_spherical in Ry, not Y00 coefficients.
            p=folder/(case+'.vsp');blocks=re.split(r'(ATOMNUMBER\s*=\s*\d+)',p.read_text())
            for a in range(f.structure.nat):
                i=2+2*a;vals=num.findall(blocks[i]);assert len(vals)==f.structure.jrj[a]
                r=f.contexts[a][0].radius;shifted=iter([float(x.replace('D','E'))+h*y for x,y in zip(vals,r)])
                blocks[i]=num.sub(lambda m:f'{next(shifted):19.12E}',blocks[i])
            p.write_text(''.join(blocks))
            p=folder/(case+'.vns');lines=p.read_text().splitlines();j=next(i for i,l in enumerate(lines) if 'NUMBER OF PW' in l);assert int(lines[j].split()[0])==len(theta)
            for i,t in enumerate(theta):
                l=lines[j+1+i];values=[float(x.replace('D','E')) for x in num.findall(l)]
                value=complex(*values)+h*t;lines[j+1+i]=l[:18]+f'{value.real:19.12E}{value.imag:19.12E}'
            p.write_text('\n'.join(lines)+'\n')
            p=folder/(case+('.in1c' if case=='SrVO3' else '.in1'));lines=p.read_text().splitlines()
            lines[0]=re.sub(r'(EF=)([-+\d.eE]+)',lambda m:m[1]+f'{float(m[2])+h:.16g}',lines[0])
            for i,l in enumerate(lines[2:],2):
                if l.startswith('K-VECTORS'):
                    lines[i]=re.sub(r'(UNIT:\s*\d+\s+)([-+\d.eE]+)',lambda m:m[1]+f'{float(m[2])+h:.16g}',l)
                elif 'CONT' in l or 'STOP' in l:
                    lines[i]=re.sub(r'^(\s*\d+\s+)([-+\d.eE]+)',lambda m:m[1]+f'{float(m[2])+h:.6f}',l)
                elif l.strip():
                    lines[i]=re.sub(r'^(\s*)([-+\d.eE]+)',lambda m:m[1]+f'{float(m[2])+h:.16g}',l)
            p.write_text('\n'.join(lines)+'\n')
            p=folder/(case+'.indmf2');lines=p.read_text().splitlines();j=next(i for i,l in enumerate(lines) if l.strip() and not l.startswith('#'));fields=lines[j].split();fields[1]=f'{(f.model.meta.mu+h)*conversion:.16g}';fields[2]='0';lines[j]=' '.join(fields);p.write_text('\n'.join(lines)+'\n');(folder/'EF.dat').write_text(fields[1]+'\n')
            for name in ['SOlapm.dat','SOlapm.datup','SOlapm.datdn']:(folder/name).unlink(missing_ok=True)
            for label,cmd,extra in [('lapw1',[str(exe),'lapw1.def'],{'WIEN_LAPW1_MATRIX_DUMP':str(folder/'native')}),('dmft2',['/Users/haule/dbin/dmft2_response','dmft2.def'],{'EDMFT_RESPONSE_EXPORT':str(folder/'response')})]:
                with (folder/(label+'.log')).open('w') as log:subprocess.run(cmd,cwd=folder,env=env|extra,stdout=log,stderr=subprocess.STDOUT,check=True)
                if label=='dmft2' and 'DMFT2 END' not in (folder/(label+'.log')).read_text():raise RuntimeError(str(folder))
            (folder/'done').touch()
        hs=[];os_=[];rerr=[]
        for k,b in enumerate(f.bases):
            c=np.column_stack((b.retained,b.complement))
            hh=read_basis_matrix(folder/f'native.ham.k{k+1:08d}.bin',b'EDMFTHV1');oo=read_basis_matrix(folder/f'native.k{k+1:08d}.bin',b'EDMFTOV1')
            h0=read_basis_matrix(central/f'native.ham.k{k+1:08d}.bin',b'EDMFTHV1');o0=read_basis_matrix(central/f'native.k{k+1:08d}.bin',b'EDMFTOV1')
            assert np.array_equal(hh.labels,h0.labels)
            hs.append(c.conj().T@((hh.matrix-h0.matrix)@c[:,:nr])/h)
            os_.append(c.conj().T@((oo.matrix-o0.matrix)@c[:,:nr])/h)
        hv=np.concatenate([a.ravel() for a in hs]);ov=np.concatenate([a.ravel() for a in os_])
        density=Density(folder/(case+'.clmval'),folder/(case+'.struct')).vector()
        ends[h]=density
        row=dict(native_Hprime_minus_identity=float(np.linalg.norm(hv-ident)),native_minus_analytic_Hprime=float(np.linalg.norm(hv-analytic)),native_Oprime=float(np.linalg.norm(ov)),density_change_norm=float(np.linalg.norm(density-f.reference.vector())),density_change_relative=float(np.linalg.norm(density-f.reference.vector())/np.linalg.norm(f.reference.vector())))
        report['steps'][h]=row;print(case,h,row,flush=True)
    report['density_central_derivative_norm']=float(np.linalg.norm((ends[.125]-ends[-.125])/.25))
    (root/(case+'_gauge_native.json')).write_text(json.dumps(report,indent=2)+'\n')
