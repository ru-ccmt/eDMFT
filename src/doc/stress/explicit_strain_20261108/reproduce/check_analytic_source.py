"""Analytic hydrostatic source vs independent native LAPW0 without cell rounding."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,time
import numpy as np
root=Path(__file__).resolve().parent
if root.name=='reproduce':root=root.parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_potential_core_state import load_potential_core_state
from lapw_strain_potential import NativePotentialFileLayout,read_potential_file
from lapw_hydrostatic_potential import hydrostatic_explicit_potential
from wien_response_density import Density

def run(case,total=False):
    previous=Path('/private/tmp/edmft_screened_native_20261108/full_bz')/case
    center=previous/'central';work=root/case/('total_hydro' if total else 'pure_hydro');work.mkdir(exist_ok=True)
    pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
    lengths=pc.structure.jrj
    layout=NativePotentialFileLayout(pc.hartree,read_potential_file(center/(case+'.vsp'),lengths),
        read_potential_file(center/(case+'.vns'),lengths),complex_interstitial=(case=='SrVO3'))
    started=time.monotonic();source,parts=hydrostatic_explicit_potential(pc);seconds=time.monotonic()-started
    selected=layout.selector.matvec(source)
    central=layout.selector.matvec(parts['central_potential'])
    covectors=np.load(Path('/private/tmp/edmft_projector_functional_20261108/full_bz')/case/'covectors.npz')
    g=layout.covector(covectors['potential'])
    q=np.load(previous/'analytic_density.npy');induced=layout.selector.matvec(pc.potential_operator.matvec(q))
    if total:selected=selected+induced
    n0=Density(center/(case+'.clmsum'),center/(case+'.struct'))
    env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))}
    env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    exe=root/'native/SRC_lapw0/lapw0'
    report=dict(case=case,scope='analytic hydrostatic potential plus prescribed density tangent, no total stress' if total else 'analytic hydrostatic potential source, no total stress',
        seconds=seconds,central_relative=float(np.linalg.norm(central-layout.central)/np.linalg.norm(layout.central)),
        projector_potential_Ry=float(g@selected),includes_induced=total,steps={},timings=[])
    for step in [.00025,.000125]:
        endpoints=[]
        for sign in [-1,1]:
            folder=work/str(sign*step);folder.mkdir(exist_ok=True)
            for name in [case+'.struct',case+'.in0',case+'.inm','lapw0.def']:
                shutil.copy2(center/name,folder/name)
            (folder/'lapw0.strain_scale').write_text(f'{1+sign*step:.17g}\n')
            density=n0.copy().add(n0.from_vector(q),sign*step) if total else n0
            density.transport(1+sign*step).write(folder/(case+'.clmsum'))
            t=time.monotonic()
            with (folder/'lapw0.log').open('w') as f:
                subprocess.run([str(exe),'lapw0.def'],cwd=folder,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
            if (folder/'lapw0.error').read_text().strip():raise RuntimeError((folder/'lapw0.error').read_text())
            assert 'DIAGNOSTIC uniform lattice scale:' in (folder/(case+'.output0')).read_text()
            report['timings'].append(dict(endpoint=folder.name,seconds=time.monotonic()-t))
            endpoints.append(layout.read(folder/(case+'.vsp'),folder/(case+'.vns')))
        numeric=(endpoints[1]-endpoints[0])/(2*step)
        row=dict(relative=float(np.linalg.norm(selected-numeric)/np.linalg.norm(numeric)),sectors={},
            scalar_native_Ry=float(g@numeric),scalar_error_Ry=float(g@(selected-numeric)))
        for name,sl in [('spherical',slice(0,layout.coordinates.spherical_size)),
                        ('nonspherical',slice(layout.coordinates.spherical_size,layout.coordinates.pw_offset)),
                        ('interstitial',slice(layout.coordinates.pw_offset,None))]:
            row['sectors'][name]=dict(relative=float(np.linalg.norm((selected-numeric)[sl])/np.linalg.norm(numeric[sl])),
                absolute=float(np.linalg.norm((selected-numeric)[sl])))
        report['steps'][str(step)]=row
        np.savez(work/f'comparison_{step}.npz',analytic=selected,native=numeric,source=source,
                 induced_from_previous_q=induced,gradient=g,**parts)
        print(case,step,json.dumps(row),flush=True)
    a=np.load(work/'comparison_0.00025.npz')['native'];b=np.load(work/'comparison_0.000125.npz')['native']
    extrap=(4*b-a)/3
    report['richardson_relative']=float(np.linalg.norm(selected-extrap)/np.linalg.norm(selected))
    report['richardson_scalar_error_Ry']=float(g@(selected-extrap))
    report['native_binary_sha256']=hashlib.sha256(exe.read_bytes()).hexdigest()
    (work/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(case,'DONE',report['richardson_relative'],flush=True)

if __name__=='__main__':run(sys.argv[1],total='--total' in sys.argv[2:])
