"""Independent LAPW0 total-potential chain-rule checks on archived full meshes.

Only scratch LAPW0 runs are launched: no SCF loop, native KS solve or material
directory change. Reuse previous explicit hydrostatic source endpoints and q.
"""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,time
import numpy as np
root=Path(__file__).resolve().parent
if root.name=='reproduce':root=root.parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_potential_core_state import load_potential_core_state
from lapw_strain_potential import NativePotentialFileLayout,read_potential_file,unequal_centered_derivative
from wien_response_density import Density
from wstruct import Struct

def run(case):
    previous=Path('/private/tmp/edmft_screened_native_20261108/full_bz')/case
    center=previous/'central'
    work=root/case;work.mkdir(exist_ok=True)
    pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
    lengths=pc.structure.jrj
    layout=NativePotentialFileLayout(pc.hartree,
        read_potential_file(center/(case+'.vsp'),lengths),
        read_potential_file(center/(case+'.vns'),lengths),complex_interstitial=(case=='SrVO3'))
    q=np.load(previous/'analytic_density.npy')
    induced=layout.selector.matvec(pc.potential_operator.matvec(q))
    covectors=np.load(Path('/private/tmp/edmft_projector_functional_20261108/full_bz')/case/'covectors.npz')
    g=layout.covector(covectors['potential'])
    n0=Density(center/(case+'.clmsum'),center/(case+'.struct'))
    np.testing.assert_array_equal(n0.vector(),pc.input_density.vector())
    assert len(q)==len(n0.vector())
    s0=Struct();s0.ReadStruct(str(center/(case+'.struct')))
    cell0=np.array([s0.a,s0.b,s0.c])
    env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))}
    env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    exe=Path('/private/tmp/edmft_xc_core_response_20261108/native/SRC_lapw0/lapw0')
    report=dict(case=case,scope='native file-coordinate explicit/induced potential split; not total stress',
                fixed_rmt=True,radial_energy='fixed',mesh={'SrVO3':75,'FeSe':140}[case],steps={},timings=[])
    for step in [.00025,.000125]:
        explicit_values=[];total_values=[];deltas=[];cells=[]
        for label,sign in [('m',-1),('p',1)]:
            source=previous/f'source_{label}_{step}'
            s=Struct();s.ReadStruct(str(source/(case+'.struct')))
            for attr in ['r0','rmt','jrj','pos','mult','rotloc','timat','tau']:
                np.testing.assert_array_equal(getattr(s,attr),getattr(s0,attr))
            cell=np.array([s.a,s.b,s.c]);scale=float(np.mean(cell/cell0));eta=scale-1
            assert np.sign(eta)==sign
            deltas.append(eta);cells.append(cell)
            explicit_values.append(layout.read(source/(case+'.vsp'),source/(case+'.vns')))
            folder=work/f'total_{label}_{step}';folder.mkdir(exist_ok=True)
            for name in [case+'.struct',case+'.in0',case+'.inm','lapw0.def']:
                shutil.copy2(source/name,folder/name)
            n0.copy().add(n0.from_vector(q),eta).transport(scale).write(folder/(case+'.clmsum'))
            started=time.monotonic()
            with (folder/'lapw0.log').open('w') as log:
                subprocess.run([str(exe),'lapw0.def'],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
            text=(folder/'lapw0.log').read_text()
            if any(x in text for x in ['Fortran runtime error','ERROR STOP','Error termination']):
                raise RuntimeError(text[-1000:])
            if (folder/'lapw0.error').read_text().strip():
                raise RuntimeError((folder/'lapw0.error').read_text())
            report['timings'].append(dict(folder=folder.name,seconds=time.monotonic()-started))
            total_values.append(layout.read(folder/(case+'.vsp'),folder/(case+'.vns')))
        explicit=unequal_centered_derivative(explicit_values[0],layout.central,explicit_values[1],*deltas)
        numeric=unequal_centered_derivative(total_values[0],layout.central,total_values[1],*deltas)
        analytic=explicit+induced
        row=dict(actual_steps=deltas,strain_direction=((cells[1]-cells[0])/(cell0*(deltas[1]-deltas[0]))).tolist(),
                 total_relative=float(np.linalg.norm(analytic-numeric)/np.linalg.norm(numeric)),
                 induced_relative=float(np.linalg.norm(analytic-numeric)/np.linalg.norm(induced)),sectors={})
        for name,sl in [('spherical',slice(0,layout.coordinates.spherical_size)),
                        ('nonspherical',slice(layout.coordinates.spherical_size,layout.coordinates.pw_offset)),
                        ('interstitial',slice(layout.coordinates.pw_offset,None))]:
            row['sectors'][name]=dict(relative=float(np.linalg.norm(analytic[sl]-numeric[sl])/np.linalg.norm(numeric[sl])),
                explicit_norm=float(np.linalg.norm(explicit[sl])),induced_norm=float(np.linalg.norm(induced[sl])),
                total_norm=float(np.linalg.norm(numeric[sl])))
        row['projector_potential_only_Ry']=dict(explicit=float(g@explicit),induced=float(g@induced),
            total=float(g@analytic),native_total=float(g@numeric),difference=float(g@(analytic-numeric)))
        report['steps'][str(step)]=row
        np.savez(work/f'potential_split_{step}.npz',explicit=explicit,induced=induced,total=analytic,native=numeric,gradient=g)
        print(case,step,json.dumps(row),flush=True)
    x,y=[report['steps'][str(h)] for h in [.00025,.000125]]
    report['error_reduction']=x['total_relative']/y['total_relative']
    report['explicit_step_relative']=float(np.linalg.norm(np.load(work/'potential_split_0.00025.npz')['explicit']-
        np.load(work/'potential_split_0.000125.npz')['explicit'])/np.linalg.norm(np.load(work/'potential_split_0.000125.npz')['explicit']))
    report['inputs']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [center/(case+'.struct'),center/(case+'.clmsum'),exe]}
    (work/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(case,'DONE',report['error_reduction'],flush=True)

if __name__=='__main__':run(sys.argv[1])
