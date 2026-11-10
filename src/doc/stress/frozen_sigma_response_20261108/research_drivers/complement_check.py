"""Converge the charge-response complement without changing the DMFT window."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys,json,shutil,time
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from setup import ROOT,SOURCES,EXE
from strain_source import WIEN,command
from test_kernel import load
from frozen_sigma_response import FrozenSigmaResponse
from wien_response_density import Density
from potential_probe import vectors,align_rows
from screen_response import assemble

def run(case,nb):
    old=SOURCES[case];root=ROOT/case/f'nb{nb}';center=root/'center';center.mkdir(parents=True,exist_ok=True)
    # All files are independent copies except vector files produced below.
    for p in (ROOT/case/'central').iterdir():
        if p.is_file() and (p.name.startswith(case+'.') and p.suffix in ['.struct','.in0','.in1','.in1c','.in2','.in2c','.inc','.inm','.klist','.kgen','.rotlm','.vsp','.vns','.clmsum','.clmcor','.indmfl','.indmfi','.indmf2'] or p.name in ['EF.dat','projectorw.dat','Eorb_imp.dat','SOlapm.dat','dmft2.def'] or p.name.startswith('sig.inp')):shutil.copy2(p,center/p.name)
    in1=center/(case+('.in1c' if case=='SrVO3' else '.in1'));a=in1.read_text().splitlines();a[-1]=f'K-VECTORS FROM UNIT:4  -10.0  2.0  {nb}  FIXED_NBANDS';in1.write_text('\n'.join(a)+'\n')
    df=Path('/private/tmp/projector_actual_stress_20261107/complete_center/lapw1.def') if case=='SrVO3' else old/'lapw1.def'
    (center/'lapw1.def').write_text(df.read_text().replace('_x',''));p=center/'dmft2.def';p.write_text(p.read_text().replace('_x',''))
    exe=str(WIEN/('lapw1c' if case=='SrVO3' else 'lapw1'))
    command(center,'lapw1',[exe,'lapw1.def'])
    command(center,'export',['/opt/homebrew/bin/mpirun','-n','4',str(EXE),'dmft2.def'],{'EDMFT_RESPONSE_EXPORT':str(center/'response')})
    f,d=FrozenSigmaResponse.from_export(center/'response');nat=len(Density(center/(case+'.clmsum'),center/(case+'.struct')).grids)
    full=Path('/private/tmp/projector_actual_stress_20261107/complete_center/SrVO3.vector_x') if case=='SrVO3' else Path('/private/tmp/projector_cross_material_20261107/FeSe/complete_center/FeSe.vector')
    dual=[];metric=0.
    for dd,(e,v,l,r) in zip(d,vectors(full,nat)):
        v=align_rows(v,l,dd);inv=solve(v,np.eye(len(v)),check_finite=False);t=inv@dd.coefficients[:,:,0];metric=max(metric,float(np.max(abs(t.conj().T@t-np.eye(nb)))));dual.append(t.conj().T@inv)
    if metric>1e-8:raise ValueError(('metric',metric))
    operators=[]
    for name in ['potential_base','reference_induced']:
        dest=root/name;dest.mkdir(exist_ok=True)
        for fn in [case+'.struct',in1.name,case+'.klist',case+'.kgen',case+'.rotlm','lapw1.def']:shutil.copy2(center/fn,dest/fn)
        for suff in ['vsp','vns']:shutil.copy2(ROOT/case/name/(case+'.'+suff),dest/(case+'.'+suff))
        command(dest,'lapw1',[exe,'lapw1.def'],{'WIEN_LAPW1_BASIS_VSP':str(center/(case+'.vsp'))})
        hs=[]
        for dd,dddual,(e,v,l,r) in zip(d,dual,vectors(dest/(case+'.vector'),nat)):
            v=align_rows(v,l,dd);t=dddual@v;hs.append((t*e)@t.conj().T)
        operators.append(np.array(hs))
    dh=(operators[1]-operators[0])/1e-5
    dr,info=f.response(dh,True);f.write_density(root/'rho',dr,d)
    # Reconstruct against exactly the extended central eigenvector gauge.
    SOURCES[case]=center
    val=assemble(case,f,d,dr,f'nb{nb}_density')
    template=Density(center/(case+'.clmsum'),center/(case+'.struct'))
    source=Density(ROOT/case/'source_0.0002500.density',center/(case+'.struct'))
    ref=Density(ROOT/case/'reference_q.density',center/(case+'.struct'))
    core=Density(ROOT/case/'reference_induced'/(case+'.clmcor'),center/(case+'.struct')).add(Density(ROOT/case/'potential_base'/(case+'.clmcor_reference'),center/(case+'.struct')),-1).scale(1e5)
    pred=template.zeros().add(source).add(val).add(core);error=pred.copy().add(ref,-1)
    result={'bands':nb,'metric_error':metric,'resolvent':info,'reference_norm':float(np.linalg.norm(ref.vector())),'fixed_point_residual_norm':float(np.linalg.norm(error.vector()))}
    pred.write(root/'predicted_q.density');np.save(root/'potential.npy',dh);(root/'results.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':run(sys.argv[1],int(sys.argv[2]))
