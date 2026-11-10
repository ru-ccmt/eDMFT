"""Check the tangent of the complete one-electron charge map.

No impurity call and no SCF loop.  Each potential is evaluated once, with
fresh projector normalization and fixed particle number.  Differencing these
single passes includes every LAPW eigenvector and, optionally, radial-basis
response.  This is an independent check on the stored-band analytic tangent.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys,json,shutil,re,time
from pathlib import Path
import numpy as np
from setup import ROOT,SOURCES,EXE
from strain_source import WIEN,command
from wien_response_density import Density
from check_charge import integrate

def evaluate(case,potential,name,frozen=False):
    src=SOURCES[case];d=ROOT/case/name;d.mkdir(exist_ok=True,parents=True)
    for p in (ROOT/case/'central').iterdir():
        if p.is_file() and ((p.name.startswith(case+'.') and p.suffix in ['.struct','.in0','.in1','.in1c','.in2','.in2c','.inc','.inm','.klist','.kgen','.rotlm','.clmsum','.indmfl','.indmfi','.indmf2']) or p.name in ['EF.dat','projectorw.dat','Eorb_imp.dat','dmft2.def'] or re.fullmatch(r'sig\.inp\d*',p.name)):
            shutil.copy2(p,d/p.name)
    for suff in ['vsp','vns']:shutil.copy2(potential/(case+'.'+suff),d/(case+'.'+suff))
    df=Path('/private/tmp/projector_actual_stress_20261107/complete_center/lapw1.def') if case=='SrVO3' else src/'lapw1.def'
    (d/'lapw1.def').write_text(df.read_text().replace('_x',''))
    p=d/'dmft2.def';p.write_text(p.read_text().replace('_x',''))
    p=d/(case+'.indmf2');a=p.read_text().splitlines();i=next(i for i,l in enumerate(a) if l.strip() and not l.startswith('#'));a[i]=f'c {(src/"EF.dat").read_text().split()[0]} 1 1.0 .T. 0.02';p.write_text('\n'.join(a)+'\n')
    if (d/'SOlapm.dat').exists():(d/'SOlapm.dat').unlink()
    times={}
    extra={'WIEN_LAPW1_BASIS_VSP':str(src/(case+'.vsp'))} if frozen else {}
    times['lapw1']=command(d,'lapw1',[str(WIEN/('lapw1c' if case=='SrVO3' else 'lapw1')),'lapw1.def'],extra)
    if frozen:
        # dmft2 must reconstruct the SAME frozen LAPW basis as LAPW1.
        # Merely freezing LAPW1 while letting ATPAR use the trial V would
        # mix coefficient and radial gauges and invalidate this comparison.
        shutil.copy2(d/(case+'.vsp'),d/(case+'.vsp_actual'))
        shutil.copy2(src/(case+'.vsp'),d/(case+'.vsp'))
    times['dmft2']=command(d,'dmft2',['/opt/homebrew/bin/mpirun','-n','4',str(EXE),'dmft2.def'])
    if frozen:shutil.copy2(d/(case+'.vsp_actual'),d/(case+'.vsp'))
    times['lcore']=command(d,'lcore',[str(WIEN/'x'),'lcore','-f',case])
    rho=Density(src/(case+'.clmsum'),src/(case+'.struct'))
    out=rho.zeros().add(Density(d/(case+'.clmval'),d/(case+'.struct'))).add(Density(d/(case+'.clmcor'),d/(case+'.struct')))
    target=round(integrate(rho,src/(case+'.struct'))[2]);mt,ir,_=integrate(out,d/(case+'.struct'))
    for g in out.pw:out.pw[g]*=(target-mt)/ir
    out.write(d/(case+'.clmout'));(d/'timing.json').write_text(json.dumps(times,indent=2))
    return out

def run(case,frozen=False):
    root=ROOT/case;label='frozen' if frozen else 'radial'
    before=evaluate(case,root/'potential_base','full_map_'+label+'_base',frozen)
    after=evaluate(case,root/'reference_induced','full_map_'+label+'_induced',frozen)
    response=after.add(before,-1).scale(1e5)
    source=Density(root/'source_0.0002500.density',SOURCES[case]/(case+'.struct'))
    ref=Density(root/'reference_q.density',SOURCES[case]/(case+'.struct'))
    predicted=source.add(response);err=predicted.copy().add(ref,-1)
    predicted.write(root/('full_map_'+label+'_predicted.density'))
    result={'case':case,'frozen_radials':frozen,'reference_norm':float(np.linalg.norm(ref.vector())),'residual_norm':float(np.linalg.norm(err.vector()))}
    (root/('full_map_'+label+'_result.json')).write_text(json.dumps(result,indent=2));print(result,flush=True)
if __name__=='__main__':run(sys.argv[1],len(sys.argv)>2 and sys.argv[2]=='frozen')
