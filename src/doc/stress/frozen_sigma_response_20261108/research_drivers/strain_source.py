"""Obtain the strain source from two SINGLE, unrelaxed charge evaluations.

These are not strained SCF calculations: each sign runs LAPW0, LAPW1,
dmft2, and lcore exactly once on the affinely transported central density.
The antisymmetric output is the source s=b+chi*v_aff.  Later screening uses
the analytic central resolvent response, not additional strained SCF runs.
"""
import os,sys,re,json,shutil,time,subprocess,concurrent.futures
from pathlib import Path
import numpy as np
from setup import ROOT,SOURCES,ENV,EXE
from wien_response_density import Density
from check_charge import integrate

WIEN=Path('/Users/haule/wien2k_24')
def command(d,label,cmd,extra={}):
    t=time.monotonic()
    with (d/(label+'.log')).open('w') as f:p=subprocess.run(cmd,cwd=d,env=ENV|extra,stdout=f,stderr=subprocess.STDOUT)
    log=(d/(label+'.log')).read_text(errors='replace')
    if p.returncode or any(word in log for word in ['Fortran runtime error','Error termination','error: command','ERROR STOP']):
        raise RuntimeError((label,p.returncode,d))
    return time.monotonic()-t

def one(case,h,name):
    src=SOURCES[case];d=ROOT/case/name;d.mkdir(parents=True,exist_ok=True)
    for p in src.iterdir():
        if not p.is_file():continue
        if (p.name.startswith(case+'.') and p.suffix in ['.struct','.in0','.in1','.in1c','.in2','.in2c','.inc','.inm','.klist','.kgen','.rotlm','.clmsum','.clmcor','.indmfl','.indmfi','.indmf2']) or p.name in ['EF.dat','projectorw.dat','Eorb_imp.dat','dmft2.def'] or re.fullmatch(r'sig\.inp\d+',p.name):shutil.copy2(p,d/p.name)
    a=(d/(case+'.struct')).read_text().splitlines();cell=np.array([float(a[3][i:i+10]) for i in [0,10,20]])
    a[3]=''.join(f'{x:10.6f}' for x in cell*(1+h))+a[3][30:];newcell=np.array([float(a[3][i:i+10]) for i in [0,10,20]]);lam=float(np.mean(newcell/cell))
    (d/(case+'.struct')).write_text('\n'.join(a)+'\n')
    rho=Density(src/(case+'.clmsum'),src/(case+'.struct'));rho.transport(lam).write(d/(case+'.clmsum'))
    p=d/(case+'.in0');a=p.read_text().splitlines();a[0]='TOT '+a[0][4:];a[1]=a[1].replace('NR2V','R2V ');p.write_text('\n'.join(a)+'\n')
    for suffix in ['in2','in2c']:
        p=d/(case+'.'+suffix)
        if p.exists():a=p.read_text().splitlines();a[0]='TOT '+a[0][4:];p.write_text('\n'.join(a)+'\n')
    p=d/(case+'.indmf2');a=p.read_text().splitlines();i=next(i for i,l in enumerate(a) if l.strip() and not l.startswith('#'));a[i]=f'c {(src/"EF.dat").read_text().split()[0]} 1 1.0 .T. 0.02';p.write_text('\n'.join(a)+'\n')
    lapwdef=Path('/private/tmp/projector_actual_stress_20261107/complete_center/lapw1.def') if case=='SrVO3' else src/'lapw1.def'
    (d/'lapw1.def').write_text(lapwdef.read_text().replace('_x',''))
    p=d/'dmft2.def';p.write_text(p.read_text().replace('_x',''))
    if (d/'SOlapm.dat').exists():(d/'SOlapm.dat').unlink()
    timing={}
    timing['lapw0']=command(d,'lapw0',[str(WIEN/'x'),'lapw0','-f',case])
    timing['lapw1']=command(d,'lapw1',[str(WIEN/('lapw1c' if case=='SrVO3' else 'lapw1')),'lapw1.def'])
    timing['dmft2']=command(d,'dmft2',['/opt/homebrew/bin/mpirun','-n','4',str(EXE),'dmft2.def'])
    if not (d/(case+'.clmval')).exists():raise RuntimeError('dmft2 did not produce density')
    timing['lcore']=command(d,'lcore',[str(WIEN/'x'),'lcore','-f',case])
    val=Density(d/(case+'.clmval'),d/(case+'.struct'));core=Density(d/(case+'.clmcor'),d/(case+'.struct'))
    total=rho.zeros().add(val).add(core)
    # Match MIXER's default normalization: charge leakage is put into the
    # interstitial Fourier coefficients, not the core or MT valence density.
    target=round(integrate(rho,src/(case+'.struct'))[2])
    mt,ir,_=integrate(total,d/(case+'.struct'))
    for g in total.pw:total.pw[g]*=(target-mt)/ir
    total.write(d/(case+'.clmout'))
    pullback=total.transport(lam,pullback=True)
    mt,ir,_=integrate(pullback,src/(case+'.struct'))
    for g in pullback.pw:pullback.pw[g]*=(target-mt)/ir
    pullback.write(d/(case+'.qout'))
    result={'actual_step':lam-1,'timings':timing,'EF_eV':float((d/'EF.dat').read_text().split()[0])}
    (d/'source.json').write_text(json.dumps(result,indent=2));print(case,name,result,flush=True)
    return pullback,lam-1

def run(case,h):
    label=f'{h:.7f}'
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as p:
        jobs=[p.submit(one,case,h,'source_p_'+label),p.submit(one,case,-h,'source_m_'+label)]
        (plus,hp),(minus,hm)=[j.result() for j in jobs]
    source=plus.add(minus,-1).scale(1/(hp-hm));source.write(ROOT/case/f'source_{label}.density')
    return source
if __name__=='__main__':run(sys.argv[1],float(sys.argv[2]) if len(sys.argv)>2 else 0.00025)
