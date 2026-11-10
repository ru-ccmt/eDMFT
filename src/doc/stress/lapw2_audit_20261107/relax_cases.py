"""Remove the small charge residual after switching from dftKS to WIEN lapw1."""
from run_cases import ROOT,WIEN,CASES,run,snapshot
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil,re,json

def setmode(folder,case,iscomplex,stress):
    for ext,mode in [('in0','STR' if stress else 'TOT'),
                     ('in2c' if iscomplex else 'in2','STRF' if stress else 'TOT')]:
        p=folder/(case+'.'+ext);rows=p.read_text().splitlines(True)
        rows[0]=rows[0].replace(rows[0].split()[0],mode,1);p.write_text(''.join(rows))

def cycle(folder,case,iscomplex,label,stress=False):
    run([str(WIEN/'x'),'lapw0','-f',case],folder,label+'_lapw0')
    run([str(WIEN/'x'),'lapw1','-f',case]+(['-c'] if iscomplex else []),folder,label+'_lapw1')
    run([str(ROOT/'SRC_lapw2'/('lapw2c' if iscomplex else 'lapw2')),'lapw2.def'],
        folder,label+'_lapw2',{'WIEN_STRESS_DEBUG':'1'} if stress else None)
    run([str(WIEN/'x'),'lcore','-f',case],folder,label+'_lcore')
    run(['/Users/haule/miniconda3-clean/envs/work311/bin/python','/Users/haule/dbin/merge_scf.py','-f',case,'--replace'],folder,label+'_merge')
    run([str(WIEN/'x'),'mixer','-f',case],folder,label+'_mixer')

def relax(item):
    case,source,iscomplex=item
    folder=ROOT/'runs'/(case+'_'+source.name)
    center=not (source.name.endswith('m001') or source.name.endswith('p001'))
    ready=['corrected_total','corrected_trace']+(['before','corrected_full'] if center else [])
    if (folder/'audit/relaxation.json').exists() and all((folder/'audit'/d/(case+'.scfm')).exists() for d in ready):
        return str(folder)
    # An old converged MSR1 history can suppress the new residual caused by
    # changing the lapw1 implementation. Restart history, retaining a backup.
    hist=folder/'audit/mixing_history_before_restart'
    hist.mkdir(exist_ok=True)
    for p in folder.glob('*.broyd*'):
        if not (hist/p.name).exists():shutil.move(str(p),str(hist/p.name))
    initial=folder/'audit/initial_cycle';initial.mkdir(exist_ok=True)
    for name in ['before','corrected_trace','corrected_total','corrected_full']:
        p=folder/'audit'/name
        if p.exists() and not (initial/name).exists():shutil.move(str(p),str(initial/name))
    setmode(folder,case,iscomplex,False)
    prev=None;history=[]
    for it in range(1,41):
        cycle(folder,case,iscomplex,f'relax_{it:02d}')
        s=(folder/(case+'.scfm')).read_text()
        dis=float(re.findall(r'^:DIS.*$',s,re.M)[-1].split()[-1])
        ene=float(re.findall(r'^:ENE.*$',s,re.M)[-1].split()[-1])
        de=None if prev is None else ene-prev
        history.append(dict(iteration=it,charge_distance=dis,energy=ene,energy_change=de))
        print(f'CONVERGENCE {folder.name} {history[-1]}',flush=True)
        if prev is not None and dis<=2e-7 and abs(de)<=5e-8:break
        prev=ene
    else:raise RuntimeError(f'{folder}: charge did not converge')
    (folder/'audit/relaxation.json').write_text(json.dumps(history,indent=2)+'\n')
    setmode(folder,case,iscomplex,True)
    cycle(folder,case,iscomplex,'final_stress',True)
    snapshot(folder,case,'corrected_trace');snapshot(folder,case,'corrected_total')
    if center:
        exe='lapw2c' if iscomplex else 'lapw2'
        run([str(WIEN/exe),'lapw2.def'],folder,'final_before',{'WIEN_STRESS_DEBUG':'1'})
        snapshot(folder,case,'before')
        run([str(ROOT/'SRC_lapw2'/exe),'lapw2.def'],folder,'final_full',
            {'WIEN_STRESS_G2_MODE':'FULL','WIEN_STRESS_DEBUG':'1'})
        snapshot(folder,case,'corrected_full')
        for ext in ['scf2','output2','clmval']:
            shutil.copy2(folder/'audit/corrected_trace'/(case+'.'+ext),folder/(case+'.'+ext))
    return str(folder)

if __name__=='__main__':
    assert (ROOT/'runs.json').exists(),'Initial regression still running'
    with ThreadPoolExecutor(max_workers=2) as pool:
        result=list(pool.map(relax,CASES))
    (ROOT/'relaxed_runs.json').write_text(json.dumps(result,indent=2)+'\n')
