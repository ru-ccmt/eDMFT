"""Isolated DFT regression: original inputs remain untouched."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess, shutil, os, time, json, hashlib

ROOT = Path(__file__).resolve().parent
WIEN = Path('/Users/haule/wien2k_24')
MATERIALS = Path('/Users/haule/projects/materials')
CASES = [
 ('FeSe', MATERIALS/'FeSe_again/FeSe_stress'/d, False)
 for d in ['DFT_RELA', 'DFT_RELA_m001', 'DFT_RELA_p001']
] + [
 ('SrVO3', MATERIALS/'SrVO3'/d, True)
 for d in ['DFT_center', 'DFT_m001', 'DFT_p001']
]
env = dict(os.environ, WIENROOT=str(WIEN), OMP_NUM_THREADS='2',
           VECLIB_MAXIMUM_THREADS='2', OPENBLAS_NUM_THREADS='1')
env['PATH'] = str(WIEN)+':/Users/haule/dbin:/opt/homebrew/bin:'+env['PATH']
for key in list(env):
    if key.startswith('WIEN_STRESS_'):
        del env[key]

def run(command, folder, label, extra=None):
    start=time.monotonic()
    print(f'START {folder.name} {label}', flush=True)
    with (folder/(label+'.log')).open('w') as log:
        result=subprocess.run(command,cwd=folder,env=env|dict(extra or {}),
                              stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'{folder}: {label} exit {result.returncode}')
    print(f'DONE {folder.name} {label} {time.monotonic()-start:.1f}s',flush=True)

def snapshot(folder,case,label):
    out=folder/'audit'/label
    out.mkdir(parents=True,exist_ok=True)
    for suffix in ['scf0','scf1','scf2','scfc','scfm','output2','clmval','clmcor','clmsum']:
        f=folder/(case+'.'+suffix)
        if f.exists():shutil.copy2(f,out/f.name)

def one(item):
    case,source,iscomplex=item
    folder=ROOT/'runs'/(case+'_'+source.name)
    shutil.copytree(source,folder,symlinks=False)
    snapshot(folder,case,'original')
    (folder/'audit/source.txt').write_text(str(source)+'\n')
    # Prevent old DMFT or optional semicore records contaminating a pure DFT sum.
    for f in ['Eorb_imp.dat',case+'.scf1s',case+'.scf2s',case+'.qdmft']:
        if (folder/f).exists():
            (folder/f).rename(folder/'audit'/('unused_'+f))
    # All components must be present for a total stress. Only isolated input
    # copies are changed; XC, integration, RMT and basis settings are preserved.
    for suffix,mode in [('in0','STR'),('in2c' if iscomplex else 'in2','STRF')]:
        p=folder/(case+'.'+suffix)
        rows=p.read_text().splitlines(True)
        oldmode=rows[0].split()[0]
        rows[0]=rows[0].replace(oldmode,mode,1)
        p.write_text(''.join(rows))
    run([str(WIEN/'x'),'lapw0','-f',case],folder,'lapw0_new')
    run([str(WIEN/'x'),'lapw1','-f',case]+(['-c'] if iscomplex else []),folder,'lapw1_new')
    exe='lapw2c' if iscomplex else 'lapw2'
    center=not (source.name.endswith('m001') or source.name.endswith('p001'))
    if center:
        run([str(WIEN/exe),'lapw2.def'],folder,'lapw2_before',{'WIEN_STRESS_DEBUG':'1'})
        snapshot(folder,case,'before')
    run([str(ROOT/'SRC_lapw2'/exe),'lapw2.def'],folder,'lapw2_corrected',{'WIEN_STRESS_DEBUG':'1'})
    snapshot(folder,case,'corrected_trace')
    run([str(WIEN/'x'),'lcore','-f',case],folder,'lcore_new')
    run(['/Users/haule/miniconda3-clean/envs/work311/bin/python',
         '/Users/haule/dbin/merge_scf.py','-f',case,'--replace'],folder,'merge_new')
    run([str(WIEN/'x'),'mixer','-f',case],folder,'mixer_new')
    snapshot(folder,case,'corrected_total')
    if center:
        run([str(ROOT/'SRC_lapw2'/exe),'lapw2.def'],folder,'lapw2_full',
            {'WIEN_STRESS_G2_MODE':'FULL','WIEN_STRESS_DEBUG':'1'})
        snapshot(folder,case,'corrected_full')
        # Leave the copied run at its tested production TRACE result.
        for suffix in ['scf2','output2','clmval']:
            shutil.copy2(folder/'audit/corrected_trace'/(case+'.'+suffix),folder/(case+'.'+suffix))
    return str(folder)

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(one,CASES))
    (ROOT/'runs.json').write_text(json.dumps(results,indent=2)+'\n')
