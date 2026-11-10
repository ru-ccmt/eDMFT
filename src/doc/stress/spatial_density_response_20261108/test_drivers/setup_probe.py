from pathlib import Path
import shutil,os,subprocess,sys
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from dmft2_step import local_definition

root=Path(__file__).resolve().parent
for case in ['SrVO3','FeSe']:
    source=Path('/private/tmp/edmft_linear_response_20261108')/case/'central'
    vectors=Path('/private/tmp/edmft_native_operator_20261108/native_probe')/case/'center'
    dest=root/case/'central';dest.mkdir(parents=True,exist_ok=False)
    for p in source.iterdir():
        if p.is_file() and (p.suffix in ['.struct','.in1','.in1c','.in2','.in2c','.indmf2','.indmfl','.indmfi','.inso','.kgen','.rotlm','.vsp','.vns'] or p.name.startswith('sig.inp') or p.name in ['EF.dat','projectorw.dat','Eorb_imp.dat','SOlapm.dat']):
            shutil.copy2(p,dest/p.name)
    for p in vectors.iterdir():
        if p.suffix in ['.struct','.vsp','.vns','.vector','.energy','.in1','.in1c','.klist','.nsh'] or p.name.startswith('native.'):
            shutil.copy2(p,dest/p.name)
    (dest/'dmft2.def').write_text(local_definition((source/'dmft2.def').read_text()))
    p=dest/(case+'.indmf2');s=p.read_text().splitlines();t=s[1].split();t[2]='0';s[1]=' '.join(t);p.write_text('\n'.join(s)+'\n')
    for suffix in ['.in2','.in2c']:
        p=dest/(case+suffix)
        if p.exists() and p.stat().st_size:
            s=p.read_text().splitlines();s[0]='TOT';p.write_text('\n'.join(s)+'\n')
    env=os.environ.copy()
    for k in list(env):
        if k.startswith(('EDMFT_','WIEN_LAPW1','DFTKS_')):del env[k]
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',EDMFT_RESPONSE_EXPORT=str(dest/'response'))
    with (dest/'run.log').open('w') as log:
        p=subprocess.run(['/Users/haule/dbin/dmft2_response','dmft2.def'],cwd=dest,env=env,stdout=log,stderr=subprocess.STDOUT)
    p.check_returncode()
    if 'DMFT2 END' not in (dest/'run.log').read_text(): raise RuntimeError('native failure '+case)
    print(case,'prepared',flush=True)
