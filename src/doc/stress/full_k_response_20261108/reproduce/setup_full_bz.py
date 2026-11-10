"""Isolated native three-k endpoints; no material directory is modified."""
from pathlib import Path
import sys,os,shutil,subprocess,re,json
root=Path(__file__).resolve().parent/'full_bz';case=sys.argv[1]
old=Path('/private/tmp/edmft_xc_core_response_20261108')/case
source=old/'charge/0';nr={'FeSe':66,'SrVO3':43}[case]
material={'FeSe':Path('/Users/haule/projects/materials/FeSe_again/FeSe_stress/FeSe_RELA'),
          'SrVO3':Path('/Users/haule/projects/materials/SrVO3/V_plus_force')}[case]
lines=(material/(case+'.klist')).read_text().splitlines()
chosen=[line for line in lines if line.strip() and line.strip()!='END']
klist='\n'.join(f'{i+1:10d}'+l[10:] for i,l in enumerate(chosen))+'\nEND\n'
env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))}
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
exe=Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if case=='SrVO3' else 'lapw1')
for step in [0,.001,-.001,.0005,-.0005]:
    folder=root/case/'charge'/str(step);folder.mkdir(parents=True,exist_ok=True)
    if (folder/'done').exists():continue
    for p in source.iterdir():
        if not p.is_file():continue
        if not (p.name.startswith(case+'.') or p.name in ['EF.dat','Eorb_imp.dat','projectorw.dat','dmft2.def','lapw1.def'] or p.name.startswith('sig.inp')):continue
        if p.suffix in ['.vector','.clmval','.output1','.output2'] or p.name.startswith(case+'.scf'):continue
        shutil.copy2(p,folder/p.name)
    for name in [case+'.vsp',case+'.vns',case+'.clmsum']:
        shutil.copy2(old/'charge'/str(step)/name,folder/name)
    (folder/(case+'.klist')).write_text(klist)
    in1=folder/(case+('.in1c' if case=='SrVO3' else '.in1'))
    text=in1.read_text();text=re.sub(r'\d+\s+FIXED_NBANDS',f'{nr} FIXED_NBANDS',text)
    assert f'{nr} FIXED_NBANDS' in text;in1.write_text(text)
    # Fresh global normalization for precisely these three k points, at each
    # endpoint; no central/old full-BZ SOlapm is reused.
    for name in ['SOlapm.dat','SOlapm.datup','SOlapm.datdn']:
        (folder/name).unlink(missing_ok=True)
    for label,cmd,extra in [('lapw1',[str(exe),'lapw1.def'],{'WIEN_LAPW1_MATRIX_DUMP':str(folder/'native')}),
                            ('dmft2',['/Users/haule/dbin/dmft2_response','dmft2.def'],{'EDMFT_RESPONSE_EXPORT':str(folder/'response')})]:
        with (folder/(label+'.log')).open('w') as log:
            subprocess.run(cmd,cwd=folder,env=env|extra,stdout=log,stderr=subprocess.STDOUT,check=True)
        if label=='dmft2' and 'DMFT2 END' not in (folder/(label+'.log')).read_text():raise RuntimeError(str(folder))
    (folder/'done').touch();print(case,step,'done',flush=True)
