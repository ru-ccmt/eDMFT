"""Independent native endpoints for each new potential sector and their sum."""
from pathlib import Path
import os,sys,shutil,subprocess,re
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from wstruct import Struct
from wien_response_density import FLOAT
case=sys.argv[1]
source=Path('/private/tmp/edmft_coupled_valence_response_20261108')/dict(FeSe='FeSe_b509_one',SrVO3='SrVO3_b703_one')[case]/'central'
s=Struct();s.ReadStruct(str(source/(case+'.struct')))
mode=case=='SrVO3'

def perturb_spherical(text,scale):
    patches=[]
    for atom in re.finditer(r'ATOMNUMBER\s*=\s*(\d+)',text):
        a=int(atom[1])-1;n=s.jrj[a];r=s.r0[a]*np.exp(np.log(s.rmt[a]/s.r0[a])*np.arange(n)/(n-1))
        end=text.find('ATOMNUMBER',atom.end());end=len(text) if end<0 else end
        channel=re.search(r'VLM\(R\) FOR L\s*0\s+M=\s*0',text[atom.end():end]);offset=atom.end()+channel.end()
        dv=.07*r*np.sin(1.3*r)
        for i,m in enumerate(list(FLOAT.finditer(text[offset:end]))[:n]):
            patches.append((offset+m.start(),offset+m.end(),float(m.group().replace('D','E'))+scale*dv[i]))
    for a,b,v in reversed(patches):text=text[:a]+f'{v:{b-a}.12E}'+text[b:]
    return text

def perturb_nonspherical(text,scale):
    patches=[]
    for atom in re.finditer(r'ATOMNUMBER\s*=\s*(\d+)',text):
        a=int(atom[1])-1
        end=text.find('ATOMNUMBER',atom.end());end=len(text) if end<0 else end
        for ch in re.finditer(r'VLM\(R\) FOR L\s*(-?\d+)\s+M=\s*(\d+)',text[atom.end():end]):
            offset=atom.end()+ch.end()
            for m in list(FLOAT.finditer(text[offset:end]))[:s.jrj[a]]:
                patches.append((offset+m.start(),offset+m.end(),float(m.group().replace('D','E'))*(1+.1*scale)))
    for a,b,v in reversed(patches):text=text[:a]+f'{v:{b-a}.12E}'+text[b:]
    return text

def perturb_ir(text,scale):
    lines=text.splitlines()
    j=next(i for i,l in enumerate(lines) if 'NUMBER OF PW' in l);n=int(lines[j].split()[0])
    for i in range(j+1,j+1+n):
        # Native potential records are 3X,3I5,2E19.12. Scale both components,
        # preserving every space-group and real-field constraint of the base.
        v=float(lines[i][18:37]);w=float(lines[i][37:56])
        lines[i]=lines[i][:18]+f'{v*(1+.1*scale):19.12E}{w*(1+.1*scale):19.12E}'
    return '\n'.join(lines)+'\n'

env=os.environ.copy()
for key in list(env):
    if key.startswith(('EDMFT_','DFTKS_','WIEN_')):del env[key]
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
exe=Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if mode else 'lapw1')
sectors=sys.argv[2].split(',') if len(sys.argv)>2 else ['nonspherical','interstitial','joint']
steps=list(map(float,sys.argv[3].split(','))) if len(sys.argv)>3 else [.001,.0005]
for sector in sectors:
    for step in steps:
        for sign in [-1,1]:
            folder=ROOT/case/sector/(('minus_' if sign<0 else 'plus_')+str(step))
            if (folder/'done').exists():continue
            folder.mkdir(parents=True,exist_ok=True)
            for p in source.iterdir():
                if p.is_file() and (p.name.startswith(case+'.') or p.name in ['EF.dat','Eorb_imp.dat','projectorw.dat','dmft2.def','lapw1.def'] or p.name.startswith('sig.inp')):
                    # Outputs are produced afresh by native programs below.
                    if p.suffix in ['.vector','.clmval','.output1','.output2'] or p.name.startswith(case+'.scf'):continue
                    shutil.copy2(p,folder/p.name)
            vns=(source/(case+'.vns')).read_text();vsp=(source/(case+'.vsp')).read_text()
            if sector in ['nonspherical','joint']:vns=perturb_nonspherical(vns,sign*step)
            if sector in ['interstitial','joint']:vns=perturb_ir(vns,sign*step)
            if sector=='joint':vsp=perturb_spherical(vsp,sign*step)
            (folder/(case+'.vns')).write_text(vns);(folder/(case+'.vsp')).write_text(vsp)
            # Fresh normalization at every endpoint, fixed Sigma/mu/energies.
            for name,cmd,extra in [('lapw1',[str(exe),'lapw1.def'],{'WIEN_LAPW1_MATRIX_DUMP':str(folder/'native')}),
                                  ('dmft2',['/Users/haule/dbin/dmft2_response','dmft2.def'],{'EDMFT_RESPONSE_EXPORT':str(folder/'response')})]:
                with (folder/(name+'.log')).open('w') as log:
                    proc=subprocess.run(cmd,cwd=folder,env=env|extra,stdout=log,stderr=subprocess.STDOUT)
                proc.check_returncode()
                if name=='dmft2' and 'DMFT2 END' not in (folder/(name+'.log')).read_text():raise RuntimeError(str(folder))
            (folder/'done').touch();print(case,sector,step,sign,'complete',flush=True)
