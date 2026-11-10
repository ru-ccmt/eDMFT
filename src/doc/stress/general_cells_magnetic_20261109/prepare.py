from pathlib import Path
import shutil,sys,json,os,subprocess,time
ROOT=Path('/private/tmp/edmft_general_cells_20261109')
SRC=Path('/Users/haule/projects/eDMFT_new/src/python')
for sub in ['base','stage']:
 p=ROOT/sub;p.mkdir(exist_ok=True)
 for f in SRC.glob('*.py'):shutil.copy2(f,p/f.name)
 shutil.copytree(SRC/'tests',p/'tests',dirs_exist_ok=True)
sys.path.insert(0,str(SRC))
from dmft2_step import local_definition
jobs=[]
for mat,ds,c in [('MnO',['MnO_center','MnO_m_0.001','MnO_p_0.001'],'MnO'),('MnTe',['dmft','dmft_m_0.001','dmft_p_0.001'],'0_800_MnTe')]:
 for d in ds:
  src=Path('/Users/haule/projects/materials')/mat/d; dst=ROOT/'native'/mat/d
  dst.mkdir(parents=True,exist_ok=False)
  suffixes=['struct','in0','in1','in2','inc','inm','klist','kgen','rotlm','indmfl','indmfldn','indmfi','indmf2','indmf2dn','clmsum']
  names=[c+'.'+x for x in suffixes]+['EF.dat','projectorw.dat','Eorb_imp.dat']+[f.name for f in src.glob('sig.inp*')]
  for name in names:
   f=src/name
   if f.is_file():shutil.copy2(f,dst/name)
  for label in ['lapw1','dmft2']+(['dmft2dn'] if mat=='MnTe' else []):
   text=local_definition((src/(label+'.def')).read_text());(dst/(label+'.def')).write_text(text)
  for s in ['in0','in2']:
   f=dst/(c+'.'+s);lines=f.read_text().splitlines();fields=lines[0].split();fields[0]='TOT';lines[0]=' '.join(fields);f.write_text('\n'.join(lines)+'\n')
  for s in ['indmf2']+(['indmf2dn'] if mat=='MnTe' else []):
   f=dst/(c+'.'+s);lines=f.read_text().splitlines();i=next(j for j,l in enumerate(lines) if l.strip() and not l.startswith('#'));words=lines[i].split();words[1]=(dst/'EF.dat').read_text().split()[0];words[2]='0';lines[i]=' '.join(words);f.write_text('\n'.join(lines)+'\n')
  jobs.append(dict(case=c,material=mat,name=d,source=str(src),directory=str(dst)))
(ROOT/'native_jobs.json').write_text(json.dumps(jobs,indent=2)+'\n')
print('Created stage and six independent native snapshots')
