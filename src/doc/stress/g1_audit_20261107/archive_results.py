from pathlib import Path
import shutil,json,hashlib,os
root=Path(__file__).resolve().parent
base=Path('/Users/haule/projects/eDMFT_new/src/doc/stress')
dest=base/'g1_audit_20261107';dest.mkdir(exist_ok=True)
# Keep the note self-contained: its final .tex includes the material addendum.
for suffix in ['tex','pdf']:
 p=root/('edmft_stress_eq19_to_g1_derivation_20261107.'+suffix)
 shutil.copy2(p,base/p.name)
for p in root.iterdir():
 if p.is_file() and p.suffix in ['.py','.json','.log','.patch','.txt']:
  shutil.copy2(p,dest/p.name)
for directory in ['before','prepared']:
 shutil.copytree(root/directory,dest/directory,dirs_exist_ok=True)
for f in (root/'runs').iterdir():
 if not f.is_dir():continue
 out=dest/'results'/f.name;out.mkdir(parents=True,exist_ok=True)
 for p in f.iterdir():
  if p.is_file() and (p.suffix in ['.def','.struct','.in0','.in1','.in1c','.in2','.in2c','.indmfl','.indmf2','.inc','.inm','.scf0','.scf1','.scf2','.scfc','.scfm','.log'] or p.name in ['EF.dat','SOlapm.dat']):shutil.copy2(p,out/p.name)
 for p in (f/'audit').rglob('*'):
  if p.is_file() and p.suffix!='.clmval':
   t=out/'audit'/p.relative_to(f/'audit');t.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,t)
 # Density equality is recorded by hash; large vectors and duplicate density
 # data stay in scratch. All source/test inputs needed for unit checks remain.
files=[]
for p in sorted(dest.rglob('*')):
 if p.is_file() and p.name!='MANIFEST.sha256':files.append(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(dest)))
(dest/'MANIFEST.sha256').write_text('\n'.join(files)+'\n')
print('Archived',len(files),'files in',dest)
