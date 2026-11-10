from pathlib import Path
import os,sys,shutil,subprocess,re
import numpy as np
root=Path(__file__).resolve().parent
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from wien_response_density import FLOAT
from wstruct import Struct
case,part=sys.argv[1:3]
original=Path('/Users/haule/projects/materials')/dict(FeSe='FeSe_again/FeSe_stress/FeSe_RELA',SrVO3='SrVO3/V_plus_force')[case]
fixture=Path('/private/tmp/edmft_coupled_valence_response_20261108')/dict(FeSe='FeSe_b509_one',SrVO3='SrVO3_b703_one')[case]/'central'
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
for a in [0,.01,-.01,.005,-.005]:
 dst=root/case/part/str(a);dst.mkdir(parents=True,exist_ok=True)
 if part=='xc':
  source=Path('/private/tmp/edmft_hartree_periodic_20261108')/case/str(a)
  for name in ['lapw0.def','xc_funcs.h',case+'.struct',case+'.in0',case+'.inm',case+'.clmsum']:
   shutil.copy2(source/name,dst/name)
  exe=root/'native/SRC_lapw0/lapw0';definition='lapw0.def'
 else:
  for name in [case+'.struct',case+'.vsp',case+'.vns']:shutil.copy2(fixture/name,dst/name)
  for name in [case+'.inc','lcore.def']:shutil.copy2(original/name,dst/name)
  s=Struct();s.ReadStruct(str(dst/(case+'.struct')))
  text=(dst/(case+'.vsp')).read_text();changes=[]
  for atom in re.finditer(r'ATOMNUMBER\s*=\s*(\d+)',text):
   ia=int(atom[1])-1;n=s.jrj[ia];r=np.geomspace(s.r0[ia],s.rmt[ia],n)
   end=text.find('ATOMNUMBER',atom.end());end=len(text) if end<0 else end
   ch=re.search(r'VLM\(R\) FOR L\s*0\s+M=\s*0',text[atom.end():end]);offset=atom.end()+ch.end()
   for j,m in enumerate(list(FLOAT.finditer(text[offset:end]))[:n]):
    changes.append((offset+m.start(),offset+m.end(),float(m.group().replace('D','E'))+a*.07*r[j]*np.sin(1.3*r[j])))
  for lo,hi,val in changes[::-1]:text=text[:lo]+f'{val:{hi-lo}.12E}'+text[hi:]
  (dst/(case+'.vsp')).write_text(text)
  (dst/'core_native.bin').unlink(missing_ok=True)
  exe=root/'native/SRC_lcore/lcore';definition='lcore.def'
 with (dst/'run.log').open('w') as log:
  p=subprocess.run([str(exe),definition],cwd=dst,env=env,stdout=log,stderr=subprocess.STDOUT)
 p.check_returncode()
 if not (dst/('xc_native.bin' if part=='xc' else 'core_native.bin')).exists():raise RuntimeError(str(dst))
 print(case,part,a,'complete',flush=True)
