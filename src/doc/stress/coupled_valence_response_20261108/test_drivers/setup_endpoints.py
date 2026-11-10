"""Isolated physical-potential endpoints with explicitly fixed radial energies."""
from pathlib import Path
import os,sys,shutil,re,subprocess
import numpy as np
from scipy.io import FortranFile
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from wien_response_density import FLOAT
from dmft2_step import local_definition
from wstruct import Struct
case=sys.argv[1];nbands=int(sys.argv[2]) if len(sys.argv)>2 else None
one_k=len(sys.argv)>3 and sys.argv[3]=='one'
source=Path('/private/tmp/edmft_spatial_response_20261108')/case/'central'
root=ROOT/((case if nbands is None else case+'_b'+str(nbands))+('_one' if one_k else ''));root.mkdir(exist_ok=True)
s=Struct();s.ReadStruct(str(source/(case+'.struct')))
cplx=case=='SrVO3'
in1=source/(case+('.in1c' if cplx else '.in1'));lines=in1.read_text().splitlines();position=2
with FortranFile(source/(case+'.vector')) as f:
 for a in range(s.nat):
  e=f.read_record('f8');elo=f.read_record('f8').reshape(4,-1,order='F');apw=e>150;e[apw]-=200
  head=lines[position].split();nr=int(head[1]);lines[position]=f'{e[-1]:.17g} {nr} {head[2]}';position+=1;counts={}
  for j in range(nr):
   fields=lines[position].split();l=int(fields[0]);count=counts.get(l,0);counts[l]=count+1
   value=e[l] if count==0 else elo[l,count-1+int(apw[l])]
   line=f' {l:1d}{value:10.6f}{0.:10.5f}{fields[3]:4s}{int(fields[4]):2d}'
   assert len(line)==28,(line,len(line))
   lines[position]=line;position+=1
if nbands:
 fields=lines[-1].split();assert fields[6]=='FIXED_NBANDS';fields[5]=str(nbands);lines[-1]=' '.join(fields)
frozen='\n'.join(lines)+'\n'
env=os.environ.copy()
for key in list(env):
 if key.startswith(('EDMFT_','DFTKS_','WIEN_LAPW1_','WIEN_STRESS_')):del env[key]
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
exe=Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if cplx else 'lapw1')
base=Path('/private/tmp/edmft_native_operator_20261108/native_probe')/case/'center/lapw1.def'
def perturb(text,amplitude):
 patches=[]
 for block in re.finditer(r'ATOMNUMBER\s*=\s*(\d+)',text):
  a=int(block[1])-1;n=s.jrj[a];r=s.r0[a]*np.exp(np.log(s.rmt[a]/s.r0[a])/(n-1)*np.arange(n));dv=.07*r*np.sin(1.3*r)
  end=text.find('ATOMNUMBER',block.end());end=len(text) if end<0 else end
  channel=re.search(r'VLM\(R\) FOR L\s*0\s+M=\s*0',text[block.end():end]);offset=block.end()+channel.end()
  for i,match in enumerate(list(FLOAT.finditer(text[offset:end]))[:n]):
   value=float(match.group().replace('D','E'))+amplitude*dv[i]
   patches.append((offset+match.start(),offset+match.end(),value))
 for a,b,value in reversed(patches):text=text[:a]+f'{value:{b-a}.12E}'+text[b:]
 return text
for step,sign in [(0,0)]+[(h,t) for h in ([.001,.0005,.00025] if one_k else [.002,.001,.0005]) for t in [-1,1]]:
 name='central' if sign==0 else ('m_' if sign<0 else 'p_')+f'{step:.0e}'
 d=root/name
 if (d/'done').exists():continue
 d.mkdir(exist_ok=True)
 for p in source.iterdir():
  if p.is_file() and (p.name.startswith(case+'.') or p.name in ['EF.dat','Eorb_imp.dat','projectorw.dat','dmft2.def'] or p.name.startswith('sig.inp')):
   shutil.copy2(p,d/p.name)
 (d/in1.name).write_text(frozen);(d/'lapw1.def').write_text(local_definition(base.read_text()))
 if one_k:
  p=d/(case+'.klist');p.write_text(p.read_text().splitlines()[0]+'\nEND\n')
 (d/(case+'.vsp')).write_text(perturb((source/(case+'.vsp')).read_text(),sign*step))
 if (d/'SOlapm.dat').exists():(d/'SOlapm.dat').unlink()
 for label,command,extra in [('lapw1',[str(exe),'lapw1.def'],{'WIEN_LAPW1_MATRIX_DUMP':str(d/'native')}),
            ('dmft2',['/Users/haule/dbin/dmft2_response','dmft2.def'],{'EDMFT_RESPONSE_EXPORT':str(d/'response')})]:
  with (d/(label+'.log')).open('w') as log:
   proc=subprocess.run(command,cwd=d,env=env|extra,stdout=log,stderr=subprocess.STDOUT)
  proc.check_returncode()
  if label=='dmft2' and 'DMFT2 END' not in (d/(label+'.log')).read_text():raise RuntimeError(str(d))
 (d/'done').touch();print(case,name,'complete',flush=True)
