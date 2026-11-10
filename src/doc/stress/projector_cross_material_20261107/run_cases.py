from pathlib import Path
import os,sys,re,shutil,subprocess,time,json,concurrent.futures
import numpy as np
from scipy.interpolate import CubicSpline
root=Path(__file__).parent
specs={'FeSe':('/Users/haule/projects/materials/FeSe_again/FeSe_stress','FeSe_RELA','FeSe_RELA_strain_m001','FeSe_RELA_strain_p001'), 'MnO':('/Users/haule/projects/materials/MnO','MnO_center','MnO_m_0.001','MnO_p_0.001')}
env=os.environ|{'WIENROOT':'/Users/haule/wien2k_24','WIEN_DMFT_ROOT':'/Users/haule/dbin','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'}
for k in list(env):
 if k.startswith(('EDMFT_','WIEN_STRESS_','DFTKS_','WIEN_LAPW1_')):del env[k]

def command(d,label,cmd,extra={}):
 t=time.monotonic();print('START',d.parent.name,d.name,label,flush=True)
 with (d/(label+'.log')).open('w') as f:p=subprocess.run(cmd,cwd=d,env=env|extra,stdout=f,stderr=subprocess.STDOUT)
 result={'rc':p.returncode,'seconds':time.monotonic()-t};print('END',d.parent.name,d.name,label,result,flush=True)
 if p.returncode:raise RuntimeError((d,label,result))
 return result

def copy_inputs(src,d,case):
 d.mkdir(exist_ok=True,parents=True)
 for p in src.iterdir():
  if not p.is_file():continue
  if (p.name.startswith(case+'.') and p.suffix in ['.struct','.in0','.in1','.in1c','.in2','.in2c','.inc','.inm','.klist','.kgen','.rotlm','.vsp','.vns','.clmsum','.clmcor','.scf0','.scfc','.indmfl','.indmfi','.indmf2']) or p.name in ['EF.dat','SOlapm.dat','projectorw.dat','Eorb_imp.dat','sig.inp'] or re.fullmatch(r'sig\.inp\d+',p.name):shutil.copy2(p,d/p.name)

def defs(d,case):
 t=Path('/private/tmp/projector_symmetry_fix_20261107/FeSe/lapw1.def').read_text().replace('FeSe',case).replace('_x','')
 (d/'lapw1.def').write_text(t)

def advect(text,struct,h,spherical):
 rad=[]
 for l in struct.splitlines():
  m=re.search(r'NPT=\s*(\d+)\s+R0=\s*([\d.]+)\s+RMT=\s*([\d.]+)',l)
  if m:rad.append(np.linspace(np.log(float(m[2])),np.log(float(m[3])),int(m[1])))
 a=text.splitlines();atom=-1;i=0
 while i<len(a):
  if 'ATOMNUMBER' in a[i]:atom=int(a[i].split('=')[1])-1
  if 'VLM(R) FOR L' not in a[i]:i+=1;continue
  i+=1
  while not a[i].strip():i+=1
  start=i;v=[]
  while len(v)<len(rad[atom]):v.extend(map(float,re.findall(r'[-+]?\d+\.\d+E[-+]\d+',a[i])));i+=1
  assert len(v)==len(rad[atom]);v=np.array(v);dv=CubicSpline(rad[atom],v)(rad[atom],1)-(v if spherical else 0)
  out=v+h*dv;a[start:i]=['   '+''.join(f'{x:19.12E}' for x in out[j:j+4]) for j in range(0,len(out),4)]
 return '\n'.join(a)+'\n'

def run(case):
 base,cs,ms,ps=specs[case];base=Path(base);r=root/case;r.mkdir(exist_ok=True);src=base/cs;central=r/'center';copy_inputs(src,central,case);defs(central,case)
 timings={}
 # Preserve the central physical potential and density as supplied.
 timings['central_lapw1']=command(central,'lapw1',['/Users/haule/wien2k_24/lapw1','lapw1.def'],{'WIEN_LAPW1_BASIS_VSP':str(central/(case+'.vsp'))})
 saved=(central/(case+'.indmf2')).read_bytes()
 timings['prepare']=command(central,'prepare',['/Users/haule/dbin/x_dmft.py','-d','dmft2'])
 (central/(case+'.indmf2')).write_bytes(saved)
 p=central/'dmft2.def';t=p.read_text().replace('_x','');p.write_text(t)
 timings['dmft2']=command(central,'dmft2',['/opt/homebrew/bin/mpirun','-n','4','/Users/haule/dbin/dmft2','dmft2.def'],{'EDMFT_STRESS_STERNHEIMER_DUMP':'1','EDMFT_STRESS_RESPONSE_DUMP':'1','EDMFT_ONEPARTICLE_DUMP':'1'})
 jobs=[]
 for name in ['complete_center','minus','plus','advection_minus','advection_plus']:
  d=r/name;copy_inputs(src,d,case);defs(d,case)
  if name=='complete_center':
   p=d/(case+'.in1');a=p.read_text().splitlines();a[-1]='K-VECTORS FROM UNIT:4  -10.0  1.0e10  0';p.write_text('\n'.join(a)+'\n')
  elif name in ['minus','plus']:
   es=base/(ms if name=='minus' else ps)
   for suf in ['vsp','vns']:shutil.copy2(es/(case+'.'+suf),d/(case+'.'+suf))
  else:
   h=(-1 if name.endswith('minus') else 1)*1.5e-5
   for suf in ['vsp','vns']:
    p=d/(case+'.'+suf);p.write_text(advect(p.read_text(),(d/(case+'.struct')).read_text(),h,suf=='vsp'))
  jobs.append(d)
 # Sequential within material; two materials execute concurrently.
 for d in jobs:timings[d.name]=command(d,'lapw1',['/Users/haule/wien2k_24/lapw1','lapw1.def'],{'WIEN_LAPW1_BASIS_VSP':str(central/(case+'.vsp'))})
 for name,srcname in [('center',cs),('minus',ms),('plus',ps)]:
  d=r/('r2v_'+name);copy_inputs(base/srcname,d,case)
  p=d/(case+'.in0');a=p.read_text().splitlines();a[0]='TOT '+a[0][4:];a[1]=a[1].replace('NR2V','R2V ');p.write_text('\n'.join(a)+'\n')
  timings['r2v_'+name]=command(d,'lapw0',['/Users/haule/wien2k_24/x','lapw0','-f',case])
 (r/'timings.json').write_text(json.dumps(timings,indent=2))
 return case
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as p:print(list(p.map(run,specs)),flush=True)
