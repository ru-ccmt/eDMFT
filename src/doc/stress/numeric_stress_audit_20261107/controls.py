from pathlib import Path
import os,subprocess,shutil,re,json,concurrent.futures
root=Path('/private/tmp/numeric_stress_audit_20261107');env=os.environ|{'WIENROOT':'/Users/haule/wien2k_24','WIEN_DMFT_ROOT':'/Users/haule/dbin','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'}
for k in list(env):
 if k.startswith(('EDMFT_STRESS_','WIEN_STRESS_','EDMFT_PROJECTOR_')):del env[k]
def tag(p,t):return float(re.findall('^'+re.escape(t)+r'\s*(.*)$',p.read_text(),re.M)[-1].split('=')[-1].split()[0])
def run(name):
 d=root/name;out={};base=d/'baseline'
 def one(label):
  with (d/(label+'.log')).open('w') as f:p=subprocess.run(['/opt/homebrew/bin/mpirun','-n','4','/Users/haule/dbin/dmft2','dmft2.def'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
  assert p.returncode==0,(name,label)
  shutil.copy2(d/'SrVO3.scf2',d/(label+'.scf2'))
  a={t:tag(d/'SrVO3.scf2',t) for t in [':SUM  :',':VSUM :',':DRHO']};a['mu_eV']=float((d/'EF.dat').read_text().split()[0]);print(name,label,a,flush=True);return a
 for f in ['sig.inp1','SOlapm.dat','SrVO3.indmf2','EF.dat']:shutil.copy2(base/f,d/f)
 (d/'SOlapm.dat').unlink();out['recomputed_normalization']=one('new_normalization')
 for f in ['sig.inp1','SOlapm.dat','SrVO3.indmf2','EF.dat']:shutil.copy2(base/f,d/f)
 if name=='minus':
  mu=float((base/'EF.dat').read_text().split()[0]);orig=(base/'SrVO3.indmf2').read_text().splitlines();parts=orig[-1].split();parts[2]='0'
  for shift in [-0.001,0.001]:
   parts[1]=str(mu+shift);(d/'EF.dat').write_text(str(mu+shift)+'\n');(d/'SrVO3.indmf2').write_text(orig[0]+'\n'+' '.join(parts)+'\n');out['mu_shift_'+str(shift)]=one('mu_'+str(shift))
 for f in ['sig.inp1','SOlapm.dat','SrVO3.indmf2','EF.dat']:shutil.copy2(base/f,d/f)
 return name,out
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:x=dict(ex.map(run,['minus','plus']))
(root/'controls.json').write_text(json.dumps(x,indent=2)+'\n')
