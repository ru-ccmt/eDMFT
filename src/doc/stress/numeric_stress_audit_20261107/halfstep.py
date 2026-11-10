from pathlib import Path
import os,shutil,subprocess,time,json,re,concurrent.futures,hashlib
root=Path('/private/tmp/numeric_stress_audit_20261107');case='SrVO3';parent=Path('/Users/haule/projects/materials/SrVO3');env=os.environ|{'WIENROOT':'/Users/haule/wien2k_24','WIEN_DMFT_ROOT':'/Users/haule/dbin','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'}
for k in list(env):
 if k.startswith(('EDMFT_STRESS_','WIEN_STRESS_','EDMFT_PROJECTOR_')):del env[k]
def tag(p,t):return float(re.findall('^'+re.escape(t)+r'\s*(.*)$',p.read_text(),re.M)[-1].split('=')[-1].split()[0])
def run(spec):
 name,srcname,h=spec;src=parent/srcname;d=root/name;d.mkdir(exist_ok=True)
 for suf in ['struct','clmsum','in0','in1c','in2c','inc','inm','klist','kgen','rotlm','vsp','vns','clmcor','indmfl','indmfi','indmf2']:
  shutil.copy2(src/(case+'.'+suf),d/(case+'.'+suf))
 for f in ['EF.dat','Eorb_imp.dat','projectorw.dat','sig.inp','sig.inp1']:shutil.copy2(src/f,d/f)
 (d/'mpi_prefix.dat').write_text('/opt/homebrew/bin/mpirun -n 4\n')
 ls=(parent/'V_plus_force'/(case+'.struct')).read_text().splitlines();length=7.260328*(1+h);ls[3]=''.join(f'{z:10.6f}' for z in [length,length,length,90,90,90]);(d/(case+'.struct')).write_text('\n'.join(ls)+'\n')
 p=d/(case+'.in1c');ls=p.read_text().splitlines();ls[1]=re.sub(r'^\s*\S+',f'{8/(1+h):10.6f}',ls[1]);p.write_text('\n'.join(ls)+'\n')
 p=d/(case+'.in2c');ls=p.read_text().splitlines();ls[0]='TOT'+ls[0][3:]
 for i,l in enumerate(ls):
  if 'GMAX' in l:ls[i]=f'{12/(1+h):10.6f}       GMAX'
 p.write_text('\n'.join(ls)+'\n')
 def command(label,cmd):
  with (d/(label+'.log')).open('w') as f:proc=subprocess.run(cmd,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
  if proc.returncode:raise RuntimeError((name,label,proc.returncode))
 command('prepare_dmft2',['/Users/haule/dbin/x_dmft.py','-d','-p','dmft2'])
 shutil.copy2(src/(case+'.indmf2'),d/(case+'.indmf2'))
 history=[];stable=0
 for it in range(1,41):
  t=time.monotonic()
  command('lapw0',['/Users/haule/wien2k_24/x','lapw0','-f',case])
  command('lapw1',['/Users/haule/dbin/x_dmft.py','-p','lapw1'])
  if (d/'SOlapm.dat').exists():raise RuntimeError('Unexpected saved normalization: this test computes it afresh every cycle.')
  mu=float((d/'EF.dat').read_text().split()[0]);p=d/(case+'.indmf2');ls=p.read_text().splitlines();parts=ls[-1].split();parts[1]=str(mu);parts[2]='1';p.write_text(ls[0]+'\n'+' '.join(parts)+'\n')
  command('dmft2',['/opt/homebrew/bin/mpirun','-n','4','/Users/haule/dbin/dmft2','dmft2.def'])
  command('lcore',['/Users/haule/wien2k_24/x','lcore','-f',case])
  command('merge',['/Users/haule/dbin/merge_scf.py','--replace'])
  command('mixer',['/Users/haule/wien2k_24/x','mixer','-f',case])
  en=tag(d/(case+'.scfm'),':ENE  :');su=tag(d/(case+'.scf2'),':SUM  :');vs=tag(d/(case+'.scf2'),':VSUM :');K=en-su+vs-tag(d/'Eorb_imp.dat',':EORB')
  s=(d/(case+'.scfm')).read_text();dis=float(re.findall(r'^:DIS.*?([0-9.]+)\s*$',s,re.M)[-1]);dk=abs(K-history[-1]['K']) if history else None
  history.append(dict(iteration=it,ENE=en,SUM=su,VSUM=vs,K=K,DIS=dis,deltaK=dk,mu=float((d/'EF.dat').read_text().split()[0]),seconds=time.monotonic()-t))
  (d/'history.json').write_text(json.dumps(history,indent=2)+'\n')
  print(name,it,'K',f'{K:.12f}','DIS',dis,'dK',dk,flush=True)
  stable=stable+1 if dis<=1e-7 and dk is not None and dk<=1.1e-8 else 0
  if stable>=2:break
 result={'directory':str(d),'strain_nominal':h,'strain_actual':float(ls[0].split()[0]) if False else float((d/(case+'.struct')).read_text().splitlines()[3][:10])/7.260328-1,'history':history,'converged':stable>=2}
 (d/'results.json').write_text(json.dumps(result,indent=2)+'\n');return name,result
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:x=dict(pool.map(run,[('half_minus','V_strain_m001',-.0005),('half_plus','V_strain_p001',.0005)]))
fd=(x['half_plus']['history'][-1]['K']-x['half_minus']['history'][-1]['K'])/(3*382.70904*(x['half_plus']['strain_actual']-x['half_minus']['strain_actual']))*14710.507641078
x['fd_GPa']=fd;(root/'halfstep_results.json').write_text(json.dumps(x,indent=2)+'\n');print('HALFSTEP FD',fd,flush=True)
