from pathlib import Path
import os,shutil,subprocess,time,json,re,concurrent.futures,hashlib
root=Path('/private/tmp/numeric_stress_audit_20261107');case='SrVO3';parent=Path('/Users/haule/projects/materials/SrVO3')
env=os.environ|{'WIENROOT':'/Users/haule/wien2k_24','WIEN_DMFT_ROOT':'/Users/haule/dbin','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'}
for k in list(env):
 if k.startswith(('EDMFT_STRESS_','WIEN_STRESS_','EDMFT_PROJECTOR_')):del env[k]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def tag(p,t):return float(re.findall(r'^'+re.escape(t)+r'\s*(.*)$',p.read_text(),re.M)[-1].split('=')[-1].split()[0])
def run(spec):
 name,srcname=spec;src=parent/srcname;d=root/name;d.mkdir(exist_ok=True)
 sufs=['struct','clmsum','clmup','clmdn','in0','in1','in1c','in2','in2c','inc','inm','klist','kgen','rotlm','vsp','vns','clmcor','r2v','indmfl','indmfi','indmf2']
 hashes={}
 for filename in [case+'.'+s for s in sufs]+['EF.dat','Eorb_imp.dat','projectorw.dat','SOlapm.dat','sig.inp','sig.inp1']:
  p=src/filename
  if p.exists():shutil.copy2(p,d/filename);hashes[filename]=sha(p)
 (d/'mpi_prefix.dat').write_text('/opt/homebrew/bin/mpirun -n 4\n')
 def command(label,cmd):
  print('START',name,label,flush=True);start=time.monotonic()
  with (d/(label+'.log')).open('w') as log:r=subprocess.run(cmd,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
  print('END',name,label,r.returncode,round(time.monotonic()-start,2),flush=True)
  if r.returncode:raise RuntimeError((name,label,r.returncode))
 def capture(label):
  out=d/label;out.mkdir(exist_ok=True)
  for n in [case+'.scf0',case+'.scfc',case+'.scf1',case+'.scf2',case+'.scfm',case+'.scf',case+'.indmf2','SOlapm.dat','EF.dat','sig.inp1']:
   if (d/n).exists():shutil.copy2(d/n,out/n)
  a={t:tag(d/(case+'.scf2'),t) for t in [':SUM  :',':VSUM :',':YSUM :',':DRHO']}
  if (d/(case+'.scfm')).exists():
   a['ENE']=tag(d/(case+'.scfm'),':ENE  :');a['K']=a['ENE']-a[':SUM  :']+a[':VSUM :']-tag(d/'Eorb_imp.dat',':EORB')
  a['SOlapm_hash']=sha(d/'SOlapm.dat');a['EF_eV']=float((d/'EF.dat').read_text().split()[0]);return a
 command('lapw0',['/Users/haule/wien2k_24/x','lapw0','-f',case])
 command('lapw1',['/Users/haule/dbin/x_dmft.py','-p','lapw1'])
 command('prepare_dmft2',['/Users/haule/dbin/x_dmft.py','-d','-p','dmft2'])
 shutil.copy2(src/(case+'.indmf2'),d/(case+'.indmf2'))
 for line in (d/'dmft2.def').read_text().splitlines():
  n=line.split("'")[1];assert not Path(n).is_absolute(),n
 command('dmft2_baseline',['/opt/homebrew/bin/mpirun','-n','4','/Users/haule/dbin/dmft2','dmft2.def'])
 command('lcore',['/Users/haule/wien2k_24/x','lcore','-f',case])
 command('merge',['/Users/haule/dbin/merge_scf.py','--replace'])
 command('mixer',['/Users/haule/wien2k_24/x','mixer','-f',case])
 result={'source':str(src),'input_hashes':hashes,'baseline':capture('baseline')}
 # Keep density-energy contribution from baseline exactly fixed for the mesh diagnostic.
 den_energy=result['baseline']['ENE']-result['baseline'][':SUM  :']-tag(d/'Eorb_imp.dat',':EORB')
 command('split_all_input_frequencies',['/Users/haule/dbin/ssplit.py','-n','100000','-t','0'])
 command('dmft2_fullmesh',['/opt/homebrew/bin/mpirun','-n','4','/Users/haule/dbin/dmft2','dmft2.def'])
 result['fullmesh']=capture('fullmesh');result['fullmesh']['K_fixed_density']=den_energy+result['fullmesh'][':VSUM :']
 result['source_inputs_preserved']=all(sha(src/n)==h for n,h in hashes.items())
 (d/'results.json').write_text(json.dumps(result,indent=2)+'\n')
 return name,result
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:result=dict(pool.map(run,[('minus','V_strain_m001'),('plus','V_strain_p001')]))
(root/'replay_results.json').write_text(json.dumps(result,indent=2)+'\n')
fac=14710.507641078/(.006*382.70904)
for mode,key in [('baseline','K'),('fullmesh','K_fixed_density')]:print('FD',mode,(result['plus'][mode][key]-result['minus'][mode][key])*fac,flush=True)
