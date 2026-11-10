from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import os,subprocess,shutil,json,re,hashlib
root=Path(__file__).resolve().parent;case='0_800_MnTe';W=Path('/Users/haule/wien2k_24')
env=dict(os.environ,WIENROOT=str(W),OMP_NUM_THREADS='2',VECLIB_MAXIMUM_THREADS='2',OPENBLAS_NUM_THREADS='1');env['PATH']=str(W)+':/Users/haule/dbin:/opt/homebrew/bin:'+env['PATH']
for k in list(env):
 if k.startswith(('WIEN_STRESS_','EDMFT_STRESS_')):del env[k]
def number(p,key):return float(re.findall('^'+re.escape(key)+r'\s*(.*)$',p.read_text(),re.M)[-1].split('=')[-1].split()[0])
def one(name):
 src=root/'DFT_runs'/name;dst=root/'DFT_TEMPS'/name;shutil.copytree(src,dst)
 shutil.copy2(dst/(case+'.clmsum_old'),dst/(case+'.clmsum'))
 f=dst/(case+'.in2');s=f.read_text();assert re.search(r'^TEMP\s',s,re.M);f.write_text(re.sub(r'^TEMP\s+', 'TEMPS    ',s,count=1,flags=re.M))
 for label,cmd in [('lapw2',[str(W/'x'),'lapw2','-f',case]),('merge',['/Users/haule/miniconda3-clean/envs/work311/bin/python','/Users/haule/dbin/merge_scf.py','-f',case,'--replace']),('mixer',[str(W/'x'),'mixer','-f',case])]:
  with (dst/(label+'_temps.log')).open('w') as log:subprocess.run(cmd,cwd=dst,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 old=number(src/(case+'.scfm'),':ENE  :');new=number(dst/(case+'.scfm'),':ENE  :')
 correction=float(re.search(r'-\(T\*S\)/2\s*=\s*([-\d.]+)',(src/(case+'.scf2')).read_text()).group(1))
 same=(src/(case+'.clmval')).read_bytes()==(dst/(case+'.clmval')).read_bytes();assert same
 assert abs(new-old-correction)<3e-8,(name,old,new,correction)
 ans={'name':name,'energy_TEMP_Ry':old,'energy_TEMPS_Ry':new,'extra_entropy_Ry':correction,'charge_identical':same}
 print(ans,flush=True);return ans
names=['0_800_MnTe_MnTe_m_0.001','0_800_MnTe_0_800_MnTe','0_800_MnTe_MnTe_p_0.001']
with ThreadPoolExecutor(max_workers=2) as ex:ans=list(ex.map(one,names))
(root/'temps_check.json').write_text(json.dumps(ans,indent=2)+'\n')
