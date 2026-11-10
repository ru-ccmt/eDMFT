from pathlib import Path
import os,shutil,subprocess,re,json
root=Path(__file__).resolve().parent;case='0_800_MnTe';W=Path('/Users/haule/wien2k_24')
env=dict(os.environ,WIENROOT=str(W),OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1');env['PATH']=str(W)+':/Users/haule/dbin:/opt/homebrew/bin:'+env['PATH']
def num(p,key):return float(re.findall('^'+re.escape(key)+r'\s*(.*)$',p.read_text(),re.M)[-1].split('=')[-1].split()[0])
out=[]
for name in ['0_800_MnTe_MnTe_m_0.001','0_800_MnTe_0_800_MnTe','0_800_MnTe_MnTe_p_0.001']:
 src=root/'DFT_runs'/name;dst=root/'DFT_TEMPS'/name
 # An identical repeat density makes the reused MSR1 history singular.
 # Reset only the isolated scratch history, preserving the files for review.
 if ':ENE  :' not in (dst/(case+'.scfm')).read_text():
  archive=dst/'repeated_history';archive.mkdir(exist_ok=True)
  for p in dst.glob(case+'.broyd*'):p.rename(archive/p.name)
  for p in [dst/'mixer_temps.log',dst/(case+'.scfm')]:shutil.copy2(p,archive/p.name)
  shutil.copy2(src/(case+'.clmsum_old'),dst/(case+'.clmsum'))
  for label,cmd in [('merge',['/Users/haule/miniconda3-clean/envs/work311/bin/python','/Users/haule/dbin/merge_scf.py','-f',case,'--replace']),('mixer',[str(W/'x'),'mixer','-f',case])]:
   with (dst/(label+'_temps_restart.log')).open('w') as log:subprocess.run(cmd,cwd=dst,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 old=num(src/(case+'.scfm'),':ENE  :');new=num(dst/(case+'.scfm'),':ENE  :')
 correction=float(re.search(r'-\(T\*S\)/2\s*=\s*([-\d.]+)',(src/(case+'.scf2')).read_text()).group(1))
 same=(src/(case+'.clmval')).read_bytes()==(dst/(case+'.clmval')).read_bytes();assert same
 assert abs(new-old-correction)<3e-8,(name,new-old,correction)
 ans={'name':name,'energy_TEMP_Ry':old,'energy_TEMPS_Ry':new,'extra_entropy_Ry':correction,'charge_identical':same,'scratch_history_restarted':(dst/'repeated_history').exists()};out.append(ans);print(ans)
(root/'temps_check.json').write_text(json.dumps(out,indent=2)+'\n')
