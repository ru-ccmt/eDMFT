from pathlib import Path
import subprocess,os,json,re
import numpy as np
root=Path(__file__).resolve().parent
report=json.loads((root/'material_comparison.json').read_text())
env=dict(os.environ,WIENROOT='/Users/haule/wien2k_24',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
env['PATH']='/Users/haule/wien2k_24:/Users/haule/dbin:/opt/homebrew/bin:'+env['PATH']
out={}
for case in ['SrVO3','FeSe']:
 folder=root/'runs'/('DFT_'+case)
 for label,cmd in [('merge',['/Users/haule/miniconda3-clean/envs/work311/bin/python','/Users/haule/dbin/merge_scf.py','-f',case,'--replace']),('mixer',['/Users/haule/wien2k_24/x','mixer','-f',case])]:
  with (folder/(label+'_g1.log')).open('w') as f:subprocess.run(cmd,cwd=folder,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 s=(folder/(case+'.scfm')).read_text()
 t=np.array([[float(x) for x in re.findall(r'^:STRESS_GPa'+f'{i:03d}:'+r'\s*(.*)',s,re.M)[-1].split()[:3]] for i in [1,2,3]])
 err=float(abs(t-np.array(report['DFT_'+case]['tensor_gpa'])).max());assert err<1e-5
 out[case]={'mixer_tensor_gpa':t.tolist(),'max_rounding_error_gpa':err}
(root/'mixer_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)
