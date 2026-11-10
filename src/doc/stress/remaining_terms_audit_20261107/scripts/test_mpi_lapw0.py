from pathlib import Path
import subprocess,os,shutil,re,numpy as np,json
r=Path(__file__).resolve().parent;results={}
for label,nproc in [('SrVO3_DMFT',4),('FeSe_DMFT',2)]:
 src=r/label;dst=r/'mpi_runs'/label;dst.mkdir(parents=True,exist_ok=True);case=label.split('_')[0]
 for p in src.iterdir():
  if p.is_file() and (p.name.startswith(case+'.') or p.name=='lapw0.def'):shutil.copy2(p,dst/p.name)
 with (dst/'run.log').open('w') as f:subprocess.run(['/opt/homebrew/bin/mpirun','-n',str(nproc),str(r/'SRC_lapw0_mpi/lapw0_mpi'),'lapw0.def'],cwd=dst,env=os.environ|{'OMP_NUM_THREADS':'1','WIEN_STRESS_DEBUG':'1'},stdout=f,stderr=subprocess.STDOUT,check=True)
 def tensor(p):
  s=p.read_text();return np.array([list(map(float,re.findall('^:STR_DENS'+f'{i:03d}:'+r'\s*(.*)',s,re.M)[-1].split())) for i in [1,2,3]])
 delta=tensor(dst/(case+'.scf0'))-tensor(src/'stage'/(case+'.scf0'))
 results[label]={'ranks':nproc,'max_difference_GPa':float(np.max(abs(delta))*14710.5164057)};print(label,results[label],flush=True)
 assert results[label]['max_difference_GPa']<1e-6
(r/'mpi_results.json').write_text(json.dumps(results,indent=2)+'\n')
