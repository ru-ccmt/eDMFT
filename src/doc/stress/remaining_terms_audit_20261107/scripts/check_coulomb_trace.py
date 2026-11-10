from pathlib import Path
import os,subprocess,re,json,numpy as np
root=Path(__file__).resolve().parent;results={}
for label in ['SrVO3_DFT','SrVO3_DMFT','FeSe_DFT','FeSe_DMFT']:
 d=root/label;case=label.split('_')[0];p=d/(case+'.in0');text=p.read_text()
 try:
  p.write_text('KXC '+text[4:])
  with (d/'kxc.log').open('w') as f: subprocess.run([str(root/'SRC_lapw0/lapw0'),'lapw0.def'],cwd=d,env=os.environ|{'WIEN_STRESS_DEBUG':'1','OMP_NUM_THREADS':'1'},stdout=f,stderr=subprocess.STDOUT,check=True)
  s=(d/(case+'.scf0')).read_text();(d/'kxc.scf0').write_text(s)
 finally:p.write_text(text)
 def val(tag):return float(re.findall('^:'+tag+r'\s*:\s*.*=\s*([-+0-9.Ee]+)',s,re.M)[-1])
 den=val('DEN');ekin=val('INTRV');exc=val('EXC');ees=den-ekin-exc
 stress=(d/'stage'/(case+'.scf0')).read_text();x=[]
 for row in re.findall(r'^:DBG_STR_DENS:\s*(.*)',stress,re.M):
  t=row.split()
  if t[0]==t[1]:x.append(sum(map(float,t[2:5])))
 vol=float(re.findall(r'^:VOL\s*:\s*.*?([0-9]+\.[0-9]+)',stress,re.M)[-1])
 results[label]={'E_es_Ry':ees,'density_coulomb_trace_Ry':sum(x)*vol,'identity_residual_mean_GPa':(sum(x)+ees/vol)*14710.5164057/3,'E_xc_Ry':exc,'volume':vol}
 print(label,results[label],flush=True)
(root/'coulomb_trace_results.json').write_text(json.dumps(results,indent=2)+'\n')
