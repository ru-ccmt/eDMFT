from pathlib import Path
import re,subprocess,os,json,shutil,numpy as np
root=Path(__file__).resolve().parent
results={}
for label in ['SrVO3_DFT','SrVO3_DMFT','FeSe_DFT','FeSe_DMFT']:
 src=root/label;case=label.split('_')[0];dst=root/'xc_scaling'/label;dst.mkdir(parents=True,exist_ok=True)
 for p in src.iterdir():
  if p.is_file() and (p.name.startswith(case+'.') or p.name=='lapw0.def'):shutil.copy2(p,dst/p.name)
 original=(src/(case+'.struct')).read_text().splitlines();rho=(src/(case+'.clmsum')).read_text().splitlines()
 npt=[int(t[15:20]) for t in original if 'NPT=' in t];pw=next(i for i,s in enumerate(rho) if 'NUMBER OF PW' in s);npw=int(rho[pw].split()[0])
 p=dst/(case+'.in0');text=p.read_text();p.write_text('KXC '+text[4:]);excs={};ees={}
 for h in [-.004,-.002,.002,.004]:
  scale=1+h;st=original.copy();abc=[float(st[3][10*i:10*i+10]) for i in range(6)];st[3]=''.join(f'{v*(scale if i<3 else 1):10.6f}' for i,v in enumerate(abc))
  for i,t in enumerate(st):
   if 'NPT=' in t:st[i]=t[:25]+f'{float(t[25:35])*scale:.9f}'[1:]+t[35:40]+f'{float(t[40:50])*scale:10.5f}'+t[50:]
  (dst/(case+'.struct')).write_text('\n'.join(st)+'\n');rr=rho.copy();atom=0
  for i,t in enumerate(rho[:pw]):
   if 'ATOMNUMBER' in t:atom=int(t.split('=')[-1])-1
   if 'CLM(R)' in t:
    k=i+2;left=npt[atom]
    while left:
     count=min(4,left);v=[float(rr[k][3+19*j:3+19*(j+1)]) for j in range(count)];rr[k]='   '+''.join(f'{x/scale:19.12E}' for x in v);k+=1;left-=count
  for i in range(pw+1,pw+npw+1):
   t=rho[i];rr[i]=t[:18]+''.join(f'{float(t[18+19*j:18+19*(j+1)])/scale**3:19.12E}' for j in range(2))
  (dst/(case+'.clmsum')).write_text('\n'.join(rr)+'\n')
  with (dst/f'lapw0_{h}.log').open('w') as f:subprocess.run([str(root/'SRC_lapw0/lapw0'),'lapw0.def'],cwd=dst,env=os.environ|{'OMP_NUM_THREADS':'1'},stdout=f,stderr=subprocess.STDOUT,check=True)
  s=(dst/(case+'.scf0')).read_text();(dst/f'scf0_{h}').write_text(s)
  def val(tag):return float(re.findall('^:'+tag+r'\s*:\s*.*=\s*([-+0-9.Ee]+)',s,re.M)[-1])
  excs[h]=val('EXC');ees[h]=val('DEN')-val('INTRV')-val('EXC')
 stress=(src/'stage'/(case+'.scf0')).read_text();rows=[s.split() for s in re.findall(r'^:DBG_STR_DENS:\s*(.*)',stress,re.M)];xcstress=sum(sum(map(float,t[5:7])) for t in rows if t[0]==t[1]);vol=float(re.findall(r'^:VOL\s*:\s*.*?([0-9]+\.[0-9]+)',stress,re.M)[-1]);fd=(excs[-.004]-8*excs[-.002]+8*excs[.002]-excs[.004])/.024
 fd2=(excs[.002]-excs[-.002])/.004
 results[label]={'xc_fd_Ry':fd,'xc_fd_2point_Ry':fd2,'xc_stress_trace_Ry':xcstress*vol,'mean_error_GPa':(xcstress-fd/vol)*14710.5164057/3,'Exc':excs,'Ees':ees}
 print(label,results[label],flush=True)
(root/'xc_scaling_results.json').write_text(json.dumps(results,indent=2)+'\n')
