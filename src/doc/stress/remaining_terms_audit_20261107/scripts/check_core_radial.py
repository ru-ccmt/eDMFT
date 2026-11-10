from pathlib import Path
import numpy as np,re,json
from scipy.interpolate import CubicSpline
from scipy.integrate import simpson
root=Path('/private/tmp/stress_remaining_audit_20261107');out={}
for d in root.iterdir():
 if not d.is_dir() or not (d/'audit_lcore.log').exists():continue
 case='SrVO3' if d.name.startswith('SrVO3') else 'FeSe'
 struct=(d/(case+'.struct')).read_text()
 jm=[int(x) for x in re.findall(r'NPT=\s*(\d+)',struct)];mult=[int(x) for x in re.findall(r'MULT=\s*(\d+)',struct)]
 vblocks=re.split(r'# atom\s+\d+\s*\n',(d/(case+'.lcorepot')).read_text())[1:]
 qblocks=re.split(r'ATOMNUMBER\s*=\s*\d+', (d/(case+'.rsplcore')).read_text())[1:]
 lines=(d/(case+'.outputc')).read_text().splitlines();stored={}
 for l in lines:
  if l.startswith(':DBG_STR_CORE_ATOM:'):
   x=l.split();a,i=int(x[1]),int(x[2]);
   if i==1:stored[a]=float(x[3])
 ans=[]
 for i,(vblock,qblock) in enumerate(zip(vblocks,qblocks)):
  rv=np.array([[float(x) for x in line.split()] for line in vblock.splitlines() if len(line.split())==2]);r=rv[:,0];v=2*rv[:,1]
  qtext=qblock.split('RLM(R) FOR L= 0   M= 0')[1]
  q=np.array([float(x) for x in re.findall(r'[-+]?\d+\.\d+(?:E[-+]\d+)',qtext)])
  assert len(q)==len(r),(i,len(q),len(r))
  n=jm[i]; rmt=r[n-1];rr=r[:n];vv=v[:n];qq=q[:n];z=float(qblock.splitlines()[0].split()[-1]);x=np.log(rr)
  # Spline the regular remainder in logarithmic radius; add nuclear pole exactly.
  smooth=vv+2*z/rr
  dv_dx=CubicSpline(x,smooth)(x,1)+2*z/rr
  exactlog=simpson(qq*dv_dx*rr,x=x)*mult[i]/3
  # Unchanged integrator endpoint stencil applied to the same derivative.
  f=qq*dv_dx*rr;weights=np.ones(n);weights[:4]=[17/48,59/48,43/48,49/48];weights[-4:]=weights[:4][::-1]
  samequad=float(np.dot(weights,f)*(x[1]-x[0]))*mult[i]/3
  fullvir=simpson(q*(CubicSpline(np.log(r),v+2*z/r)(np.log(r),1)+2*z/r)*r,x=np.log(r))*mult[i]/3
  ans.append({'atom':i+1,'stored_Ry':stored[i+1],'log_regularized_Ry':samequad,'difference_Ry':samequad-stored[i+1],'quad_sensitivity_Ry':exactlog-samequad,'core_charge_outside_MT':simpson(q[n-1:]*r[n-1:],x=np.log(r[n-1:])), 'q_at_MT':qq[-1],'tail_virial_Ry':fullvir-exactlog})
 vol=float(re.findall(r':VOL.*=\s*([\d.]+)',(d/(case+'.scf')).read_text())[-1]);err=sum(a['difference_Ry'] for a in ans)*14710.5164057/vol
 out[d.name]={'atoms':ans,'mean_stress_difference_GPa':err}
 print(d.name,'independent derivative minus stored, GPa',err,'outside charges',[a['core_charge_outside_MT'] for a in ans])
(root/'core_radial_results.json').write_text(json.dumps(out,indent=2)+'\n')
