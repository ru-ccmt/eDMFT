"""Independent full native scalar along a fixed density tangent at small strain."""
from pathlib import Path
import sys,os,json,re,time,shutil,subprocess
import numpy as np
root=Path(__file__).resolve().parent;sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_potential_core_state import load_potential_core_state
from lapw_density_functional import PreparedDensityFunctional
from lapw_functional_cancellation import CoreEigenvalueFunctional
from lapw_strain_field import diagonal_strain_potential
from wien_response_density import Density
case='FeSe';src=root/case/'refined';pc=load_potential_core_state(src/'potential_core.npz');q=np.load(root/case/'halfstep/density_tangent.npy');s=pc.structure;h=float(sys.argv[1]) if len(sys.argv)>1 else .000125;tag='fresh_direction' if h==.000125 else 'fresh_direction_'+str(h)
obs=PreparedDensityFunctional(pc);core=CoreEigenvalueFunctional.from_state(pc);ve,_=diagonal_strain_potential(pc,np.eye(3));v=ve+pc.potential_operator.matvec(q)
analytic={'DEN':float(obs.explicit(np.eye(3))[0]+obs.action(q)),'core':float(core.action(v))}
# Complete derivative already saved using this same native density tangent.
analytic['total']=json.loads((root/case/'halfstep/density_readout.json').read_text())['derivative_with_native_density_tangent_Ry'];analytic['lattice']=analytic['total']-analytic['core']-analytic['DEN']
weights=np.array([s.mult[a]*st.occ for a,atom in pc.core.atoms.items() for st in atom.states]);shift=sum(s.mult[a]*sum(st.occ for st in atom.states)*pc.metadata['core'][a]['shift_Ry'] for a,atom in pc.core.atoms.items());rows=[]
env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))};env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',WIEN_STRESS_DEBUG='1')
for label,sign in [('minus',-1),('plus',1)]:
 out=root/case/tag/label;out.mkdir(parents=True,exist_ok=True)
 if not (out/'RESULT.json').exists():
  names=[case+'.'+z for z in ['struct','in0','in1','in2','inm','inc','indmfl','indmfi','indmf2','klist','kgen','rotlm']]+['EF.dat','Eorb_imp.dat','projectorw.dat','lapw0.def','lapw1.def','dmft2.def','lcore.def']+[p.name for p in src.glob('sig.inp*')]
  for name in names:
   if (src/name).exists():shutil.copy2(src/name,out/name)
  lines=(src/(case+'.struct')).read_text().splitlines();cell=np.array([s.a,s.b,s.c]);lines[3]=''.join(f'{x*(1+sign*h):10.6f}' for x in cell)+lines[3][30:];new=np.array([float(lines[3][i:i+10]) for i in [0,10,20]]);scale=float(np.mean(new/cell));(out/(case+'.struct')).write_text('\n'.join(lines)+'\n')
  rho=pc.input_density.from_vector(pc.input_density.vector()+(scale-1)*q).transport(scale);template=Density(src/(case+'.clmsum'),src/(case+'.struct'));template.rad={k:rho.rad[k].copy() for k in template.rad};template.pw={k:rho.pw[k] for k in template.pw};template.write(out/(case+'.clmsum'))
  p=out/(case+'.in1');lines=p.read_text().splitlines();w=lines[1].split();w[0]=f'{float(w[0])/scale:.8f}';lines[1]=' '.join(w);p.write_text('\n'.join(lines)+'\n')
  p=out/(case+'.in2');lines=p.read_text().splitlines()
  for i,line in enumerate(lines):
   if 'GMAX' in line:lines[i]=f'{float(line.split()[0])/scale:.8f} GMAX'
  p.write_text('\n'.join(lines)+'\n')
  p=out/(case+'.indmf2');lines=p.read_text().splitlines();i=next(i for i,l in enumerate(lines) if l.strip() and not l.lstrip().startswith('#'));w=lines[i].split();w[2]='1';lines[i]=' '.join(w);p.write_text('\n'.join(lines)+'\n')
  shutil.copy2('/private/tmp/edmft_screened_native_20261108/full_bz/FeSe/central/lapw1.fixgvec',out/'lapw1.fixgvec');(out/'lapw1.fixgvec.read').touch()
  workers={'lapw0':'/private/tmp/edmft_xc_core_response_20261108/native/SRC_lapw0/lapw0','lapw1':'/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1/lapw1','dmft2':'/Users/haule/dbin/dmft2_response','lcore':'/private/tmp/edmft_xc_core_response_20261108/native/SRC_lcore/lcore'}
  for name,exe in workers.items():
   with (out/(name+'.log')).open('w') as f:subprocess.run([exe,name+'.def'],cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
  den=float(re.findall(r'^:DBG_EN_DEN_TOTAL:\s*(\S+)',(out/(case+'.scf0')).read_text(),re.M)[-1]);L=float(re.findall(r'^:VSUM.*?=\s*(\S+)',(out/(case+'.scf2')).read_text(),re.M)[-1]);es=[]
  with (out/'core_native.bin').open('rb') as f:
   while True:
    head=np.fromfile(f,'i4',7)
    if not len(head):break
    es.append(float(np.fromfile(f,'f8',4)[0]));f.seek(int(4*head[3]*8),1)
  C=float(2*weights@np.array(es)+shift);data=dict(offset=scale-1,DEN=den,core=C,lattice=L,total=den+C+L);(out/'RESULT.json').write_text(json.dumps(data,indent=2)+'\n')
 rows.append(json.loads((out/'RESULT.json').read_text()));print(label,rows[-1],flush=True)
fd={key:(rows[1][key]-rows[0][key])/(rows[1]['offset']-rows[0]['offset']) for key in analytic};report=dict(analytic=analytic,native_fd=fd,errors={key:fd[key]-analytic[key] for key in fd},rows=rows);(root/case/tag/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
