from material_augmentation import *
import os,shutil,time
case=sys.argv[1];mat=Material(case);s=mat.s
root=R/case/'moving_probe';root.mkdir(exist_ok=True)
dv={}
for ia in range(s.nat):
 n=s.jrj[ia];r=s.r0[ia]*np.exp(np.log(s.rmt[ia]/s.r0[ia])/(n-1)*np.arange(n));dv[ia]=.07*r*np.sin(1.3*r)
 for l in range(5):mat.radial(ia,l,dv[ia])
# Actual fixed on-disk projector 5. Keep its radial MT boundary value finite.
path=Path('/Users/haule/projects/materials/SrVO3/V_plus_force/projectorw.dat') if case=='SrVO3' else Path('/Users/haule/projects/materials/FeSe_again/FeSe_stress/FeSe_RELA/projectorw.dat')
projectors={}
with path.open() as f:
 count=int(f.readline().split()[1])
 for _ in range(count):
  nr,nr0,ia,l=map(int,f.readline().split()[1:]);rows=np.array([list(map(float,f.readline().split())) for _ in range(nr)])
  if l==2:
   if nr0!=s.jrj[ia-1]:raise ValueError('projector grid length')
   np.testing.assert_allclose(rows[:nr0,0],mat.context[ia-1,2]['r'],rtol=1e-8,atol=1e-9)
   projectors[ia-1]=rows[:nr0,1:3]
print('prepared',case,'projector atoms',list(projectors),flush=True)
# One central gauge is sufficient for testing all native angular sectors.
ik=0;ev,c,labels=mat.cv[ik];analytic_h=np.zeros((len(ev),len(ev)),complex);analytic_o=analytic_h.copy();p0=[];dp0=[]
layout_records=[];latom=0
for ia in range(s.nat):
 layout,k,overlap,proj,dk,do,dp=mat.small_kernels(ia,mat.vns[ia],projectors.get(ia))
 for mu in range(s.mult[ia]):
  z,dz,*indices=mat.coefficients(ia,latom,ik)
  aug=AugmentationContractions(z,k,overlap,proj)
  dh,ds,dphi=aug.action(dz,dk,do,dp);analytic_h+=dh;analytic_o+=ds
  if ia in projectors:p0.append(aug.value()[2]);dp0.append(dphi)
  layout_records.append((ia,latom,aug));latom+=1
p0=np.array(p0);dp0=np.array(dp0)
env=os.environ.copy()
for key in list(env):
 if key.startswith(('WIEN_','EDMFT_','DFTKS_')):del env[key]
env.update(OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
exe=R/'native/SRC_lapw1'/('lapw1c' if mat.complex else 'lapw1')

def perturb_file(text,amplitude,basis):
 patches=[]
 for block in re.finditer(r'ATOMNUMBER\s*=\s*(\d+)',text):
  ia=int(block[1])-1;end=text.find('ATOMNUMBER',block.end());end=len(text) if end<0 else end
  section=text[block.end():end]
  for channel in re.finditer(r'VLM\(R\) FOR L\s*(-?\d+)\s+M=\s*(\d+)',section):
   offset=block.end()+channel.end()
   for i,match in enumerate(list(FLOAT.finditer(text[offset:end]))[:s.jrj[ia]]):
    value=float(match.group().replace('D','E'))+amplitude*dv[ia][i] if basis else 0.
    patches.append((offset+match.start(),offset+match.end(),value))
 for a,b,value in reversed(patches):text=text[:a]+f'{value:{b-a}.12E}'+text[b:]
 return text

def run(step,sign,zero):
 d=root/(('zero_' if zero else 'full_')+('m' if sign<0 else 'p')+f'_{step:.0e}');d.mkdir(exist_ok=True)
 for p in mat.source.iterdir():
  if p.is_file() and (p.suffix in {'.struct','.vsp','.vns','.in1','.in1c','.klist'} or p.name=='lapw1.def'):shutil.copy2(p,d/p.name)
 p=d/(case+'.klist');p.write_text(p.read_text().splitlines()[0]+'\nEND\n')
 p=d/'basis.vsp';p.write_text(perturb_file((d/(case+'.vsp')).read_text(),sign*step,True))
 if zero:
  p=d/(case+'.vns');p.write_text(perturb_file(p.read_text(),0,False))
 if not (d/'done').exists():
  with (d/'run.log').open('w') as log:proc=subprocess.run([str(exe),'lapw1.def'],cwd=d,env=env|{'WIEN_LAPW1_MATRIX_DUMP':str(d/'native'),'WIEN_LAPW1_BASIS_VSP':str(d/'basis.vsp')},stdout=log,stderr=subprocess.STDOUT)
  proc.check_returncode();(d/'done').touch()
 with FortranFile(d/(case+'.vector')) as f:
  for ia in range(s.nat):
   np.testing.assert_allclose(f.read_record('f8'),mat.header[ia],atol=1e-13,rtol=0)
   np.testing.assert_allclose(f.read_record('f8').reshape(-1,4).T,mat.lo[ia],atol=1e-13,rtol=0)
 h=read_basis_matrix(d/'native.ham.k00000001.bin',b'EDMFTHV1');o=read_basis_matrix(d/'native.k00000001.bin',b'EDMFTOV1')
 np.testing.assert_array_equal(h.labels,labels);np.testing.assert_array_equal(o.labels,labels)
 return h.matrix,o.matrix

rows=[]
for step in [2e-3,1e-3,5e-4]:
 hs=[];os=[]
 for sign in [-1,1]:
  full=run(step,sign,False);zero=run(step,sign,True);hs.append(full[0]-zero[0]);os.append(full[1])
 fd=c.conj().T@((hs[1]-hs[0])/(2*step))@c
 # Overlap here is restricted to l<=4 and cannot be compared to the all-l
 # native overlap derivative. Full spherical assembly is a separate check.
 row=dict(step=step,H_nonspherical_relative=float(np.linalg.norm(fd-analytic_h)/np.linalg.norm(fd)),H_nonspherical_max=float(np.max(abs(fd-analytic_h))))
 pend=[]
 for sign in [-1,1]:
  for ia in range(s.nat):
   mat.v[ia][0,0]+=sign*step*dv[ia]
   for l in range(5):mat.radial(ia,l)
   mat.v[ia][0,0]-=sign*step*dv[ia]
  pvalue=[]
  for ia,latom,_ in layout_records:
   if ia not in projectors:continue
   layout,k,overlap,proj,*_=mat.small_kernels(ia,mat.vns[ia],projectors[ia]);z,*_=mat.coefficients(ia,latom,ik)
   pvalue.append(z.conj().T@proj)
  pend.append(np.array(pvalue))
 fdp=(pend[1]-pend[0])/(2*step)
 row.update(raw_projector_relative=float(np.linalg.norm(fdp-dp0)/np.linalg.norm(dp0)),raw_projector_max=float(np.max(abs(fdp-dp0))))
 rows.append(row);print(case,row,flush=True)
(root/'report.json').write_text(json.dumps(rows,indent=2)+'\n')
np.savez(root/'central_response.npz',h_nonspherical=analytic_h,overlap_l4=analytic_o,raw_projector=p0,raw_projector_dot=dp0)
