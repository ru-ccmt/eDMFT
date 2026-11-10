import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys,struct,re,json,time
sys.dont_write_bytecode=True
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from pathlib import Path
import numpy as np
from scipy.io import FortranFile,FortranEOFError
from scipy.linalg import solve
from scipy.special import spherical_jn
from check_projector_sternheimer_fd import read_dump,read_normalization,read_response_matrix
from check_interstitial_overlap_stress import read_atoms
from wstruct import Struct
root=Path(__file__).parent

def vectors(p,nat):
 f=FortranFile(p);rad=[f.read_record('f8') for j in range(2*nat)]
 while True:
  try:head=f.read_record('u1').tobytes()
  except FortranEOFError:return
  x,y,z,name,n,nb,w=struct.unpack('<ddd10siid',head)
  labels=f.read_record('i4').reshape(3,n,order='F');e=[];v=[]
  for j in range(nb):
   index,ej=struct.unpack('<id',f.read_record('u1').tobytes());e.append(ej)
   vv=f.read_record('u1').tobytes();dtype='f8' if len(vv)==8*n else 'c16';v.append(np.frombuffer(vv,dtype=dtype))
  yield np.array(e),np.array(v).T,labels,rad

def tensor(p,tag):
 text=p.read_text();return np.array([[float(x) for x in re.findall('^'+re.escape(tag)+f'{i:03d}'+r':?\s+([^\n]+)',text,re.M)[-1].split()[:3]] for i in range(1,4)])

def run(case):
 r=root/case;c=r/'center';s=Struct();s.ReadStruct(str(c/(case+'.struct')));atoms=read_atoms(c/(case+'.struct'))
 cell=np.array([s.a,s.b,s.c]);vol=np.prod(cell)/(4 if s.lattyp=='F' else 1);fac=14710.507641078/(3*vol)
 ffttext=(r/'r2v_center'/(case+'.output0')).read_text();shape=tuple(map(int,re.search(r'global FFT-parameters:\s*(\d+)\s+(\d+)\s+(\d+)',ffttext).groups()))
 xyz=np.array(np.meshgrid(*[np.fft.fftfreq(n)*n for n in shape],indexing='ij'));qnorm=np.sqrt(np.sum((xyz*(2*np.pi/cell)[:,None,None,None])**2,axis=0))
 def stars(p,transposed=False):
  lines=p.read_text().splitlines();j=next(i for i,x in enumerate(lines) if 'TOTAL POTENTIAL IN INTERSTITIAL' in x)+2;n=int(lines[j].split()[0]);a=np.zeros(shape,complex);mask=np.zeros(shape,bool)
  for line in lines[j+1:j+1+n]:
   g=np.array([int(line[i:i+5]) for i in [3,8,13]]);z=float(line[18:37])+1j*float(line[37:56]);star={}
   for rot,tau in zip(s.timat,s.tau):
    gg=tuple(((rot.T if transposed else rot)@g)%shape);star.setdefault(gg,[]).append(np.exp(2j*np.pi*np.dot(tau,g)))
   for gg,ph in star.items():
    phase=np.mean(ph);a[gg]=z*(phase if p.suffix=='.vtotal' else phase.conjugate())/(len(star) if p.suffix=='.vtotal' else 1);mask[gg]=True
  return a,mask
 def theta(h):
  t=np.zeros(shape,complex)
  for atom in atoms:
   x=qnorm*atom.rmt/(1+h);rr=np.ones_like(x);np.divide(3*spherical_jn(1,x),x,out=rr,where=x!=0)
   for p in atom.positions:t-=4*np.pi*atom.rmt**3/(3*vol*(1+h)**3)*rr*np.exp(-2j*np.pi*np.einsum('a,aijk->ijk',p,xyz))
  t[0,0,0]+=1;return t
 def product(a,b):return np.fft.fftn(np.fft.ifftn(a).real*np.fft.ifftn(b).real)*np.prod(shape)
 th=theta(0);recon={}
 for trans in [False,True]:
  vc,mask=stars(r/'r2v_center'/(case+'.vtotal'),trans);wc,_=stars(r/'r2v_center'/(case+'.vns'),trans)
  recon[str(trans)]=float(np.max(abs((product(vc,th)-wc.conj())[mask])))
 transpose=min([False,True],key=lambda t:recon[str(t)])
 assert recon[str(transpose)]<1e-8,(case,recon)
 vc,mask=stars(r/'r2v_center'/(case+'.vtotal'),transpose)
 vp,_=stars(r/'r2v_plus'/(case+'.vtotal'),transpose);vm,_=stars(r/'r2v_minus'/(case+'.vtotal'),transpose)
 wp,_=stars(r/'plus'/(case+'.vns'),transpose);wm,_=stars(r/'minus'/(case+'.vns'),transpose)
 spec=__import__('run_cases').specs[case];base=Path(spec[0])
 cp=np.array([float(x) for x in (base/spec[3]/(case+'.struct')).read_text().splitlines()[3].split()[:3]]);cm=np.array([float(x) for x in (base/spec[2]/(case+'.struct')).read_text().splitlines()[3].split()[:3]])
 h=(np.mean(cp/cell)-np.mean(cm/cell))/2
 true=product((vp-vm)/(2*h),th);true[~mask]=0
 raw=((wp-wm)/(2*h)).conj();boundary=product(vc,(theta(1e-5)-theta(-1e-5))/(2e-5));boundary[~mask]=0
 print(case,'FFT checks',recon,'step',h,flush=True)
 paths=sorted(c.glob('dmft_sternheimer.k????????.bin'));nk=len(paths);assert nk>0
 vecs={n:vectors(r/n/(case+'.vector'),s.nat) for n in ['complete_center','center','minus','plus','advection_minus','advection_plus']}
 matrices={n:[] for n in ['raw','advection','interstitial_raw','interstitial_true','interstitial_boundary','response','norm_response','energies']}
 checks={'metric':0.,'radial_headers':0.,'center_H':0.,'min_window_gap':1e99,'degenerate_pairs_skipped':0}
 for path in paths:
  d=read_dump(path);k=d.kpoint;nb=d.nstored;no=d.nbasis-d.nlo
  norm=read_normalization(c/f'dmft_sternheimer_norm.k{k:08d}.bin');rr=read_response_matrix(c/f'stress_response.k{k:08d}.bin').density;rn=np.zeros_like(rr)
  W=np.arange(d.dm_nemin-d.nemin,d.dm_nemax-d.nemin+1);out=np.setdiff1d(np.arange(nb),W)
  for j,iw in enumerate(W):
   den=d.energies[iw]-d.energies[out];checks['min_window_gap']=min(checks['min_window_gap'],float(np.min(abs(den))))
   use=abs(den)>1e-10;checks['degenerate_pairs_skipped']+=int(np.sum(~use));value=np.zeros(len(out),complex)
   value[use]=np.einsum('aij,ij->a',norm.projector[out[use]],norm.gradient[j])/den[use]
   rn[out,iw]=value;rn[iw,out]=value.conj()
  e0,v0,l0,rad0=next(vecs['complete_center']);assert len(e0)==len(v0);assert np.array_equal(l0,d.labels)
  inv=solve(v0,np.eye(len(v0)),check_finite=False);t0=inv@d.coefficients[:,:,0];dual=t0.conj().T@inv
  checks['metric']=max(checks['metric'],float(np.max(abs(t0.conj().T@t0-np.eye(nb)))))
  H={}
  for name in vecs:
   if name=='complete_center':continue
   e,v,l,rad=next(vecs[name]);assert np.array_equal(l,l0);t=dual@v;H[name]=(t*e)@t.conj().T
   checks['radial_headers']=max(checks['radial_headers'],max(float(np.max(abs(x-y))) for x,y in zip(rad,rad0)))
  checks['center_H']=max(checks['center_H'],float(np.max(abs(H['center']-np.diag(d.energies)))))
  matrices['raw'].append((H['plus']-H['minus'])/(2*h));matrices['advection'].append((H['advection_plus']-H['advection_minus'])/(3e-5))
  labels=d.labels[:,:no];diff=(labels[:,:,None]-labels[:,None,:])%np.array(shape)[:,None,None];A=d.coefficients[:no,:,0]
  for name,op in [('interstitial_raw',raw),('interstitial_true',true),('interstitial_boundary',boundary)]:matrices[name].append(A.conj().T@op[tuple(diff)]@A)
  matrices['response'].append(rr);matrices['norm_response'].append(rn);matrices['energies'].append(d.energies)
  if k%30==0:print(case,'k',k,checks,flush=True)
 matrices={n:np.array(v) for n,v in matrices.items()}
 matrices['MT_physical']=matrices['raw']-matrices['interstitial_raw'];matrices['MT_only']=matrices['MT_physical']+matrices['advection'];matrices['total']=matrices['MT_only']+matrices['interstitial_true']
 vals={}
 for name in ['raw','MT_physical','advection','interstitial_true','interstitial_boundary','MT_only','total']:
  vv={key:float(np.einsum('kij,kji->',matrices[key],matrices[name]).real) for key in ['response','norm_response']};vv['total_Ry']=vv['response']+vv['norm_response'];vv['GPa']=vv['total_Ry']*fac;vals[name]=vv
 analytic=float(np.trace(tensor(c/(case+'.scf0'),':STR_DENS')-tensor(c/(case+'.scfc'),':STR_CORE')-tensor(c/(case+'.scf2'),':STR_VAL'))*14710.507641078/3)
 result={'case':case,'nk':nk,'bands':nb,'window':[d.dm_nemin,d.dm_nemax],'volume':vol,'strain':h,'checks':checks,'FFT_reconstruction':recon,'components':vals,'analytic_GPa':analytic,'corrected_GPa':analytic+vals['total']['GPa'],'MT_only_GPa':analytic+vals['MT_only']['GPa']}
 (r/'results.json').write_text(json.dumps(result,indent=2,default=lambda x:x.item()));np.savez_compressed(r/'matrices.npz',**matrices);print(json.dumps(result,indent=2,default=lambda x:x.item()),flush=True)
if __name__=='__main__':
 for case in sys.argv[1:] or ['MnO','FeSe']:run(case)
