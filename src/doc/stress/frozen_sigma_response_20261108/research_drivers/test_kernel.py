import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys,time,json
from pathlib import Path
import numpy as np
sys.path[:0]=['/Users/haule/projects/eDMFT_new/src/python',str(Path(__file__).parent)]
from frozen_sigma_response import *
from check_projector_sternheimer_fd import read_dump,read_normalization
from check_oneparticle_stress_fd import read_one_particle
ROOT=Path(__file__).parent
srcs={'SrVO3':Path('/private/tmp/projector_stationarity_probe_20261107/SrVO3'),'FeSe':Path('/private/tmp/projector_cross_material_20261107/FeSe/center')}
def load(case):
 src=srcs[case];m=Metadata.read(ROOT/case/'central/response.meta');files=sorted(src.glob('dmft_sternheimer.k????????.bin'));d=[read_dump(p) for p in files];raw=[read_normalization(src/f'dmft_sternheimer_norm.k{i+1:08d}.bin').projector for i in range(len(d))]
 w=np.array([float(s.split()[5]) for s in (src/(case+'.klist')).read_text().splitlines() if s.strip() and not s.startswith('END')]);w/=w.sum()
 f=FrozenSigmaResponse(m,np.array([x.energies for x in d]),np.array(raw),np.array([x.projector for x in d]),w,np.arange(d[0].dm_nemin-d[0].nemin,d[0].dm_nemax-d[0].nemin+1))
 return f,d
if __name__=='__main__':
 case=sys.argv[1] if len(sys.argv)>1 else 'SrVO3';f,d=load(case);out={}
 print('loaded',case,'shape',f.eps.shape,'Xerror',np.max(abs(f.X-f.meta.normalization)),'offset',np.linalg.norm(f.offset),flush=True)
 cache=ROOT/case/'rho_python.npy'
 if cache.exists():rho=np.load(cache)
 else:rho=f.density();np.save(cache,rho)
 ref=[]
 for k in range(len(d)):
  rr=read_one_particle(srcs[case]/f'dmft_rho_rhoe.k{k+1:08d}.bin');r=np.zeros_like(rho[k]);r[:rr.nfull,:rr.nfull]=rr.rho;ref.append(r)
 ref=np.array(ref);out['central_density_maxerror']=float(np.max(abs(ref-rho)));out['central_N']=float(np.trace(rho,axis1=1,axis2=2).real.sum());print('central',out,flush=True)
 f.write_density(ROOT/case/'rho_central',rho,d)
 if case=='SrVO3':op=np.load('/private/tmp/projector_stationarity_probe_20261107/band_matrix_evidence.npz')['potential']
 else:op=np.load('/private/tmp/projector_cross_material_20261107/FeSe/matrices.npz')['total']
 dr,info=f.response(op,fixed_number=False);out['response']=info;np.save(ROOT/case/'test_response.npy',dr);f.write_density(ROOT/case/'rho_test',dr,d);print('analytic',info,flush=True)
 for h in [1e-4,5e-5]:
  p=f.density(op,h,f.embedding_endpoint(op,h));m=f.density(op,-h,f.embedding_endpoint(op,-h));fd=(p-m)/(2*h);err=np.max(abs(fd-dr));out[str(h)]={'maxerror':float(err),'relative':float(np.linalg.norm(fd-dr)/np.linalg.norm(dr))};print('FD',h,out[str(h)],flush=True)
 (ROOT/case/'kernel_checks.json').write_text(json.dumps(out,indent=2))
