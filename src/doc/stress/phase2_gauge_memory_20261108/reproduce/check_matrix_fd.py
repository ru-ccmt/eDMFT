from pathlib import Path
import sys,json,os
import numpy as np
root=Path(__file__).resolve().parent
baseline=os.environ.get('EDMFT_MATRIX_BASELINE')=='1'
sys.path[:0]=([] if baseline else [str(root/'stage')])+['/Users/haule/projects/eDMFT_new/src/python']
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_native_snapshot import read_basis_matrix
case=sys.argv[1];src=Path('/private/tmp/edmft_full_k_response_20261108')/case/'charge/0'
kw={} if baseline else dict(matrix_quadrature='wien2k24')
lib='/Users/haule/dbin/liblapw_response.so' if baseline else str(root/'build/liblapw_response.so')
with load_native_valence_state(src,case,complex_vectors=case=='SrVO3',energy_policy='fixed',library=lib,**kw) as f:
 state=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
 select,_=state.hartree.native_potential_output(tuple(f.bands.spherical.slices),{a:op.layout.harmonics for a,(op,ks) in f.bands.ns.items()},f.bands.ir[0].representatives,complex_interstitial=case=='SrVO3')
 d=np.load('/private/tmp/edmft_hartree_periodic_20261108/'+case+'_direction.npz')['density']
 h,o,raw,hq,oq=f.bands.output.unpack(f.bands.action(select.matvec(state.potential_operator.matvec(d))))
 norm=lambda xs:float(np.sqrt(sum(np.linalg.norm(x)**2 for x in xs)))
 rows={};saved={}
 for step in [.001,.0005]:
  row={}
  for name,magic,left,right in [('H',b'EDMFTHV1',h,hq),('O',b'EDMFTOV1',o,oq)]:
   errors=[];refs=[]
   for k,b in enumerate(f.bases):
    file=f"native{'.ham' if name=='H' else ''}.k{k+1:08d}.bin"
    plus=read_basis_matrix(src.parent/str(step)/file,magic);minus=read_basis_matrix(src.parent/str(-step)/file,magic)
    assert np.array_equal(plus.labels,minus.labels)
    rr,qr=b.project_matrix((plus.matrix-minus.matrix)/(2*step))
    errors += [left[k]-rr,np.sqrt(2)*(right[k]-qr)];refs += [rr,np.sqrt(2)*qr]
   saved[step,name]=refs
   row[name]=dict(relative=norm(errors)/norm(refs),error_norm=norm(errors),reference_norm=norm(refs))
  rows[step]=row;print(case,step,row,flush=True)
 extrap={}
 for name,left,right in [('H',h,hq),('O',o,oq)]:
  refs=[(4*a-b)/3 for a,b in zip(saved[.0005,name],saved[.001,name])]
  analytic=[y for a,b in zip(left,right) for y in (a,np.sqrt(2)*b)]
  extrap[name]=norm([a-b for a,b in zip(analytic,refs)])/norm(refs)
 rows['Richardson']=extrap
 print(case,'Richardson',extrap,flush=True)
 (root/(case+('_matrix_baseline' if baseline else '_matrix_fd')+'.json')).write_text(json.dumps(rows,indent=2))
