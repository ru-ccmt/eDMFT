from build_material import *
import json,time
from lapw_response_chain import generalized_electronic_operator
from check_oneparticle_stress_fd import read_basis_matrix
case=sys.argv[1];label=sys.argv[2] if len(sys.argv)>2 else case
keep=int(sys.argv[3]) if len(sys.argv)>3 else None
start=time.perf_counter();fixture=build(case,ROOT/label/'central',keep)
model=fixture['model'];mat=fixture['mat'];s=mat.s;bands=fixture['bands'];spatial=fixture['spatial'];direct=fixture['direct'];ctx=fixture['contexts']
dv={a:.07*cs[0].radius*np.sin(1.3*cs[0].radius) for a,cs in ctx.items()};v=bands.pack(dv)
h,o,p=bands.coordinates.unpack(bands.action(v));rho_dot,info=model.generalized_response(h,o,p,fixed_number=False)
analytic=spatial.action(rho_dot)+direct.action(v)
e,coords=generalized_electronic_operator(model,fixed_number=False);r,_=spatial.fixed_basis_operator()
op=r@e@bands.operator+direct.operator
base=spatial.value();native=fixture['reference'].vector()
report=dict(case=case,setup_seconds=time.perf_counter()-start,
   central_density_relative=float(np.linalg.norm(base-native)/np.linalg.norm(native)),
   normalization_offset_norm=float(np.linalg.norm(model.offset)),
   coupled_transpose=op.check_transpose(samples=1,tolerance=1e-10),
   bands_transpose=bands.operator.check_transpose(samples=1,tolerance=1e-10),steps=[])
print(case,'central/transpose',report,flush=True)
for step in [.002,.001,.0005,.00025]:
 if not (ROOT/label/('m_'+f'{step:.0e}')/'done').exists():continue
 ds=[];hs=[];os=[]
 for sign in ['m','p']:
  folder=ROOT/label/(sign+'_'+f'{step:.0e}')
  ds.append(Density(folder/(case+'.clmval'),folder/(case+'.struct')).vector())
  hm=[];om=[]
  for k,record in enumerate(fixture['records']):
   native_h=read_basis_matrix(folder/f'native.ham.k{k+1:08d}.bin',b'EDMFTHV1')
   native_o=read_basis_matrix(folder/f'native.k{k+1:08d}.bin',b'EDMFTOV1')
   np.testing.assert_array_equal(native_h.labels,record.labels)
   c=record.coefficients[:,:,0];hm.append(c.conj().T@native_h.matrix@c);om.append(c.conj().T@native_o.matrix@c)
  hs.append(np.array(hm));os.append(np.array(om))
 fd=(ds[1]-ds[0])/(2*step);dh=(hs[1]-hs[0])/(2*step);do=(os[1]-os[0])/(2*step)
 row=dict(step=step,density_relative=float(np.linalg.norm(fd-analytic)/np.linalg.norm(fd)),density_maximum=float(np.max(abs(fd-analytic))),
          H_relative=float(np.linalg.norm(dh-h)/np.linalg.norm(dh)),O_relative=float(np.linalg.norm(do-o)/np.linalg.norm(do)),
          FD_norm=float(np.linalg.norm(fd)),analytic_norm=float(np.linalg.norm(analytic)))
 report['steps'].append(row);print(case,row,flush=True)
output=label+('_keep'+str(keep) if keep else '')
(ROOT/(output+'_coupled.json')).write_text(json.dumps(report,indent=2)+'\n')
np.savez(ROOT/(output+'_coupled.npz'),h=h,o=o,p=p,analytic=analytic,last_fd=fd,central_density=base,native_density=native)
for cs in ctx.values():
 for c in cs.values():c.close()
