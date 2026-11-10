"""Native coupled density with only retained eigenvectors used in assembly."""
from build_reduced_material import *
import time,json
from lapw_complement_response import complement_electronic_operator
from stress_response_adjoint import PairedOperator

case=sys.argv[1];nr={'FeSe':66,'SrVO3':43}[case]
label={'FeSe':'FeSe_b509_one','SrVO3':'SrVO3_b703_one'}[case]
old=Path('/private/tmp/edmft_coupled_valence_response_20261108')
start=time.perf_counter();fixture=build(case,old/label/'central',nr)
b=fixture['bands'];spatial=fixture['spatial'];direct=fixture['direct'];model=fixture['complement']
e,inputs,output=complement_electronic_operator(model,fixed_number=False)
shape=b.coordinates.band.shape

def compress(vector):
 h,o,p=b.coordinates.unpack(vector)
 return inputs.pack(h[:,:nr,:nr],o[:,:nr,:nr],p[:,:nr],h[:,nr:,:nr],o[:,nr:,:nr])

def expand_covector(vector):
 h,o,p,hq,oq=inputs.unpack(vector)
 hh=np.zeros(shape,complex);oo=np.zeros_like(hh);pp=np.zeros(b.coordinates.raw_shape,complex)
 hh[:,:nr,:nr]=h;oo[:,:nr,:nr]=o;pp[:,:nr]=p
 hh[:,nr:,:nr]=hq;oo[:,nr:,:nr]=oq
 hh[:,:nr,nr:]=np.array(hq).swapaxes(1,2).conj();oo[:,:nr,nr:]=np.array(oq).swapaxes(1,2).conj()
 return b.coordinates.pack(hh,oo,pp)

restrict=PairedOperator(b.coordinates.size,inputs.size,compress,expand_covector,'retain RR and QR input blocks')

def density_forward(vector):
 rr,qr=output.unpack(vector);rho=np.zeros(shape,complex)
 rho[:,:nr,:nr]=rr;rho[:,nr:,:nr]=qr;rho[:,:nr,nr:]=np.array(qr).swapaxes(1,2).conj()
 return spatial.action(rho)

def density_reverse(vector):
 # Only the fixed-basis density-matrix leg belongs in R. D is added separately.
 rho,_,_=spatial.pullback(vector)
 return output.pack(rho[:,:nr,:nr],rho[:,nr:,:nr])

r=PairedOperator(output.size,len(spatial.template.vector()),density_forward,density_reverse,'spatial RR plus QR density')
op=r@e@restrict@b.operator+direct.operator
v=b.pack({a:.07*cs[0].radius*np.sin(1.3*cs[0].radius) for a,cs in fixture['contexts'].items()})
actual=op.matvec(v)
reference=np.load(old/(label+'_coupled.npz'))
report=dict(case=case,retained_bands=nr,complete_coefficient_dimension=shape[1],
       setup_and_forward_seconds=time.perf_counter()-start,
       central_density_relative=float(np.linalg.norm(spatial.value()-fixture['reference'].vector())/np.linalg.norm(fixture['reference'].vector())),
       full_spectral_relative=float(np.linalg.norm(actual-reference['analytic'])/np.linalg.norm(reference['analytic'])),
       native_fd_relative=float(np.linalg.norm(actual-reference['last_fd'])/np.linalg.norm(reference['last_fd'])),
       transpose=op.check_transpose(samples=1,tolerance=3e-10),
       independent_assembly=True)
def timing(operation):
 times=[]
 for _ in range(3):
  start=time.perf_counter();operation();times.append(time.perf_counter()-start)
 return float(np.median(times))
y=np.random.default_rng(163).normal(size=op.output_size)
report['coupled_forward_seconds']=timing(lambda:op.matvec(v))
report['coupled_reverse_seconds']=timing(lambda:op.rmatvec(y))
(ROOT/(case+'_reduced_material.json')).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
for cs in fixture['contexts'].values():
 for ctx in cs.values():ctx.close()
