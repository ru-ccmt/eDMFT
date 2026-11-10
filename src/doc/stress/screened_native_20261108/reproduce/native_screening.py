"""Independent cold-start native Anderson solve; no analytic action or initial q."""
import sys,json,time
from pathlib import Path
import numpy as np
from native_reference import NativeReference
from frozen_sigma_screening import solve_response
case=sys.argv[1];mesh='three_k';r=NativeReference(case,mesh);work=r.work
b=np.load(work/'source_0.000125.npy');counter=0;step=5e-6
def action(q):
    global counter
    p,_=r.evaluate('cold_p_'+str(counter),r.rho.from_vector(r.rho.vector()+step*q))
    m,_=r.evaluate('cold_m_'+str(counter),r.rho.from_vector(r.rho.vector()-step*q))
    counter+=1
    return (p-m)/(2*step)
def checkpoint(q,res,row):
    print(case,json.dumps(row),flush=True)
start=time.monotonic()
q,report=solve_response(b,action,absolute_tolerance=5e-4,max_iterations=30,callback=checkpoint)
report['seconds']=time.monotonic()-start
ref=np.load(work/'analytic_density.npy')
report['density_relative_difference']=float(np.linalg.norm(q-ref)/np.linalg.norm(ref))
report['density_difference_norm']=float(np.linalg.norm(q-ref))
report['source_norm']=float(np.linalg.norm(b))
report['scope']='independent cold-start native finite-difference charge map; same local central tangent and source'
np.save(work/'native_screened_density.npy',q)
(work/'native_screening_report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
if not report['converged']:raise RuntimeError('native screened response did not converge')
