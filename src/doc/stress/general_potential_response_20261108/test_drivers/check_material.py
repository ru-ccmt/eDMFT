"""All native potential sectors with projected complement and independent FD."""
from pathlib import Path
import sys,time,json
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
# Preload staged implementations before the archived material fixture adjusts
# its own diagnostic import path. No installed code is changed during staging.
import lapw_spherical_response,lapw_nonspherical_response
from lapw_potential_response import NativeInterstitialPotential,NativePotentialBands
from lapw_complement_response import complement_electronic_operator
from stress_response_adjoint import PairedOperator
from wien_response_density import Density
from check_oneparticle_stress_fd import read_basis_matrix
sys.path.insert(0,'/private/tmp/edmft_complement_response_20261108')
from build_reduced_material import build

case=sys.argv[1];nr={'FeSe':66,'SrVO3':43}[case]
label={'FeSe':'FeSe_b509_one','SrVO3':'SrVO3_b703_one'}[case]
source=Path('/private/tmp/edmft_coupled_valence_response_20261108')/label/'central'
start=time.perf_counter();fixture=build(case,source,nr)
b=fixture['bands'];spatial=fixture['spatial'];direct=fixture['direct'];s=fixture['mat'].s
lines=(source/(case+'.vns')).read_text().splitlines()
j=next(i for i,l in enumerate(lines) if 'NUMBER OF PW' in l);nstar=int(lines[j].split()[0]);reps=[];pw=[]
for line in lines[j+1:j+1+nstar]:
    reps.append([int(line[i:i+5]) for i in [3,8,13]])
    pw.append(complex(float(line[18:37]),float(line[37:56])))
reps=np.array(reps);pw=np.array(pw)
if not fixture['mat'].complex:pw=pw.real
irs=[]
for record in fixture['records']:
    npw=record.nbasis-record.nlo
    irs.append(NativeInterstitialPotential(record.coefficients[:npw,:,0],record.labels[:,:npw].T,
                  reps,s.timat.transpose(0,2,1),s.tau,complex_values=fixture['mat'].complex))
bands=NativePotentialBands(b,irs,retained=nr)
e,inputs,outputs=complement_electronic_operator(fixture['complement'],fixed_number=False)
assert inputs.size==bands.output.size
shape=b.coordinates.band.shape

def density_forward(vector):
    rr,qr=outputs.unpack(vector);rho=np.zeros(shape,complex)
    rho[:,:nr,:nr]=rr;rho[:,nr:,:nr]=qr;rho[:,:nr,nr:]=np.array(qr).swapaxes(1,2).conj()
    return spatial.action(rho)
def density_reverse(vector):
    rho,_,_=spatial.pullback(vector)
    return outputs.pack(rho[:,:nr,:nr],rho[:,nr:,:nr])
r=PairedOperator(outputs.size,len(spatial.template.vector()),density_forward,density_reverse,'spatial retained+complement density')
op=r@e@bands.operator+direct.operator@bands.spherical_selector
sv=b.pack({a:.07*cs[0].radius*np.sin(1.3*cs[0].radius) for a,cs in fixture['contexts'].items()})
ns={a:.1*part.v for a,(part,ks) in bands.ns.items()}
vectors={sector:bands.input.pack(sv if sector=='joint' else np.zeros_like(sv),
           ns if sector in ['joint','nonspherical'] else {a:np.zeros_like(v) for a,v in ns.items()},
           .1*pw if sector in ['joint','interstitial'] else np.zeros_like(pw))
         for sector in ['nonspherical','interstitial','joint']}
report=dict(case=case,setup_seconds=time.perf_counter()-start,retained=nr,coefficient_size=shape[1],sectors={})
for sector,v in vectors.items():
    t=time.perf_counter();analytic=op.matvec(v);elapsed=time.perf_counter()-t
    h,o,p,hq,oq=bands.output.unpack(bands.action(v))
    rows=[]
    steps=([.01,.005] if sector=='nonspherical' and (ROOT/case/sector/'plus_0.01/done').exists() else [])+[.001,.0005]
    for step in steps:
        folders=[ROOT/case/sector/(sign+'_'+str(step)) for sign in ['minus','plus']]
        if not all((d/'done').exists() for d in folders):raise RuntimeError('native endpoints incomplete: '+sector)
        ds=[Density(d/(case+'.clmval'),d/(case+'.struct')).vector() for d in folders]
        fd=(ds[1]-ds[0])/(2*step);hh=[];qq=[];oo=[];oqs=[]
        for k,record in enumerate(fixture['records']):
            hm=[];om=[];c=record.coefficients[:,:,0]
            for d in folders:
                nh=read_basis_matrix(d/f'native.ham.k{k+1:08d}.bin',b'EDMFTHV1')
                no=read_basis_matrix(d/f'native.k{k+1:08d}.bin',b'EDMFTOV1')
                np.testing.assert_array_equal(nh.labels,record.labels)
                hm.append(nh.matrix);om.append(no.matrix)
            # Difference in the shared coefficient gauge before multiplying
            # by ill-conditioned LAPW coefficients. Do not subtract two large
            # independently projected central band Hamiltonians.
            dh=c.conj().T@(((hm[1]-hm[0])/(2*step))@c[:,:nr])
            do=c.conj().T@(((om[1]-om[0])/(2*step))@c[:,:nr])
            hh.append(dh[:nr]);qq.append(dh[nr:]);oo.append(do[:nr]);oqs.append(do[nr:])
        packed=bands.output.matrix.pack(h,hq);numeric=bands.output.matrix.pack(np.array(hh),qq)
        expected_o=bands.output.matrix.pack(o,oq);numeric_o=bands.output.matrix.pack(np.array(oo),oqs)
        row=dict(step=step,density_relative=float(np.linalg.norm(analytic-fd)/np.linalg.norm(fd)),
          density_maximum=float(np.max(abs(analytic-fd))),density_fd_norm=float(np.linalg.norm(fd)),
          H_relative=float(np.linalg.norm(packed-numeric)/np.linalg.norm(numeric)),
          O_difference_norm=float(np.linalg.norm(expected_o-numeric_o)),raw_derivative_norm=float(np.linalg.norm(p)))
        rows.append(row);print(case,sector,row,flush=True)
    report['sectors'][sector]=dict(rows=rows,forward_seconds=elapsed)
    np.savez(ROOT/(case+'_'+sector+'.npz'),analytic=analytic,fd=fd,input=v)
report['potential_transpose']=bands.operator.check_transpose(samples=1,tolerance=1e-10)
report['coupled_transpose']=op.check_transpose(samples=1,tolerance=1e-10)
# Number closure in the matrix representation is checked independently of
# real-space charge quadrature, which is a separate native-adapter task.
report['central_density_relative']=float(np.linalg.norm(spatial.value()-fixture['reference'].vector())/np.linalg.norm(fixture['reference'].vector()))
(ROOT/(case+'_general_potential.json')).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
for cs in fixture['contexts'].values():
    for ctx in cs.values():ctx.close()
