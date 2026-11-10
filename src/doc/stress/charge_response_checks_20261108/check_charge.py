"""All native potential sectors with projected complement and independent FD."""
from pathlib import Path
import sys,time,json
import numpy as np
ROOT=Path('/private/tmp/edmft_general_potential_response_20261108')
DEST=Path(__file__).resolve().parent
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
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

from wien_response_operators import charge_covectors
mt,ir=charge_covectors(fixture['reference'],source/(case+'.struct'))
native_mt=mt.copy();pos=0
for (a,l,m),values in fixture['reference'].rad.items():
    n=len(values)
    if (l,m)==(0,0):native_mt[pos:pos+n]=s.mult[a]*fixture['contexts'][a][0].overlap_metric[:,0]
    pos+=n
charge=mt+ir;native_charge=native_mt+ir
rho=fixture['complement'].rho
from scipy.special import spherical_jn
# Integrate the untruncated, alias-free PW product independently of the
# finite native output star list. This distinguishes storage truncation from
# a missed electronic or direct-basis response term.
def exact_ir(matrices):
    total=0.
    for ir_map,rr in zip(spatial.interstitial,matrices):
        shape=ir_map.shape
        g=np.stack(np.meshgrid(*[np.fft.fftfreq(n)*n for n in shape],indexing='ij'),axis=-1)
        mag=np.linalg.norm(2*np.pi*g/np.array([s.a,s.b,s.c]),axis=-1)
        theta=np.all(g==0,axis=-1).astype(complex)
        site=0
        for a,mult in enumerate(s.mult):
            z=mag*s.rmt[a];bess=np.ones_like(z);mask=z!=0
            bess[mask]=3*spherical_jn(1,z[mask])/z[mask]
            for _ in range(mult):
                theta-=4*np.pi*s.rmt[a]**3/(3*ir_map.volume)*bess*np.exp(-2j*np.pi*(g@s.pos[site]));site+=1
        grid=np.einsum('...n,...n->...',ir_map._waves(ir_map.c@rr),ir_map.psi.conj()).real
        f=np.fft.fftn(grid)/ir_map.ng
        total+=np.vdot(theta,f).real
    return float(total)

central=float(np.trace(rho,axis1=1,axis2=2).sum().real)
report={'case':case,'central_band_number':central,
 'central_spline_charge':float(charge@spatial.value()),
 'central_native_radial_charge':float(native_charge@spatial.value()),'central_full_fourier_charge':float(native_mt@spatial.value()+exact_ir(spatial.rho)),'directions':{}}
for sector,v in vectors.items():
    h,o,p,hq,oq=bands.output.unpack(bands.action(v));rows=[]
    for fixed in [False,True]:
        ee,ii,oo=complement_electronic_operator(fixture['complement'],fixed_number=fixed)
        dc=ee.matvec(bands.action(v));rr,qr=oo.unpack(dc)
        dn=r.matvec(dc)+direct.operator.matvec(bands.input.unpack(v)[0])
        count=np.trace(rr,axis1=1,axis2=2).sum().real+np.einsum('kij,kji->',rho,o).real
        fullrho=np.zeros(shape,complex);fullrho[:,:nr,:nr]=rr;fullrho[:,nr:,:nr]=qr;fullrho[:,:nr,nr:]=np.array(qr).swapaxes(1,2).conj()
        rows.append({'fixed_number':fixed,'full_fourier_number_derivative':float(native_mt@dn+exact_ir(fullrho)),'band_number_derivative':float(count),
                     'spline_number_derivative':float(charge@dn),
                     'native_radial_number_derivative':float(native_charge@dn),
                     'density_norm':float(np.linalg.norm(dn))})
        np.savez(DEST/(case+'_'+sector+('_fixedN' if fixed else '_fixedMu')+'.npz'),density=dn,charge_covector=charge)
    report['directions'][sector]=rows
print(json.dumps(report,indent=2),flush=True)
(DEST/(case+'_charge_check.json')).write_text(json.dumps(report,indent=2)+'\n')
for cs in fixture['contexts'].values():
    for ctx in cs.values():ctx.close()
