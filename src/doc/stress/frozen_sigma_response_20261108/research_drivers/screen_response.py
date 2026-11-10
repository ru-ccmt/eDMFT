"""First screened pass, then optional linear screening iterations.

The nonlinear strained sources are sampled once.  All subsequent electronic
responses use analytic G dH G in the central band representation.  LAPW0 and
a frozen-basis LAPW1 probe evaluate K acting on the density direction; lcore
supplies the small core response to this induced potential.  No impurity or
strained charge-SCF calculation is performed.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys,json,time,subprocess
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from setup import ROOT,ENV,EXE,prepare
from test_kernel import load
from potential_probe import PotentialProbe
from wien_response_density import Density
from strain_source import command

def assemble(case,model,dumps,rho,name):
    prefix=ROOT/case/(name+'_rho');model.write_density(prefix,rho,dumps)
    d=prepare(case,name)
    command(d,'response',['/opt/homebrew/bin/mpirun','-n','4',str(EXE),'dmft2.def'],{'EDMFT_RESPONSE_DENSITY_PREFIX':str(prefix)})
    return Density(d/(case+'.clmval'),d/(case+'.struct'))

def charge(d,p):
    mt=0.
    for (at,l,m),v in d.rad.items():
        if l==0 and m==0:
            x=d.grids[at];mt+=d.mult[at]*CubicSpline(x,v*np.exp(x)).integrate(x[0],x[-1])
    ir=sum(z*p.theta[tuple(np.array(g)%p.shape)].conjugate()*p.vol for g,z in d.pw.items())
    return float(mt+ir.real)

def run(case,h=0.00025,iterations=1,mix=0.3,step=1e-5):
    f,dumps=load(case);p=PotentialProbe(case,dumps);root=ROOT/case
    source=Density(root/f'source_{h:.7f}.density',p.structure)
    cache=root/f'affine_operator_{h:.7f}.npy'
    if cache.exists():aff=np.load(cache)
    else:aff=p.affine_operator(h);np.save(cache,aff)
    if case=='SrVO3':z=np.load('/private/tmp/projector_stationarity_probe_20261107/band_matrix_evidence.npz');baseline=5.67688415;numeric=2.933836407
    else:z=np.load('/private/tmp/projector_cross_material_20261107/FeSe/matrices.npz');baseline=-0.1164916715;numeric=-1.0491512951
    R=z['response']+z['norm_response'];fac=14710.507641078/(3*p.vol)
    contract=lambda op:float(np.einsum('kij,kji->',R,op).real)*fac
    history={'case':case,'source_step':h,'kernel_step':step,'central_H_replay_error_Ry':p.h0error,'source_charge_derivative':charge(source,p),'affine_correction_GPa':contract(aff),'iterations':[]}
    print('START',history,flush=True)
    q=source.copy()
    for i in range(iterations):
        start=time.monotonic();kv,core=p.apply_kernel(q,f'screen_{i:02d}_potential',step=step,core=True)
        dr,info=f.response(kv,fixed_number=True)
        val=assemble(case,f,dumps,dr,f'screen_{i:02d}_density')
        response=p.rho.zeros().add(val).add(core)
        # Linearized default MIXER charge normalization.  At the converged
        # central density alpha~1.  Keep MT response intact and remove the
        # tiny integrated leakage through the interstitial density shape.
        leakage=charge(response,p)
        base_ir=sum(z*p.theta[tuple(np.array(g)%p.shape)].conjugate()*p.vol for g,z in p.rho.pw.items()).real
        for g in response.pw:response.pw[g]-=leakage/base_ir*p.rho.pw[g]
        target=source.copy().add(response)
        residual=target.copy().add(q,-1)
        correction=contract(aff+kv)
        row={'iteration':i,'correction_GPa':correction,'corrected_stress_GPa':baseline+correction,'error_vs_numeric_GPa':baseline+correction-numeric,'induced_correction_GPa':contract(kv),'q_charge':charge(q,p),'response_valence_charge':charge(val,p),'response_core_charge':charge(core,p),'target_charge':charge(target,p),'density_residual_norm':float(np.linalg.norm(residual.vector())),'seconds':time.monotonic()-start,'resolvent':info}
        history['iterations'].append(row);print('ITER',json.dumps(row),flush=True)
        q.write(root/f'screen_{i:02d}_q.density');target.write(root/f'screen_{i:02d}_target.density');np.save(root/f'screen_{i:02d}_potential.npy',kv)
        (root/'screen_results.json').write_text(json.dumps(history,indent=2))
        q.add(residual,mix)
    return history

if __name__=='__main__':run(sys.argv[1],iterations=int(sys.argv[2]) if len(sys.argv)>2 else 1)
