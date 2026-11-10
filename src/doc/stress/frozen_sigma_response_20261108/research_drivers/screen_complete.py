"""Matrix-free screening using the COMPLETE LAPW charge-map tangent.

This reference engine samples a directional derivative at the CENTRAL cell.
It is not a pair of strained SCF calculations.  Each action runs the frozen-
Sigma charge map once (LAPW0, LAPW1, dmft2, lcore); all radial and unoccupied
LAPW response is retained.  Anderson acceleration solves (1-J)q=s without
forming a susceptibility matrix.  The analytic stored-band engine can be
substituted when its radial/complement residual has been controlled.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys,json,time
from pathlib import Path
import numpy as np
from setup import ROOT,SOURCES
from full_map_check import evaluate
from strain_source import command,WIEN
from wien_response_density import Density
from test_kernel import load
from potential_probe import PotentialProbe

def run(case,iterations=15,step=1e-5,mix=.3,history=7,label='complete_screen',h=.00025,initial_path=None):
    root=ROOT/case;src=SOURCES[case];structure=src/(case+'.struct')
    rho=Density(src/(case+'.clmsum'),structure)
    source=Density(root/f'source_{h:.7f}.density',structure)
    basefile=root/'full_map_radial_base'/(case+'.clmout')
    if basefile.exists():base=Density(basefile,structure)
    else:base=evaluate(case,root/'potential_base','full_map_radial_base')
    q=source.copy();xs=[];fs=[];result={'case':case,'step':step,'source_step':h,'method':'complete frozen-Sigma charge-map tangent, Anderson','iterations':[]}
    if initial_path is not None:q=Density(initial_path,structure)
    for i in range(iterations):
        start=time.monotonic();name=f'{label}_{i:02d}';d=root/name;d.mkdir(exist_ok=True)
        # Reuse the well-defined central potential input template.
        import shutil
        for suff in ['struct','in0','inc','inm']:
            shutil.copy2(root/'potential_base'/(case+'.'+suff),d/(case+'.'+suff))
        rho.copy().add(q,step).write(d/(case+'.clmsum'))
        command(d,'lapw0',[str(WIEN/'x'),'lapw0','-f',case])
        out=evaluate(case,d,name+'_charge')
        jv=out.add(base,-1).scale(1/step)
        target=source.copy().add(jv);residual=target.copy().add(q,-1)
        x=q.vector();f=residual.vector();norm=float(np.linalg.norm(f))
        row={'iteration':i,'residual_norm':norm,'relative_residual':norm/max(1.,float(np.linalg.norm(source.vector()))),'seconds':time.monotonic()-start}
        result['iterations'].append(row);print(row,flush=True)
        q.write(root/(name+'_q.density'));target.write(root/(name+'_target.density'))
        (root/(label+'.json')).write_text(json.dumps(result,indent=2))
        if norm<.02:break
        if i==iterations-1:break
        xs.append(x.copy());fs.append(f.copy())
        xs=xs[-history-1:];fs=fs[-history-1:]
        update=x+mix*f
        if len(xs)>1:
            dx=np.column_stack([xs[j+1]-xs[j] for j in range(len(xs)-1)])
            df=np.column_stack([fs[j+1]-fs[j] for j in range(len(fs)-1)])
            coeff=np.linalg.lstsq(df,f,rcond=1e-10)[0]
            update-=(dx+mix*df)@coeff
        q=q.from_vector(update)
    # Contract only once: iterations do not need the costly frozen-basis
    # operator probe used to express delta V in the central KS gauge.
    f,dumps=load(case);p=PotentialProbe(case,dumps)
    kv,_=p.apply_kernel(q,label+'_final_potential',step=step)
    affpath=root/f'affine_operator_{h:.7f}.npy'
    if affpath.exists():aff=np.load(affpath)
    else:aff=p.affine_operator(h);np.save(affpath,aff)
    if case=='SrVO3':z=np.load('/private/tmp/projector_stationarity_probe_20261107/band_matrix_evidence.npz');analytic=5.67688415;numeric=2.933836407
    else:z=np.load('/private/tmp/projector_cross_material_20261107/FeSe/matrices.npz');analytic=-.1164916715;numeric=-1.0491512951
    R=z['response']+z['norm_response'];fac=14710.507641078/(3*p.vol)
    correction=float(np.einsum('kij,kji->',R,aff+kv).real)*fac
    ref=Density(root/'reference_q.density',structure)
    result.update(correction_GPa=correction,corrected_stress_GPa=analytic+correction,numeric_stress_GPa=numeric,error_vs_numeric_GPa=analytic+correction-numeric,reference_density_error=float(np.linalg.norm(q.copy().add(ref,-1).vector())))
    q.write(root/(label+'_final_q.density'));np.save(root/(label+'_final_potential.npy'),aff+kv)
    (root/(label+'.json')).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
    return result
if __name__=='__main__':run(sys.argv[1],int(sys.argv[2]) if len(sys.argv)>2 else 15)
