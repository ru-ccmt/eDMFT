"""Audit native DEN identity defects: native XC versus LDA, probe versus screened q.

Only independent scratch inputs are written. Full-mesh q is reused; this is
not another screened solve, strain calculation or material reconvergence.
"""
from pathlib import Path
import sys,os,json,re,shutil,subprocess,time
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_potential_core_state import load_potential_core_state
from lapw_functional_cancellation import SpatialPotentialPairing
from lapw_charge_response import native_charge_covectors
from wien_response_density import Density
from lapw_xc_response import MuffinTinXCResponse,InterstitialXCResponse,PeriodicXCResponse
old=Path('/private/tmp/edmft_xc_core_response_20261108')
env=os.environ.copy();env.update(WIEN_STRESS_DEBUG='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
for case in sys.argv[1:] or ['SrVO3','FeSe']:
    pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
    pair=SpatialPotentialPairing(pc.hartree,pc.structure.mult)
    g=pair.gradient(pair.density_vector(pc.input_density))
    src=old/case/'xc/0';base=Density(src/(case+'.clmsum'),src/(case+'.struct'))
    np.testing.assert_allclose(base.vector(),pc.input_density.vector(),atol=0,rtol=0)
    mt,ir=native_charge_covectors(base,pc.structure,inversion=pc.metadata['inversion'])
    report={}
    mt_lda=[]
    for a,mtxc in enumerate(pc.xc.mt):
        center=np.column_stack([pc.input_density.rad[key] for key in pc.hartree.keys if key[0]==a])
        mt_lda.append(MuffinTinXCResponse(mtxc.r,mtxc.ells,center,mtxc.Y,mtxc.weights,functional='LDA',harmonics=mtxc.harmonics))
    ir_lda=InterstitialXCResponse(pc.hartree.reanalysis,np.array(list(pc.input_density.pw.values())),
        [pc.structure.a,pc.structure.b,pc.structure.c],functional='LDA')
    lda=PeriodicXCResponse(pc.hartree,mt_lda,ir_lda)
    for direction,steps in [('probe',[.01,.005]),('screened_full_bz',[.0001,.00005])]:
        d=np.load('/private/tmp/edmft_hartree_periodic_20261108/'+case+'_direction.npz')['density'] if direction=='probe' else np.load('/private/tmp/edmft_strain_charge_20261109/full_bz/'+case+'/analytic_density.npy')
        row=dict(particle_number_derivative=float((mt+ir)@d),sectors={})
        for sector,op in [('Hxc',pc.potential_operator),('H+LDA',lda.hxc_operator)]:
            if sector=='H+LDA' and pc.metadata['functional']=='LDA':
                row['sectors'][sector]=row['sectors']['Hxc'].copy()
                continue
            analytic=float(-g@op.matvec(d));vals={};records=[]
            for h in steps:
                for a in [h,-h]:
                    dst=root/case/'density_energy'/direction/sector/str(a)
                    if direction=='probe' and sector=='Hxc':dst=root/case/'native_energy'/str(a)
                    else:
                        dst.mkdir(parents=True,exist_ok=True)
                        for name in ['lapw0.def','xc_funcs.h',case+'.struct',case+'.in0',case+'.inm']:
                            shutil.copy2(src/name,dst/name)
                        for line in (dst/'lapw0.def').read_text().splitlines():
                            name=line.split("'")[1]
                            if Path(name).is_absolute() or '..' in Path(name).parts:raise ValueError(name)
                        base.from_vector(base.vector()+a*d).write(dst/(case+'.clmsum'))
                        if sector=='H+LDA':
                            lines=(dst/(case+'.in0')).read_text().splitlines()
                            lines[0]='TOT XC_LDA'
                            (dst/(case+'.in0')).write_text('\n'.join(lines)+'\n')
                        with (dst/'run.log').open('w') as log:
                            subprocess.run([str(old/'native/SRC_lapw0/lapw0'),'lapw0.def'],cwd=dst,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
                    text=(dst/(case+'.scf0')).read_text()
                    vals[a]=float(re.findall(r':DBG_EN_DEN_TOTAL:\s*(\S+)',text)[-1])
                    if sector=='H+LDA':assert ':POT  : POTENTIAL OPTION EX_LDA EC_LDA VX_LDA VC_LDA' in text
                fd=(vals[h]-vals[-h])/(2*h)
                records.append(dict(step=h,fd_Ry=fd,defect_Ry=fd-analytic))
            rich=(4*records[1]['fd_Ry']-records[0]['fd_Ry'])/3
            row['sectors'][sector]=dict(stationary_identity_Ry=analytic,steps=records,richardson_Ry=rich,defect_Ry=rich-analytic)
        row['XC_minus_LDA_defect_Ry']=row['sectors']['Hxc']['defect_Ry']-row['sectors']['H+LDA']['defect_Ry']
        report[direction]=row
        print(case,direction,json.dumps(row),flush=True)
    (root/case/'density_energy_audit.json').write_text(json.dumps(report,indent=2)+'\n')
