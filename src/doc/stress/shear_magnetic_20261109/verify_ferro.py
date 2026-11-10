from pathlib import Path
import contextlib,io,json,os,re,sys
import numpy as np
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root/'stage'))
from wien_response_density import Density
from dmft2_spin import average_charge
work=root/'native_Ferro';case='0_800_MnTe'
up=Density(work/(case+'.clmval'),work/(case+'.struct'))
dn=Density(work/(case+'.clmvaldn'),work/(case+'.struct'))
expected=(up.vector()+dn.vector())/2
excess=[]
for suffix in ['', 'dn']:
    text=(work/(case+'.outputdmf2'+suffix)).read_text()
    excess.append(float(re.findall(r'rho-rho_expected=\s*([-.0-9]+)',text)[-1]))
mu=float((work/'EF.dat').read_text())
old=os.getcwd()
try:
    os.chdir(work)
    with (root/'Ferro_average.log').open('w') as out:
        average_charge(case,'/Users/haule/dbin',out)
finally:os.chdir(old)
actual=Density(work/(case+'.clmval'),work/(case+'.struct')).vector()
error=float(np.max(np.abs(actual-expected)))
np.testing.assert_allclose(actual,expected,atol=2e-6,rtol=2e-6)
assert abs(sum(excess)/2)<1e-6
report=dict(test='Synthetic unequal-spin MnTe-derived spectra; not a converged ferrimagnetic material',down_self_energy_shift_eV=.2,mu_eV=mu,native_branch_target=62.,physical_spin_populations=[(62+x)/2 for x in excess],total_number_error=sum(excess)/2,average_charge_max_coordinate_error=error,policy=json.loads((work/'preparation.json').read_text())['physical_policy'])
(root/'Ferro_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
