from pathlib import Path
import sys,warnings,json
import numpy as np
from scipy.linalg import lu_factor
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root))
original=np.linalg.slogdet

def checked(a):
 with warnings.catch_warnings():
  warnings.simplefilter('error',RuntimeWarning)
  try: return original(a)
  except RuntimeWarning as exc:
   np.save(root/'logdet_warning_matrix.npy',a)
   with np.errstate(all='ignore'): signs,values=original(a)
   reference=[]
   for m in a:
    lu,_=lu_factor(m,check_finite=True)
    reference.append(np.log(abs(np.diag(lu))).sum())
   report=dict(warning=str(exc),shape=list(a.shape),finite_input=bool(np.isfinite(a).all()),finite_result=bool(np.isfinite(values).all()),max_condition=float(max(np.linalg.cond(m) for m in a)),difference_from_scipy_lu=float(np.max(abs(values-reference))),minimum_singular_value=float(min(np.linalg.svd(m,compute_uv=False).min() for m in a)))
   (root/'logdet_warning.json').write_text(json.dumps(report,indent=2)+'\n')
   print(report,flush=True)
   raise
np.linalg.slogdet=checked
from check_material import run
run('SrVO3')
