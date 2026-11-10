"""Rerun LAPW0 in independent scratch copies, enabling existing energy diagnostics."""
from pathlib import Path
import os,shutil,subprocess,json,sys
root=Path(__file__).resolve().parent
old=Path('/private/tmp/edmft_xc_core_response_20261108')
env=os.environ.copy()
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',WIEN_STRESS_DEBUG='1')
for case in sys.argv[1:] or ['SrVO3','FeSe']:
    for h in [0,.01,-.01,.005,-.005]:
        src=old/case/'xc'/str(h);dst=root/case/'native_energy'/str(h)
        dst.mkdir(parents=True,exist_ok=True)
        for name in ['lapw0.def','xc_funcs.h',case+'.struct',case+'.in0',case+'.inm',case+'.clmsum']:
            shutil.copy2(src/name,dst/name)
        # All def paths must stay in this scratch working directory.
        for line in (dst/'lapw0.def').read_text().splitlines():
            name=line.split("'")[1]
            if Path(name).is_absolute() or '..' in Path(name).parts:raise ValueError(name)
        with (dst/'run.log').open('w') as log:
            subprocess.run([str(old/'native/SRC_lapw0/lapw0'),'lapw0.def'],cwd=dst,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        print(case,h,'done',flush=True)
