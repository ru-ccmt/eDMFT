"""Additional small-step endpoints, isolated from prior checkpoints and runs."""
from pathlib import Path
import sys,shutil,json
root=Path(__file__).resolve().parent
sys.path[:0]=['/private/tmp/edmft_screened_native_20261108','/Users/haule/projects/eDMFT_new/src/python']
from native_reference import NativeReference
for case in ['SrVO3','FeSe']:
 n=NativeReference(case,'three_k')
 old=n.work;n.work=root/'native_extra'/case
 (n.work/'central').mkdir(parents=True,exist_ok=True)
 shutil.copy2(old/'central/lapw1.fixgvec',n.work/'central/lapw1.fixgvec')
 h=.0000625
 for label,sign in [('m',-1),('p',1)]:
  n.evaluate(f'source_{label}_{h}',scale=1+sign*h)
 (n.work/'timings.json').write_text(json.dumps(n.timings,indent=2)+'\n')
