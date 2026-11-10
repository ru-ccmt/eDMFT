from pathlib import Path
import json,subprocess,hashlib
r=Path('/private/tmp/edmft_native_operator_20261108');out={}
for kind in ['h','basis']:
 rows=[]
 for p in sorted((r/'refactor_checks').glob('*.input')):
  proc=subprocess.run([str(r/('installed_check_'+kind)),str(p)],capture_output=True,text=True,check=True)
  rows.append(dict(case=p.stem,values=list(map(float,proc.stdout.split()))))
 out[kind]=dict(cases=len(rows),rows=rows,maximum_errors=[max(row['values'][i] for row in rows) for i in range(4 if kind=='h' else 5)])
paths=[Path('/Users/haule/projects/eDMFT_new/src/dmft2/atpar')/name for name in ['libatpar.a','lapw_radial_response.f90','lapw_radial_hamiltonian.f90']]
out['installed_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
(r/'installed_native.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:dict(cases=v['cases'],maximum_errors=v['maximum_errors']) for k,v in out.items() if k in ['h','basis']},indent=2))
