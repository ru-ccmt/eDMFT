from pathlib import Path
import json
import subprocess

base=Path(__file__).resolve().parent
root=Path('/Users/haule/projects/eDMFT_new')
work=base/'production_link'
work.mkdir(exist_ok=False)
atpar=root/'src/dmft2/atpar'
exe=work/'test_lapw_radial_response'
build=subprocess.run(['/opt/homebrew/bin/gfortran','-O2','-fopenmp','-fcheck=bounds,pointer,mem,do',
    '-ffpe-trap=invalid,zero,overflow','-I'+str(atpar),
    str(root/'src/dmft2/tests/test_lapw_radial_response.f90'),str(atpar/'libatpar.a'),'-o',str(exe)],
    cwd=work,capture_output=True,text=True)
(work/'build.log').write_text(build.stdout+build.stderr)
build.check_returncode()
rows=[]
for path in sorted((base/'installed_native').glob('*.input')):
    test=subprocess.run([str(exe),str(path)],cwd=work,capture_output=True,text=True)
    (work/(path.stem+'.log')).write_text(test.stdout+test.stderr)
    test.check_returncode()
    values=list(map(float,test.stdout.split()))
    rows.append(dict(case=path.stem,central_relative=values[0],basis_fd_relative=values[1],forms_fd_relative=values[2],
        basis_adjoint_relative=values[3],forms_adjoint_relative=values[4]))
assert len(rows)==32
report=dict(library=str(atpar/'libatpar.a'),tests=len(rows),rows=rows,
    maximum_errors={key:max(row[key] for row in rows) for key in list(rows[0])[1:]})
(work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2))
