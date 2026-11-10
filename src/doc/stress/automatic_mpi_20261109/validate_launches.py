from pathlib import Path
import json, os, subprocess, sys, time
root = Path(__file__).resolve().parent
sys.path[:0] = [str(root/'stage'),str(root/'stage/tests')]
from response_mpi import resolve_launches
import dmft2_step
records = []
for case, directory, ranks in [
    ('SrVO3','/Users/haule/projects/materials/SrVO3/V_plus_force',4),
    ('FeSe','/Users/haule/projects/materials/FeSe_again/FeSe_stress/FeSe_RELA',6),
]:
    os.chdir(directory)
    # Deliberately ignore the MPI-only JSON we added earlier: this validates
    # exactly the new, config-free defaults without writing the physical case.
    req = dmft2_step.request(case)
    req['config'] = {}
    prefixfile = Path('mpi_prefix.dat2') if Path('mpi_prefix.dat2').is_file() else Path('mpi_prefix.dat')
    prefix = prefixfile.read_text().splitlines()[0].strip()
    resolved = dmft2_step.preflight(case,req,'/Users/haule/dbin',prefix)
    record = dict(case=case, prefix=prefix, python=resolved['mpi'], native=resolved['native_mpi'])
    os.chdir(root)
    command = resolved['mpi']+[sys.executable,str(root/'stage/tests/mpi_response_check.py')]
    env = os.environ | {'PYTHONPATH':str(root/'stage')+':'+str(root/'stage/tests')}
    result = subprocess.run(command,env=env,capture_output=True,text=True,timeout=60)
    record.update(returncode=result.returncode,stdout=result.stdout,stderr=result.stderr)
    records.append(record)
    (root/'real_mpi_tests.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(record),flush=True)
    assert result.returncode == 0
    assert json.loads(result.stdout)['ranks'] == ranks
print('Both physical-case preflights and automatic MPI launches passed.')
