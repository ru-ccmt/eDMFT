"""Check recorded independent evidence before presenting this checkpoint."""
from pathlib import Path
import ast
import hashlib
import json

root = Path(__file__).resolve().parent
if root.name == 'reproduce':
    root = root.parent
for label in ['repository','dbin']:
    log = (root / ('tests_'+label+'.log')).read_text()
    assert 'Ran 194 tests' in log and '\nOK\n' in log
    assert 'Verified staged helper/test and installed dependency import paths.' in log
for case in ['SrVO3','FeSe']:
    r = json.loads((root / case / 'REPORT.json').read_text())
    assert abs(r['central_scalar_full_minus_retained_Ry']) < 1e-10
    assert r['joint_vs_native_pullback_error'] < 1e-12
    adj = r['screened_induced']
    assert adj['adjoint']['converged']
    assert adj['difference_Ry'] <= adj['residual_bound_Ry'] * 1.001 + 1e-12
    if case=='SrVO3':
        assert max(x['absolute_error_Ry'] for x in r['native_steps'].values()) < 2e-10
    else:
        a,b = [r['native_steps'][s]['absolute_error_Ry'] for s in ['0.001','0.0005']]
        assert 3.8 < a/b < 4.2
        assert r['richardson']['relative_error'] < 1e-7
    f = json.loads((root / 'full_bz' / case / 'REPORT.json').read_text())
    assert f['kpoints'] == {'SrVO3':75,'FeSe':140}[case]
    assert f['retained'] == {'SrVO3':43,'FeSe':66}[case]
    assert f['pullback_error_Ry'] < 1e-12
    assert not f['central']['stress_applied']
for path in (root / 'stage').rglob('*.py'):
    ast.parse(path.read_text(),filename=str(path))
for item in json.loads((root / 'installation_manifest.json').read_text()):
    src,dest = Path(item['source']),Path(item['destination'])
    sha = lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    assert sha(src)==item['source_sha256']
    assert (sha(dest) if dest.exists() else None)==item['before_sha256']
print('Independent scalar, adjoint, full-mesh pullback, 194-test suites and untouched preimages verified.')
