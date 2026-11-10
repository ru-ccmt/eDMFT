"""Verify recorded analytic/native evidence and unchanged installation targets."""
from pathlib import Path
import ast,hashlib,json
root=Path(__file__).resolve().parent
if root.name=='reproduce':root=root.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for label in ['repository','dbin']:
    log=(root/('tests_'+label+'.log')).read_text()
    assert 'Ran 201 tests' in log and '\nOK\n' in log
    assert 'Verified two staged helpers, one staged test and all installed dependency paths.' in log
for case in ['SrVO3','FeSe']:
    for kind in ['pure_hydro','total_hydro']:
        report=json.loads((root/case/kind/'REPORT.json').read_text())
        a,b=[report['steps'][h]['relative'] for h in ['0.00025','0.000125']]
        assert a/b>3.8,(case,kind,a/b)
        assert report['central_relative']<1e-11
        bound={('SrVO3','pure_hydro'):2e-8,('FeSe','pure_hydro'):2e-6,
               ('SrVO3','total_hydro'):1e-7,('FeSe','total_hydro'):8e-5}[case,kind]
        assert b<bound
        assert abs(report['richardson_scalar_error_Ry'])<1e-8
for path in (root/'stage').rglob('*.py'):ast.parse(path.read_text(),filename=str(path))
for item in json.loads((root/'installation_manifest.json').read_text()):
    source,dest=Path(item['source']),Path(item['destination'])
    assert sha(source)==item['source_sha256']
    assert (sha(dest) if dest.exists() else None)==item['before_sha256']
print('201-test suites, native source/combined-potential comparisons and untouched destination preimages verified.')
