"""Verify staged preimages, native scalar checks and regression evidence."""
from pathlib import Path
import hashlib,json,re
root=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for item in json.loads((root/'installation_manifest.json').read_text()):
    source,dest=Path(item['source']),Path(item['destination'])
    assert sha(source)==item['source_sha256']
    assert (sha(dest) if dest.exists() else None)==item['before_sha256']
for rel,digest in json.loads((root/'evidence/SHA256.json').read_text()).items():
    assert sha(root/'evidence'/rel)==digest,rel
for name in ['tests_repository.log','tests_dbin.log']:
    log=(root/name).read_text()
    assert 'Ran 239 tests' in log and log.rstrip().endswith('OK')
    assert 'skipped' not in log.lower()
for case in ['SrVO3','FeSe']:
    core=json.loads((root/case/'REPORT.json').read_text())
    assert abs(core['steps']['0.01']['core_error_Ry'])<1e-8
    strain=json.loads((root/case/'core_strain.json').read_text())
    assert abs(strain['richardson_error_Ry'])<3e-8
    xc=json.loads((root/case/'xc_scalar.json').read_text())
    assert abs(xc['central_error_Ry'])<5.1e-9
    assert all(abs(row['XC_minus_LDA_error_Ry'])<2e-7 for row in xc['directions'].values())
for p in (root/'stage').rglob('*.py'):compile(p.read_text(),str(p),'exec')
print('Verified preimages, evidence checksums, 239+239 tests, native core/strain/XC scalars. No files installed.')
