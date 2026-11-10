"""Add measured preparation timing and clarify source metadata, not numerics."""
from pathlib import Path
import json,re
root=Path(__file__).resolve().parent
for mesh in ['three_k','full_bz']:
 for case in ['SrVO3','FeSe']:
  path=root/mesh/case/'REPORT.json';r=json.loads(path.read_text())
  assert r['screened']['converged']
  log=(root/f'{case}_{mesh}.log').read_text()
  r['valence_preparation_seconds']=float(re.search(r'prepared ([0-9.eE+-]+)',log)[1])
  # The older band-geometry report refers to RAW geometry; the composed
  # charge source additionally passes through E, which includes normalization.
  r['geometry_source_includes_window_normalization']=False
  r['includes_window_normalization']=True
  r['metadata_note']='Raw geometric overlaps precede normalization; the shared electronic map used in this computed source includes window, metric and normalization derivatives.'
  path.write_text(json.dumps(r,indent=2)+'\n')
