from pathlib import Path
import json,hashlib,difflib,tarfile
root=Path(__file__).resolve().parent;repo=Path('/Users/haule/projects/eDMFT_new');dbin=Path('/Users/haule/dbin');checkpoint=repo/'src/doc/stress/xc_core_response_20261108'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
# Only instrumented differences, not native WIEN source trees or executables.
patch=[]
for directory,name in [('SRC_lcore','resld.f'),('SRC_lapw0','xcpot1.F'),('SRC_lapw0','lapw0.F')]:
 a=Path('/Users/haule/wien2k_24')/directory/name;b=root/'native'/directory/name
 patch.extend(difflib.unified_diff(a.read_text().splitlines(True),b.read_text().splitlines(True),fromfile='a/'+directory+'/'+name,tofile='b/'+directory+'/'+name))
(root/'native_reference_ONLY.patch').write_text(''.join(patch))
# Compact central geometry and core states, plus physical file endpoints.
paths=[]
for case in ['FeSe','SrVO3']:
 for part in ['core','hxc_core']:
  for step in ['0','0.01','-0.01','0.005','-0.005']:
   folder=root/case/part/step
   paths.extend(folder/name for name in ['core_native.bin',case+'.struct',case+'.inc',case+'.vsp'])
 for step in ['0','0.001','-0.001','0.0005','-0.0005']:
  folder=root/case/'xc'/step
  paths.extend(folder/name for name in [case+'.struct',case+'.in0',case+'.clmsum',case+'.vsp',case+'.vns'])
 paths.extend((root/case/'xc/0').glob('xc_geometry_*.bin'))
 paths.append(root/case/'xc/0/xc_native.bin')
 paths.extend(root.glob(case+'*.npz'))
 paths.append(Path('/private/tmp/edmft_hartree_periodic_20261108')/(case+'_direction.npz'))
with tarfile.open(root/'compact_native_evidence.tar.gz','w:gz') as tar:
 for p in paths:
  arc=str(p.relative_to(root)) if p.is_relative_to(root) else 'input_directions/'+p.name
  tar.add(p,arcname=arc,recursive=False)
(root/'EVIDENCE_SHA256.json').write_text(json.dumps({str(p):sha(p) for p in paths},indent=2)+'\n')
rows=[]
def add(src,dst):rows.append(dict(source=str(src),target=str(dst),old_sha256=sha(dst),new_sha256=sha(src)))
source_diff=[]
for name in ['lapw_hartree_periodic.py','lapw_xc_response.py','lapw_core_response.py']:
 src=root/'stage'/name;target=repo/'src/python'/name
 add(src,target);add(src,dbin/name)
 source_diff.extend(difflib.unified_diff(target.read_text().splitlines(True) if target.exists() else [],src.read_text().splitlines(True),fromfile='a/src/python/'+name,tofile='b/src/python/'+name))
for src in sorted((root/'stage/tests').glob('test_*.py')):
 dst=repo/'src/python/tests'/src.name;add(src,dst)
 source_diff.extend(difflib.unified_diff(dst.read_text().splitlines(True) if dst.exists() else [],src.read_text().splitlines(True),fromfile='a/src/python/tests/'+src.name,tofile='b/src/python/tests/'+src.name))
for src,dst in [(root/'stage/stress_response_roadmap_20261108.txt',repo/'src/doc/stress/stress_response_roadmap_20261108.txt'),(root/'stage/hartree_periodic_REPORT.txt',repo/'src/doc/stress/hartree_periodic_response_20261108/REPORT.txt')]:
 add(src,dst)
 source_diff.extend(difflib.unified_diff(dst.read_text().splitlines(True),src.read_text().splitlines(True),fromfile='a/'+str(dst.relative_to(repo)),tofile='b/'+str(dst.relative_to(repo))))
(root/'xc_core_response.patch').write_text(''.join(source_diff))
names=['REPORT.txt','full_tests.log','xc_core_tests.log','core_material.json','xc_mt_material.json','xc_periodic_material.json','hxc_core_material.json','FeSe_composition.json','SrVO3_composition.json','FeSe_fixed_n_roots.json','SrVO3_fixed_n_roots.json','SrVO3_before_file_conjugation.json','FeSe_full_l_density_control.json','EVIDENCE_SHA256.json','compact_native_evidence.tar.gz','native_reference_ONLY.patch','xc_core_response.patch','prepare_native.py','setup_endpoints.py','setup_hxc_core.py','extra_charge_endpoints.py','setup_charge.py','fixed_n_endpoints.py','check_core.py','check_xc.py','check_xc_periodic.py','check_hxc_core.py','check_composition.py','run_all_tests.py','package_patch.py','install_after_approval.py']
for name in names:add(root/name,checkpoint/name)
(root/'installation_manifest.json').write_text(json.dumps(rows,indent=2)+'\n')
# Inputs whose conventions are relied on must not change unnoticed during review.
expected={str(repo/'src/python'/name):sha(repo/'src/python'/name) for name in ['lapw_hartree_periodic.py','lapw_potential_response.py','lapw_hartree_radial.py','lapw_density_response.py','lapw_response_chain.py','frozen_sigma_tangent.py','stress_response_adjoint.py']}
(root/'expected_before_installation.json').write_text(json.dumps(expected,indent=2)+'\n')
print('Prepared',len(rows),'guarded file copies; patch',root/'xc_core_response.patch')
print('Evidence archive bytes',(root/'compact_native_evidence.tar.gz').stat().st_size)
