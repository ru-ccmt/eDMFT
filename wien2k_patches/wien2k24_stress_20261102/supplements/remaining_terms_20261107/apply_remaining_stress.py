#!/usr/bin/env python3
"""Apply the independently audited 2026-11-07 LAPW0/lcore supplement.

Only known pre-audit sources are accepted. Compiler/library settings in the
site's Makefile are preserved; only gga_metric_tensor.o is added to OBJS.
No build, installation, or reference calculation is performed by this script.
"""
import argparse, datetime, hashlib, json, pathlib, shutil, subprocess, tempfile
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('wienroot',type=pathlib.Path)
p.add_argument('--apply',action='store_true',help='apply after verifying all inputs; default checks only')
a=p.parse_args();bundle=pathlib.Path(__file__).resolve().parent
manifest=json.loads((bundle/'source_changes.json').read_text());root=a.wienroot.resolve()
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
files={k:v for k,v in manifest.items() if k!='SRC_lapw0/Makefile'}
make=root/'SRC_lapw0/Makefile';text=make.read_text()
old='kahan_summ_gxyz.o stress_gga.o';new=old+' gga_metric_tensor.o'
all_old=all(sha(root/k)==v['before_sha256'] for k,v in files.items())
all_new=all(sha(root/k)==v['after_sha256'] for k,v in files.items())
if all_new and new in text:
 print('Already contains the complete remaining-terms audit supplement.');raise SystemExit
if not all_old or new in text or text.count(old)!=1:
 bad=[k for k,v in files.items() if sha(root/k)!=v['before_sha256']]
 p.error('Unknown or partially modified sources; refusing to overwrite. Files: '+', '.join(bad or ['SRC_lapw0/Makefile']))
# Validate/apply the patch in a disposable copy before touching WIENROOT.
with tempfile.TemporaryDirectory(prefix='wien_stress_patch_') as temp:
 tmp=pathlib.Path(temp)
 for k in manifest:
  src=root/k;dst=tmp/k;dst.parent.mkdir(parents=True,exist_ok=True)
  if src.exists():shutil.copy2(src,dst)
 with (bundle/'remaining-stress-audit.patch').open() as f:
  subprocess.run(['patch','--batch','--forward','-p1'],cwd=tmp,stdin=f,check=True,stdout=subprocess.PIPE)
 for k,v in files.items():
  if sha(tmp/k)!=v['after_sha256']:raise RuntimeError('Patched source hash mismatch: '+k)
 if (tmp/'SRC_lapw0/Makefile').read_text()!=text.replace(old,new):raise RuntimeError('Unexpected Makefile change')
 print('PASS: exact source hashes and portable Makefile patch.')
 if not a.apply:raise SystemExit
 backup=root/('stress_remaining_backup_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
 for k in manifest:
  src=root/k;dst=backup/k
  if src.exists():dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
 for k in manifest:shutil.copy2(tmp/k,root/k)
 print('Applied. Backups:',backup)
 print('Next build SRC_lapw0 and SRC_lcore using your site compiler settings, then install those executables.')
