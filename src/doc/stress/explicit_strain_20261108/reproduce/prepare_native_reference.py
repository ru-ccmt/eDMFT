"""Build only a scratch native precision reference, never the installed LAPW0."""
from pathlib import Path
import hashlib,shutil,subprocess

root=Path(__file__).resolve().parent
if root.name=='reproduce':root=root.parent
original=Path('/private/tmp/edmft_xc_core_response_20261108/native/SRC_lapw0')
target=root/'native/SRC_lapw0'
if target.exists():
    raise RuntimeError('Native reference already exists; inspect it before rebuilding')
shutil.copytree(original,target)
subprocess.run(['patch','--batch','--forward','-p1','-i',str(root/'native_readstruct_precision.patch')],cwd=target,check=True)
with (root/'native_build.log').open('w') as log:
    subprocess.run(['make','seq'],cwd=target,stdout=log,stderr=subprocess.STDOUT,check=True)
print(target/'lapw0',hashlib.sha256((target/'lapw0').read_bytes()).hexdigest())
