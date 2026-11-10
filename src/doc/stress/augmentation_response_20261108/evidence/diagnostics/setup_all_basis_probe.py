from pathlib import Path
import shutil,subprocess
r=Path('/private/tmp/edmft_augmentation_response_20261108');base=Path('/private/tmp/edmft_radial_all_l_20261108/wien/SRC_lapw1')
d=r/'native/SRC_lapw1';shutil.copytree(base,d)
p=d/'atpar.F';s=p.read_text();old='VRB(1,JATOM),RNOT,DX,JRI,ZZ(JATOM),jatom,efermi)';assert old in s
s=s.replace(old,'VR(1,JATOM),RNOT,DX,JRI,ZZ(JATOM),jatom,efermi)');p.write_text(s)
for mode in ['complex','real']:
 with (r/'native'/('build_'+mode+'.log')).open('w') as log:rval=subprocess.run(['make',mode],cwd=d,stdout=log,stderr=subprocess.STDOUT)
 rval.check_returncode();print('built',mode,flush=True)
