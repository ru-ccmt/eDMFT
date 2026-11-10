from pathlib import Path
import shutil,subprocess
root=Path('/private/tmp/edmft_native_operator_20261108')
src=Path('/private/tmp/edmft_radial_all_l_20261108/wien/SRC_lapw1')
d=root/'native_probe/SRC_lapw1';shutil.copytree(src,d)
p=d/'atpar.F';s=p.read_text()
s=s.replace('VRB(1,JATOM),RNOT,DX,JRI,ZZ(JATOM),jatom,efermi)','VR(1,JATOM),RNOT,DX,JRI,ZZ(JATOM),jatom,efermi)')
s=s.replace('basis_vr_solve(:)=VRB(:,JATOM)','basis_vr_solve(:)=VRB(:,JATOM)\n            if (JATOM.ne.1.or.L.ne.5) basis_vr_solve(:)=VR(:,JATOM)')
s=s.replace('if (basis_restore_high_l.and.L.gt.LNSMAX) then','if (basis_restore_high_l.and.L.eq.5.and.JATOM.eq.1) then')
s=s.replace('CALL OUTWINB(REL,VRB(1,JATOM)','CALL OUTWINB(REL,VR(1,JATOM)')
s=s.replace('if (basis_vsp_enabled) then\n               LMMIN=0','if (.false.) then ! isolated atom1,l5 probe: no low-l change\n               LMMIN=0')
p.write_text(s)
with (root/'native_probe/build.log').open('w') as log:
 r=subprocess.run(['make','real','-j4'],cwd=d,stdout=log,stderr=subprocess.STDOUT)
r.check_returncode();print(d/'lapw1')
