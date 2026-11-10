"""Isolated read-only oracle: stop LAPW0 after Hartree, before XC.

Only a scratch copy is instrumented. The installed source/executable are never
changed. Stream exports avoid text rounding and isolate Coulomb from XC.
"""

from pathlib import Path
import shutil

root = Path(__file__).resolve().parent
src = Path("/Users/haule/wien2k_24/SRC_lapw0")
dst = root / "native/SRC_lapw0"
shutil.copytree(src, dst, dirs_exist_ok=True)
p = dst / "lapw0.F"
text = p.read_text()
marker = "! Serial output of coulomb potential"
assert text.count(marker) == 1
code = """! Scratch reference export: electronic derivative at fixed ionic geometry.
      open(997,file='hartree_native.bin',access='stream',form='unformatted',status='replace')
      write(997) nat,nkk,ncon,iff1,iff2,iff3
      write(997) vol
      write(997) kzz(:,1:nkk),inst(1:nkk),absk(1:nkk),rhok(1:nkk),potk(1:nkk)
      do jatom=1,nat
        write(997) jri(jatom),lmmax(jatom)
        write(997) lm(:,1:lmmax(jatom),jatom)
        write(997) clm(1:jri(jatom),1:lmmax(jatom),jatom)
        write(997) v(1:jri(jatom),1:lmmax(jatom),jatom)
      enddo
      call rean0(nkk,kzz,iff1,iff2,iff3,ust)
      write(997) ust(:,:,:,1)
      if (.not.allocated(rhoksp)) allocate(rhoksp(nkk,2))
      rhoksp(1:nkk,1)=potk
      rhoksp(1:nkk,2)=potk
      call setfft(nkk,iff1,iff2,iff3,rhoksp,fft)
      allocate(vkexc(nkk,2))
      call rean4(nkk,iff1,iff2,iff3,ust,fft,vkexc)
      write(997) vkexc(1:nkk,1)/vol
      close(997)
      stop 'HARTREE REFERENCE COMPLETE'
!
"""
p.write_text(text.replace(marker, code + marker))
(root / "native_instrumentation.txt").write_text(code)
