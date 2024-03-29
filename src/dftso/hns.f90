subroutine hns(ity,isi,iei,ief,tot)
  USE abcd,  ONLY: abcdlm
  USE vns,   ONLY: mvns, vaa, vab, vbb, vbd, vad, vdd, nvns
  use lolog, ONLY: mrf
  IMPLICIT NONE
  INTEGER, intent(in) :: ity, isi, iei, ief
  COMPLEX*16, intent(out) :: tot
  ! external
  REAL*8 :: gaunt1
  !
  COMPLEX*16, parameter :: imag = (0.d0,1.d0)
  COMPLEX*16 :: rsum, imag1
  INTEGER    :: index, indf, indi, irf, m, minu, mm
  !
  irf=mrf(1,ity)
  tot=(0.d0,0.d0)
  do index=1,nvns(ity)
     m=mvns(index,ity)
     mm=iabs(mvns(index,ity))
100  continue
     do indi=2,4
        indf=indi+mm
        if ((indf).ge.2.and.(indf).le.4) then
           rsum=(vaa(ity,index,isi)*conjg(abcdlm(1,indf,ief,isi))*abcdlm(1,indi,iei,isi)+  &
                vbb(ity,index,isi)*conjg(abcdlm(2,indf,ief,isi))*abcdlm(2,indi,iei,isi)+  &
                vdd(ity,index,isi)*conjg(abcdlm(irf,indf,ief,isi))*abcdlm(irf,indi,iei,isi)+  &
                vab(ity,index,isi)*conjg(abcdlm(1,indf,ief,isi))*abcdlm(2,indi,iei,isi)+  &
                vab(ity,index,isi)*conjg(abcdlm(2,indf,ief,isi))*abcdlm(1,indi,iei,isi)+  &
                vad(ity,index,isi)*conjg(abcdlm(1,indf,ief,isi))*abcdlm(irf,indi,iei,isi)+  &
                vad(ity,index,isi)*conjg(abcdlm(irf,indf,ief,isi))*abcdlm(1,indi,iei,isi)+  &
                vbd(ity,index,isi)*conjg(abcdlm(2,indf,ief,isi))*abcdlm(irf,indi,iei,isi)+  &
                vbd(ity,index,isi)*conjg(abcdlm(irf,indf,ief,isi))*abcdlm(2,indi,iei,isi))*  &
                gaunt1(1,2,1,indf-3,mm,indi-3)
                 
           IF (MM .NE. 0) THEN
              MINU = 1
              IMAG1 = (1.0D+0,0.0D+0)
              IF (M.LT.0) THEN
                 IMAG1=(-1.0D+0,0.0D+0)
                 MINU=-1
              ENDIF
              IF (MOD(MM,2).EQ.1) THEN
                 IMAG1=-IMAG1
                 MINU=-MINU
              ENDIF
              IF (MM.GT.0) MINU = 1
              rsum=rsum*imag1*dble(minu)/sqrt(2.d0)
           END IF
           tot=tot+rsum
        end if
     enddo
     if (mm.gt.0) then
        mm=-mm
        goto 100
     end if
  end do
end subroutine hns
