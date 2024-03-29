subroutine vnsrint(jatom,isi,radf,jj)
  USE param, ONLY: nrad, labc, nloat, clight
  USE structure, ONLY: jri, dx, Rmt
  USE vns, ONLY: lvns, mvns, nvns, vaa, vab, vad, vbb, vbd, vdd
  IMPLICIT NONE
  INTEGER, intent(in) :: jatom, isi, jj
  real*8, intent(in)  :: radf(nrad,0:labc,2,2,nloat)
  interface
     Function RINT13(REL,A,B,X,Y,NRAD,DX_,JRI_,R0_) result(S) ! Calculates overlap between psi_1=(A,B) and psi_2=(X,Y) functions
       REAL*8 :: S
       LOGICAL, intent(in) :: REL   ! relativistic or not
       REAL*8, intent(in)  :: A(NRAD),B(NRAD),X(NRAD),Y(NRAD)
       INTEGER, intent(in) :: NRAD, JRI_
       REAL*8, intent(in)  :: DX_, R0_
     End Function RINT13
     Function RINT13g(C1,C2,A,B,X,Y,NRAD,DX_,JRI_,R0_) result(S)
       ! instead of (C1,C2) = (1,1/137.0359895d0**2) we allow arbitrary C1 and C2 
       REAL*8 :: S
       REAL*8, intent(in) :: C1, C2
       REAL*8, intent(in) :: A(NRAD), B(NRAD), X(NRAD), Y(NRAD)
       INTEGER, intent(in):: NRAD, JRI_
       REAL*8, intent(in) :: DX_, R0_
     End Function RINT13g
  end interface
  ! locals
  LOGICAL :: REL
  real*8  :: vlm(nrad,5),a(nrad),b(nrad)
  INTEGER :: Jri_, index, itape, ity, LL, MM, LM1, j, m, lmmax
  real*8  :: dummy, dx_, rmt2, r0_
  ! We should probably make this better, namely, using
  ! atpar/readPot.f90 on master node, and then using mpi_broadcast
  ! But since this is used only in combination with p1/2 states, we
  ! are not optimizing it yet. Work for future...
  itape=21+isi
  REWIND itape
  READ (itape,5070)
  do ity=1,jatom
     READ (itape,5010)
     READ (itape,5020)LMMAX
     index = 0
     DO LM1=1,LMMAX
        READ (itape,5030)LL,MM
        if (iabs(ll).eq.2.and.ity.eq.jatom) then
           lvns(ity)=.true.
           index=index+1
           mvns(index,ity)=isign(mm,ll)
           READ(itape,5040)(VLM(J,index),J=1,JRI(JATOM))
        else
           READ(itape,5040)(dummy,J=1,JRI(JATOM))
        endif
        nvns(jatom)=index
        READ (itape,5060)
     enddo
     READ (itape,5050)
  enddo

5070 FORMAT (/,/)
5010 FORMAT (3X)
5020 FORMAT (15X,I3,/,/)
5030 FORMAT (15X,I3,5X,I2,/)
5040 FORMAT (3X,4E19.12)
5050 FORMAT (/,/,/)
5060 FORMAT (/)
  REL=.True.
  DX_ = DX(JATOM)
  Jri_ = JRI(JATOM)
  RMT2=RMT(JATOM)
  R0_ = RMT2*DEXP(DX_*(1.D0-Jri_))
  DO index=1,nvns(jatom)
     do m=1,jri_
        a(m)=radf(m,1,1,isi,1)*vlm(m,index)
        b(m)=radf(m,1,2,isi,1)*vlm(m,index)
     enddo     
     !CALL RINT13_old(1.d0,1.0/(clight**2),A,B,radf(1,1,1,isi,1),radf(1,1,2,isi,1),vaa(jatom,index,isi),JATOM,R0_,DX_,Jri_)
     !CALL RINT13_old(1.d0,1.0/(clight**2),A,B,radf(1,1,1,isi,2),radf(1,1,2,isi,2),vab(jatom,index,isi),JATOM,R0_,DX_,Jri_)
     !CALL RINT13_old(1.d0,1.0/(clight),A,B,radf(1,1,1,isi,jj),radf(1,1,2,isi,jj),vad(jatom,index,isi),JATOM,R0_,DX_,Jri_)
     vaa(jatom,index,isi) = RINT13(REL,A,B,radf(1,1,1,isi,1),radf(1,1,2,isi,1),NRAD,DX_,JRI_,R0_)
     vab(jatom,index,isi) = RINT13(REL,A,B,radf(1,1,1,isi,2),radf(1,1,2,isi,2),NRAD,DX_,JRI_,R0_)
     vad(jatom,index,isi) = RINT13g(1.d0,1.0/(clight),A,B,radf(1,1,1,isi,jj),radf(1,1,2,isi,jj),NRAD,DX_,JRI_,R0_)
     do m=1,jri_
        a(m)=radf(m,1,1,isi,2)*vlm(m,index)
        b(m)=radf(m,1,2,isi,2)*vlm(m,index)
     enddo     
     !CALL RINT13_old(1.d0,1.0/(clight**2),A,B,radf(1,1,1,isi,2),radf(1,1,2,isi,2),vbb(jatom,index,isi),JATOM,R0_,DX_,Jri_)
     !CALL RINT13_old(1.d0,1.0/(clight),A,B,radf(1,1,1,isi,jj),radf(1,1,2,isi,jj),vbd(jatom,index,isi),JATOM,R0_,DX_,Jri_)
     vbb(jatom,index,isi) = RINT13(REL,A,B,radf(1,1,1,isi,2),radf(1,1,2,isi,2),NRAD,DX_,JRI_,R0_)
     vbd(jatom,index,isi) = RINT13g(1.d0,1.0/(clight),A,B,radf(1,1,1,isi,jj),radf(1,1,2,isi,jj),NRAD,DX_,JRI_,R0_)
     do m=1,jri_
        a(m)=radf(m,1,1,isi,jj)*vlm(m,index)
        b(m)=radf(m,1,2,isi,jj)*vlm(m,index)
     enddo
     !CALL RINT13_old(1.d0,1.d0,A,B,radf(1,1,1,isi,jj),radf(1,1,2,isi,jj),vdd(jatom,index,isi),JATOM,R0_,DX_,Jri_)
     vdd(jatom,index,isi) = RINT13g(1.d0,1.d0,A,B,radf(1,1,1,isi,jj),radf(1,1,2,isi,jj),NRAD,DX_,JRI_,R0_)
  ENDDO
END subroutine vnsrint
