SUBROUTINE HashTableSize(kmax, KP_equiv, taup)
  use out, only      : NKK
  use structure, only: iord, imat, tau
  use totpot, only   : KPxyz
  IMPLICIT NONE
  INTEGER, intent(out)   :: kmax(3)
  INTEGER, intent(inout) :: KP_equiv(3,iord)  ! for temporary storage
#ifdef _REAL  
  REAL*8,  intent(inout)   :: taup(iord)           ! for temporary storage
#else
  COMPLEX*16, intent(inout):: taup(iord)
#endif  
  !        Local Scalars
  INTEGER  :: NST, i, jk, js
  !        External Subroutines
  EXTERNAL           STERN
  !
  kmax(1:3)=1
  ! This will compute the size of the array we need to store warp!
  DO jk=1,NKK
     CALL STERN(NST,IORD,IMAT,KPxyz(:,jk),TAU,KP_equiv,TAUP)
     DO js=1,NST
        do i=1,3
           kmax(i) = max( kmax(i), abs(KP_equiv(i,js)) )
        enddo
     ENDDO
  enddo
END SUBROUTINE HashTableSize


SUBROUTINE SetHashTable(WARP, warp1, warp2, warp3, KP_equiv, taup)
  use out, only      : NKK
  use structure, only: iord, imat, tau
  use totpot, only   : POTK, KPxyz
  use mpi, only : Qprint
  IMPLICIT NONE
  INTEGER, intent(in)    :: warp1, warp2, warp3
#ifdef _REAL  
  REAL*8,     intent(out) :: warp(-warp1:warp1,-warp2:warp2,-warp3:warp3)
#else
  COMPLEX*16, intent(out) :: warp(-warp1:warp1,-warp2:warp2,-warp3:warp3)
#endif  
  INTEGER, intent(inout) :: KP_equiv(3,iord)    ! for temporary storage
#ifdef _REAL  
  REAL*8,  intent(inout)   :: taup(iord)           ! for temporary storage
#else
  COMPLEX*16, intent(inout):: taup(iord)
#endif  
  !        Local Scalars
  INTEGER  :: NST, jk, js
  !        External Subroutines
  EXTERNAL           STERN

  warp(-warp1:warp1,-warp2:warp2,-warp3:warp3) = 0.0D0
  
  DO jk=1,NKK
     CALL STERN(NST,IORD,IMAT,KPxyz(:,jk),TAU,KP_equiv,TAUP)
     if(jk.lt.9 .and. Qprint) then
        WRITE (6,6000) KPxyz(1:3,jk), NST
     endif
     DO js=1,NST
#ifdef _REAL           
        WARP(KP_equiv(1,js),KP_equiv(2,js),KP_equiv(3,js)) = POTK(JK)*TAUP(JS)
#else           
        WARP(KP_equiv(1,js),KP_equiv(2,js),KP_equiv(3,js)) = POTK(JK)*conjg(TAUP(JS))
#endif           
        IF (jk.LT.9 .and. Qprint) THEN
           WRITE(6,6040) js, KP_equiv(1,js),KP_equiv(2,js),KP_equiv(3,js),TAUP(js),WARP(KP_equiv(1,js),KP_equiv(2,js),KP_equiv(3,js))
        ENDIF
     ENDDO
  ENDDO
  !
6000 FORMAT(3X,'K=',3I5,'  IND=',I2)
#ifdef _REAL
6040 FORMAT(20X,I2,'. WAVE=',3I5,'    TAUP=',F10.5,/,48X,'WARPING=',F10.5)
#else
6040 FORMAT(20X,I2,'. WAVE=',3I5,'    TAUP=',2F10.5,/,48X,'WARPING=',2F10.5)
#endif
END SUBROUTINE SetHashTable



