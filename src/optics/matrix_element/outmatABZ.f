SUBROUTINE OUTMATABZ(cornum,modo,nnk,nnj)
  use opme
  use bindex
  use xa3
  use xrpar
  use struk
  !ad
  !ad ___________________ CALCULATION OF MATRIX ELEMENTS _________________
  !ad
  INCLUDE 'param.inc'
  IMPLICIT REAL*8 (A-H,O-Z) 
  !ole ##### Begin #####
  LOGICAL     INFO_FLAG
  LOGICAL     REL,IM,LSO,SPIN,MME_FLAG
  !ole ##### End #####
  character*3 modo
  CHARACTER*3  OUTME
  CHARACTER*10 KNAME
  COMPLEX*16  O(3,6),OX1,OX2,OY1,OY2,OZ2,OZ1     
  complex*16  ionne,complxone
  integer(4),intent(in) :: cornum
  integer(4) :: corenum
  !     
  COMMON /SYM2/   TAU(3,NSYM),IORD,IMAT(3,3,NSYM) 
  COMMON /outop/ Ncol,icol(9)          
  COMMON /KPOI/ S,T,Z,NEMIN,NEMAX,KKZ,N,NNLO,KNAME
  COMMON /LEAD/ KFIRST,KLAST,KEND,KSTEP,KOUT,KSTOP
  COMMON /COM/  EMIN,EMAX,ELECN,EULIMIT,EDLIMIT,NK,IOUT,NSPIN,NAT,NBAND,ix,NB(NKPT),MINWAV,MAXWAV
  COMMON /MIM / MIMA(2) 
  !ole ##### Begin #####
  COMMON /CLOGIC/ LSO,SPIN,REL,MME_FLAG
  !ole ##### End #####
  !      COMMON /BINDEX/ N_(numeo),NN_(numeo),NIN(NUME,NUME) 
  !ad   write(6,*)'lso,spin:  ',LSO,SPIN 
  NEMIN=MIMA(1)
  NEMAX=MIMA(2)
  ionne=(0d0,1.0d0)      
  !ole ##### Begin #####
  INFO_FLAG=.FALSE.
  complxone=(1.0d0,0.0d0)
  IF (MME_FLAG) THEN
     WRITE(24,9010) NK,NEMIN,NEMAX,EMIN,EMAX,KNAME
  ELSE
     corenum=12+cornum
     WRITE(corenum,9010) NK,NEMIN,NEMAX,EMIN,EMAX,KNAME
     if(modo.EQ.'ON ') then
        WRITE(4,9010) NK,NEMIN,NEMAX,EMIN,EMAX,KNAME
     endif
  END IF
  !
  !ole ##### End #####
  !    
  NBINDEX=0
  
  DO NB2=NEMIN,NEMAX   ! 119
     NBINDEX=NBINDEX+1
     jkl=0
     do mm=-nnj,nnj,2
        jkl=jkl+1
        O(1,jkl)=(DOPMATX(nb2,jkl))
        O(2,jkl)=(DOPMATY(nb2,jkl))
        O(3,jkl)=(DOPMATZ(nb2,jkl))
     enddo
     !ad   write out complex matrix elements before symmetrization
     jkl=0
     do mm=-nnj,nnj,2
        jkl=jkl+1
        if(modo.EQ.'ON ') write(4,9040) core_name,NB2,O(1,jkl),O(2,jkl),O(3,jkl)!,E(core_name)-E(NB1)
     enddo
     !ole ##### Begin #####
     
     IF ((MME_FLAG).AND.(.NOT.INFO_FLAG)) THEN
        jkl=0
        do mm=-nnj,nnj,2
           jkl=jkl+1
           WRITE(24,9030) core_name,NB2,(DBLE(O(i,jkl)),AIMAG(O(i,jkl)),i=1,Ncol/2)
        enddo
     END IF
     !       write(77,*) NB2, E(NB2) 
     IF (.NOT.MME_FLAG) THEN
        !ole #####  End  #####
        !ad
        call outsyma(cornum,o,nb2,nnk,nnj)
        !ole ##### Begin #####
     END IF
     !ole #####  End  #####
     !119  CONTINUE
  ENDDO
!--------------------------------------------------------------------------

  RETURN
9010 FORMAT(/,2X,' KP:',I6,' NEMIN NEMAX : ',2I5,' dE:',2f5.2,' K:',a10 /)
9011 FORMAT(2X,' KP:',I6,' NEMIN NEMAX : ',2I5,'  IORD :',I4,'  J_value: ',I4)
9030 FORMAT(3X,A4,I4,6E13.6,F13.8)
9040 FORMAT(3X,A3,2X,I4,6E13.6)
1001 format(2X,I3,3X,F8.5)
END SUBROUTINE OUTMATABZ

SUBROUTINE OUTSYMA(cornum,O,nb2,nnk,nnj)
  !ad
  !ad  ________________ SYMMETRIZATION OF MATRIX ELEMENTS ________________                  
  !ad
  !ad   May 1999: updated by Jan Kunes
  !ad
  use xrpar
  use xa3
  IMPLICIT REAL*8 (A-H,O-Z)
  INCLUDE 'param.inc'
  COMMON /outop/  Ncol,icol(9)
  COMMON /SYM2/   TAU(3,NSYM),IORD,IMAT(3,3,NSYM)
  COMMON /SYMd/   det(NSYM)
  COMPLEX*16  O(3,6),o1(3,6)
  complex*16 o2(6,6)
  complex*16 OP(3,6)
  complex*16 zeroco,onecmplx
  integer(4), intent(in) :: cornum
  integer(4) :: corenum
  DIMENSION outm(9)
  !ad
  zeroco = (0.0d0,0.0d0)
  onecmplx = (0.0d0,1.0d0)
  syfac=1.D0/REAL(iord)
  do ii=1,6
     jkl=0
     do mm=-nnj,nnj,2
        jkl=jkl+1
        o2(ii,jkl)=zeroco
     enddo
  enddo
  !ad
  do l=1,iord   ! 772
     do ii=1,3
        jkl=0
        do mm=-nnj,nnj,2
           jkl=jkl+1
           o1(ii,jkl)=zeroco
           OP(ii,jkl)=zeroco
        enddo
     end do
     do   i=1,3	
        jkl=0	
        do mm=-nnj,nnj,2
           jkl=jkl+1
           o1(i,jkl)=O(i,jkl)
        end do
     end do
     
     jkl=0
     do mm=-nnj,nnj,2
        jkl=jkl+1
        OP(1,jkl)=O1(1,jkl)+onecmplx*O1(2,jkl)
        OP(2,jkl)=O1(1,jkl)-onecmplx*O1(2,jkl)
        OP(3,jkl)=zeroco
     enddo
     
     jkl=0
     do mm=-nnj,nnj,2
        jkl=jkl+1
        o2(1,jkl)=o2(1,jkl)+op(1,jkl)*conjg(op(1,jkl))
        o2(2,jkl)=o2(2,jkl)+op(2,jkl)*conjg(op(2,jkl))
        o2(3,jkl)=o2(3,jkl)+op(3,jkl)*conjg(op(3,jkl))
     enddo
     !---------------------------------------------------------
     !772  continue
  enddo
  !____________________________________________________
  jkl=0
  do i=1,9
     outm(i)=0.d0
  enddo
  !outm=0
  do mm=-nnj,nnj,2
     jkl=jkl+1
     outm(1)=outm(1)+Dble(o2(1,jkl))*syfac
     outm(2)=outm(2)+Dble(o2(2,jkl))*syfac
     outm(3)=outm(3)+Dble(o2(3,jkl))*syfac
  enddo
  !____________________________________________________
  
  corenum=12+cornum
  WRITE(corenum,9030)  core_name,NB2,(outm(icol(i)),i=1,Ncol)
  !cad
9030 FORMAT(3X,A4,I4,9E13.6)
9040 FORMAT(3X,I4,2X,E13.6,9E13.6)
9039 FORMAT(3X,I4,2X,E13.6,9E13.6)
9031 FORMAT(3X,A4,3I4,4E13.6)
1102 format(2X,4E13.6)
23 continue
END SUBROUTINE OUTSYMA
