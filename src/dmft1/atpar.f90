module matpar
  USE param
  REAL*8  :: alo(0:lomax,nloat,nrf,2)
  INTEGER :: nlo, nlov, nlon, ilo(0:lomax)
  LOGICAL :: lapw(0:lmax2)
  !REAL*8  :: RI_MAT(0:lmax2,0:lmax2,nrf,nrf,2)
  REAL*8  :: RI_MAT(nrf,nrf,0:lmax2,2)
  REAL*8  :: P(0:LMAX2,2,nrf), DP(0:LMAX2,2,nrf)
  REAL*8  :: RF1(NRAD,0:LMAX2,2,nrf),RF2(NRAD,0:LMAX2,2,nrf)
contains

  SUBROUTINE ATPAR (jatom,jtape,el_store,elo_store,is,ISPIN)
    ! There is a better atpar in dmft2 part, but somehow I did not manage to make use of it.
    ! Notice that all routines called from here are obtained from atpar library
    USE com_mpi,   ONLY: myrank, master
    USE structure, ONLY: nat, rel, jri, mult, aname, r0, dx, ZZ, Rmt
    IMPLICIT NONE
    INTEGER, intent(in) :: jatom, jtape, is, ispin
    REAL*8,  intent(in) :: el_store(0:lmax2,nat),elo_store(0:lomax,1:nloat,nat)
    !
    interface
       Function RINT13(REL,A,B,X,Y,NRAD,DX_,JRI_,R0_) result(S) ! Calculates overlap between psi_1=(A,B) and psi_2=(X,Y) functions
         REAL*8 :: S
         LOGICAL, intent(in) :: REL   ! relativistic or not
         REAL*8, intent(in)  :: A(NRAD),B(NRAD),X(NRAD),Y(NRAD)
         INTEGER, intent(in) :: NRAD, JRI_
         REAL*8, intent(in)  :: DX_, R0_
       End Function RINT13
    end interface
    LOGICAL    :: rlo(1:nloat,0:lomax)
    REAL*8     :: A(NRAD), B(NRAD), AE(NRAD), BE(NRAD)
    !REAL*8     :: emist(0:lomax,nloat),
    REAL*8     :: e(0:LMAX2),elo(0:lomax,nloat),pei(0:lmax2),Vr(nrad)
    LOGICAL    :: loor(nloat,0:lomax), qprint
    !
    INTEGER    :: idummy, i, j, imax, irf, jlo, kappa, l, m, node, nodel, nodes
    REAL*8     :: cross, dele, delei, E1, EI, FL, ovlp, r_m, trx, try
    REAL*8     :: UV, DUV, UVB, DUVB, UVE, DUVE, pi12lo, pe12lo, plo, dplo

    call get_ilo(el_store,elo_store,jatom,mult,nat,lmax2,lomax,nloat, e,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo)
    
    !---------------------------------------------------------------------  
    !.....READ TOTAL SPHERICAL POTENTIAL V(0,0) OF TAPEjtape=VSP               
    !     NORM OF V(0,0)=V00(R)*R/SQRT(4.D0*PI)     
    ! itape points to the vector file, jtape to the potential file                        
    READ(jtape,1980)
    READ(jtape,2000)IDUMMY
    READ(jtape,2031)
    READ(jtape,2022)(VR(J),J=1,jri(jatom))
    READ(jtape,2031)
    READ(jtape,2030)
    DO J=1,jri(JATOM)
       VR(J)=VR(J)/2.0D0  
    ENDDO

    P(:,is,:)=0
    DP(:,is,:)=0
    rf1(:,:,is,:)=0.d0
    rf2(:,:,is,:)=0.d0
    
    if (myrank.EQ.master) then
       if (ispin.eq.2) then
          if(is.eq.1)then
             write(6,*) 'ATPAR for spin up'
          else
             write(6,*) 'ATPAR for spin down'  
          endif
       else
          write(6,*) 'ATPAR calculation'
       endif
       
       WRITE(6,7) ANAME(JATOM)
       WRITE(6,5) E
       WRITE(6,14)
    endif
    !                                                                       
    RF1(:,:,is,:)=0
    RF2(:,:,is,:)=0
    DO l=0,LMAX2
       DELE=2.0D-3
       DELEI=0.25D0/DELE
       FL=L
       EI=E(l)/2.0d0
       !     CALCULATE ENERGY-DERIVATIVE BY FINITE DIFFERENCE                  
       !     DELE IS THE UPWARD AND DOWNWARD ENERGY SHIFT IN HARTREES          
       !                                                                       
       E1=EI-DELE
       ! OUTPUT IS: UVB,DUVB,NODEL
       CALL OUTWIN(A,B,NODEL,UVB,DUVB,REL,VR,r0(JATOM),DX(JATOM),jri(JATOM),E1,FL,zz(jatom),nrad)
       OVLP = rint13(REL,A,B,A,B,nrad,dx(jatom),jri(jatom),r0(jatom))
       TRX=1.0D0/SQRT(OVLP)
       IMAX=jri(JATOM)
       DO M=1,IMAX
          AE(M)=TRX*A(M)
          BE(M)=TRX*B(M)
       ENDDO !45 CONTINUE                                                          
       UVB=TRX*UVB                                                       
       DUVB=TRX*DUVB                                                     
       E1=EI+DELE                                                        
       ! OUTPUT IS: UVE,DUVE,NODE
       CALL OUTWIN(A,B,NODE,UVE,DUVE,REL,VR,r0(JATOM),DX(JATOM),jri(JATOM),E1,FL,zz(jatom),nrad)
       OVLP = rint13(REL,A,B,A,B,nrad,dx(jatom),jri(jatom),r0(jatom))
       TRX=1.0d0/SQRT(OVLP)
       UVE=DELEI*(TRX*UVE-UVB)
       DUVE=DELEI*(TRX*DUVE-DUVB)
       IMAX=jri(JATOM)
       DO M=1,IMAX ! 50
          AE(M)=DELEI*(TRX*A(M)-AE(M))
          BE(M)=DELEI*(TRX*B(M)-BE(M))
       ENDDO !50 CONTINUE
       ! AE and BE set! 
       !
       !     CALCULATE FUNCTION AT EI
       !
       ! OUTPUT IS: UV,DUV,NODES
       CALL OUTWIN(A,B,NODES,UV,DUV,REL,VR,r0(JATOM),DX(JATOM),jri(JATOM),EI,FL,zz(jatom),nrad)
       OVLP = rint13(REL,A,B,A,B,nrad,dx(jatom),jri(jatom),r0(jatom))
       TRX=1.0d0/SQRT(OVLP)
       P(l,is,1)=TRX*UV
       DP(l,is,1)=TRX*DUV
       IMAX=jri(JATOM)
       DO M=1,IMAX ! 60
          A(M)=TRX*A(M)
          B(M)=TRX*B(M)
       ENDDO ! 60
       ! A and B set
       !                                                                       
       !     INSURE ORTHOGONALIZATION                                          
       !                                                                       
       CROSS = rint13(REL,A,B,AE,BE,nrad,dx(jatom),jri(jatom),r0(jatom))
       TRY=-CROSS
       IMAX=jri(JATOM)
       DO M=1,IMAX ! 55
          AE(M)=(AE(M)+TRY*A(M))
          BE(M)=(BE(M)+TRY*B(M))
       ENDDO ! 55
       IMAX=jri(JATOM)
       DO I=1,IMAX  ! 80
          RF1(I,l,is,1)=A(I)                                                   
          RF2(I,l,is,1)=B(I)                                                   
          RF1(I,l,is,2)=AE(I)                                                  
          RF2(I,l,is,2)=BE(I)                                                  
       ENDDO   !80    CONTINUE                                                          
       P(l,is,2)=UVE+TRY*P(l,is,1)                                                
       DP(l,is,2)=DUVE+TRY*DP(l,is,1) 
       PEI(l) = rint13(REL,AE,BE,AE,BE,nrad,dx(jatom),jri(jatom),r0(jatom)) ! <dotu|dotu>
       if (myrank.EQ.master) WRITE(6,8) L,P(l,is,1),DP(l,is,1),P(l,is,2),DP(l,is,2)
    ENDDO ! 70
    ! Results saved for P(l,is,irf),DP(l,is,irf)
    !                         
    ! local orbitals part
    !
    DO l=0,lomax ! 170
       irf=2
       ! If lapw=True, then we compute   u_{new} = a^lo u + b^lo dotu + c^lo u^{LO}
       !   and   alo(l,is,jlo,1) contains a^lo
       !         alo(l,is,jlo,2) contains b^lo
       !         alo(l,is,jlo,jlo+2) contains c^lo
       !
       !   In this case we always have irf=jlo+2, which is either 3 or 4.
       !
       ! If lapw=False, then we compute  u_{jlo=1} = a^lo u + b^lo dotu
       !  and    alo(l,is,1,1)  contains a^lo
       !         alo(l,is,1,2)  contains b^lo
       !         alo(l,is,1,3:) is set to zero
       !
       !  then for jlo=2 or 3 we have real local orbital in APW+lo   u_{jlo=2 or 3} = a^lo u + c^lo u^{LO}
       !  and   alo(l,is,jlo,1)   contains a^lo
       !        alo(l,is,jlo,2)   is set to zero
       !        alo(l,is,jlo,jlo+1) contains c^lo
       !
       !  In this case we have irf=jlo+1, which is 3 or 4.
       !
       do jlo=1,ilo(l) ! 170
          !if (.not.loor(jlo,l)) CYCLE ! This is not real LO, but just apw+lo    ! Bug Nov 30, 2016 (Walber). This simplistic line does not work for actinides. Need to use the complicated line below.
          if (  lapw(l) .or. jlo.gt.1  ) then 
             ! We go here inside for real local orbital, i.e., everytime except when computing APW+lo.
             !  namely,  APW+lo is computed when lapw=False and jlo=1, then the first component (jlo=1) contains APW+lo functions, which are constructed from  a^{lo} u + b^{lo} dotu , but not from any real local orbital.
             ! 
             irf=irf+1
             DELE=2.0D-3                                                       
             DELEI=0.25D0/DELE                                                 
             FL=L                                                              
             EI=elo(l,jlo)/2.d0                                                         
             !                                                                       
             !     CALCULATE FUNCTION AT EI 
             IF(rlo(jlo,l)) THEN
                ei=elo(l,nloat)/2.d0
                kappa=l
                ! output: uv,duv,nodes
                CALL diracout(A,B,rel,vr,r0(jatom),dx(jatom),jri(jatom),ei,fl,kappa,uv,duv,nodes,zz(jatom),nrad)
                ! output: b
                CALL dergl(a,b,r0(jatom),dx(jatom),jri(jatom),nrad)  ! BUG corrected Nov 30, 2016 (Walber) . Added last variable nrad.
                DO m = 1, jri(jatom)
                   r_m = r0(jatom)*exp(dx(jatom)*(m-1))
                   b(m) = b(m)*r_m/(2.d0*clight+(elo(l,jlo)-2.d0*vr(m)/r_m)/(2.d0*clight))
                   b(m)=b(m)*clight
                ENDDO
             ELSE
                CALL outwin(A,B,NODES,UV,DUV,rel,vr,r0(jatom),dx(jatom),jri(jatom),ei,fl,zz(jatom),nrad)
             ENDIF
             !                                                                       
             OVLP = rint13(REL,A,B,A,B,nrad,dx(jatom),jri(jatom),r0(jatom)) 
             TRX=1.0d0/SQRT(OVLP)
             IMAX=jri(JATOM)
             rf1(1:IMAX,l,is,irf) = A(1:IMAX)*TRX
             rf2(1:IMAX,l,is,irf) = B(1:IMAX)*TRX
             P(l,is,irf)  = UV*TRX
             DP(l,is,irf) = DUV*TRX
             
             plo = p(l,is,irf)
             dplo=dp(l,is,irf)
             pi12lo = rint13(REL,rf1(:,l,is,1),rf2(:,l,is,1),rf1(:,l,is,irf),rf2(:,l,is,irf),nrad,dx(jatom),jri(jatom),r0(jatom)) 
             pe12lo = rint13(REL,rf1(:,l,is,2),rf2(:,l,is,2),rf1(:,l,is,irf),rf2(:,l,is,irf),nrad,dx(jatom),jri(jatom),r0(jatom))
             
             qprint=myrank.EQ.master
             !if (myrank.EQ.master) WRITE(6,*)
             alo(l,jlo,:,is)=0.d0
             call abc0(alo(l,jlo,1,is),alo(l,jlo,2,is),alo(l,jlo,irf,is),l,jlo,lapw(l),P(l,is,1),DP(l,is,1),P(l,is,2),DP(l,is,2),pei(l),plo,dplo,pi12lo,pe12lo,nrad,lomax,nloat,lmax2,Rmt(jatom),qprint)
          else
             ! For APW+lo we do not need any local orbital. Instead alo(l,1,1:2,is) are determined by u_{new}=a^lo u + b^lo dotu.
             qprint=myrank.EQ.master
             !if (myrank.EQ.master) WRITE(6,*)
             alo(l,jlo,:,is)=0.d0
             call abc0(alo(l,jlo,1,is),alo(l,jlo,2,is),alo(l,jlo,2+jlo,is),l,jlo,lapw(l),P(l,is,1),DP(l,is,1),P(l,is,2),DP(l,is,2),pei(l),plo,dplo,pi12lo,pe12lo,nrad,lomax,nloat,lmax2,Rmt(jatom),qprint)
          endif
       enddo
    ENDDO ! 170  continue
    !... CALCULATION OF RADIAL INTEGRALS
    IF (.NOT.(ISPIN.EQ.2.AND.IS.EQ.1)) THEN !  GOTO 469
       CALL RADINT(JATOM,ISPIN)
    ENDIF
    !469 
    RETURN                                                            
    
2022 FORMAT(3X,4E19.12)
5   FORMAT(10X,' ENERGY PARAMETERS ARE',7F12.6)                        
7   FORMAT(/10X,'ATOMIC PARAMETERS FOR ',A10/)                        
14  FORMAT(/11X,1HL,5X,4HU(R),10X,5HU'(R),9X,5HDU/DE,8X,6HDU'/DE,6X,7HNORM-U')
8   FORMAT(10X,I2,5E14.6,5X,3I2)
1980 FORMAT(3X)
2000 FORMAT(16X,I2//)
2030 FORMAT(///)
2031 FORMAT(/)
  END SUBROUTINE ATPAR


  SUBROUTINE ROUT(L1,iso)
    USE param
    !IMPLICIT REAL*8 (A-H,O-Z)
    IMPLICIT NONE
    INTEGER, intent(in) :: l1, iso
    !
    INTEGER :: irf1, irf2, is1
    write(6,*)'L1=',L1 
    WRITE(6,*)'RADIAL INTEGRALS:'
    DO IRF1=1,NRF
       DO IS1=1,iso
          DO IRF2=1,NRF
             if (abs(RI_MAT(IRF1,IRF2,L1,IS1)).gt.1d-4) WRITE(6,*)irf1,irf2,'     ',is1,RI_MAT(IRF1,IRF2,L1,IS1)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE ROUT
  
  SUBROUTINE RADINT(JATOM,ISPIN)
    USE param
    USE com_mpi
    USE structure
    USE case
    IMPLICIT NONE
    INTEGER, intent(in) :: jatom, ispin
    interface
       Function RINT13(REL,A,B,X,Y,NRAD,DX_,JRI_,R0_) result(S) ! Calculates overlap between psi_1=(A,B) and psi_2=(X,Y) functions
         REAL*8 :: S
         LOGICAL, intent(in) :: REL   ! relativistic or not
         REAL*8, intent(in)  :: A(NRAD),B(NRAD),X(NRAD),Y(NRAD)
         INTEGER, intent(in) :: NRAD, JRI_
         REAL*8, intent(in)  :: DX_, R0_
       End Function RINT13
    end interface
    INTEGER :: if1, if2, is1, l
    ! rf1 and rf2 are two parts of the solution of the Dirac equation : large and small component
    ! rf1 has the following indexes:
    !    rf1(r,l,s,irf)  where r is radial distance, s is spin,
    !                    and irf is one of the radial functions 
    !    for LAPW      we have u(E1),dot{u}(E1)
    !    for LAPW+LO   we have u(E1),dot{u}(E1), u(E2)
    !    for APW+lo    we have u(E1),dot{u}(E1)
    !    for APW+lo+LO we have u(E1),dot{u}(E1), u(E2)
    !
    ri_mat(1:nrf,1:nrf,0:lmax2,1:2)=0.d0
    do l=0,3
       do  IS1=1,ISPIN
          do  IF1=1,NRF
             do  IF2=1,NRF
                !CALL RINT13(rf1(1,l,is1,if1),rf2(1,l,is1,if1), rf1(1,l,is1,if2),rf2(1,l,is1,if2), ri_mat(l,l,if1,if2,is1),JATOM)
                ri_mat(if1,if2,l,is1) = rint13(REL,rf1(:,l,is1,if1),rf2(:,l,is1,if1),rf1(:,l,is1,if2),rf2(:,l,is1,if2),nrad,dx(jatom),jri(jatom),r0(jatom))
             enddo
          enddo
       enddo
       if (myrank.EQ.master) call rout(l,ispin)
    enddo
  END SUBROUTINE RADINT
  
end module matpar
