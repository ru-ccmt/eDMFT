SUBROUTINE get_rotij_tauij(rotij,tauij,pos,alat,iz,tau,iord,nat,ndf,mult,lattic)
  !     THE MATRICES  ROTIJ(3,3,INDEX)  TRANSFORM THE POSITION OF AN      
  !     ATOM TO ITS CORRESPONDING POSITION OF AN EQUIVALENT ATOM.
  IMPLICIT NONE
  !               
  REAL(8),   intent(out) :: rotij(3,3,ndf),tauij(3,ndf)
  REAL(8), intent(inout) :: pos(3,ndf)
  REAL*8,     intent(in) :: alat(3)
  INTEGER,    intent(IN) :: iord,nat,ndf
  INTEGER,    intent(IN) :: iz(3,3,iord),mult(nat)
  CHARACTER*4,intent(IN) :: lattic
  REAL(8),    intent(IN) :: tau(3,iord)
  ! locals
  REAL*8        :: aa,bb,cc
  CHARACTER*67  :: ERRMSG
  REAL(8)       :: x(3),x1(3),toler,toler2,one
  INTEGER       :: i,i1,j,m,jatom,ncount!,index,index1,k, l
  DATA TOLER/1.D-7/,ONE/1.D0/
  !
  ! Additions for averaging
  integer :: jtimes,nthis, latom, lfirst
  real*8  :: X00, Y00, Z00, X0, Y0, Z0,XA,YA,ZA
  real*8  :: XB, YB, ZB, TT
  logical :: found
  !
  toler2=1.5d0*toler
  !
  rotij(1:3,1:3,1:ndf)=0.0D0
  tauij(1:3,1:ndf)=0.0D0
  !
  !---------------------------------------------------------------------  
  !  Patch up symmetry issues with positions in limited precision
  aa=alat(1)
  bb=alat(2)
  cc=alat(3)
  !index = sum(mult)
  ! In this loop we increase precision of atom positions POS(:) by averaging over all possible symmetry related positions,
  ! which should give exactly the same position, but due to numerics are slihtly different.
  DO JTIMES=1,2
     !      write(88,*)'Pass ',Jtimes
     DO JATOM=1,sum(mult)
        NTHIS=0
        X00=POS(1,JATOM)
        Y00=POS(2,JATOM)
        Z00=POS(3,JATOM)
        X0=0
        Y0=0
        Z0=0                                           
        DO I=1,IORD
!!! **** maybe more readable implementation **** !!!
           !! Ra(:) is obtained by applyin symmetry operation on atom position POS(:,jatom)
           !RA(:) = matmul(pos(:,jatom),IZ(:,:,i)) + TAU(:,i)
           !RA(1) = mod(RA(1)+1.d0, 1.d0)
           !RA(2) = mod(RA(2)+1.d0, 1.d0)
           !RA(3) = mod(RA(3)+1.d0, 1.d0)
           !do j=1,3
           !   if (RA(j).gt..0.999999) RA(j) = RA(j) -1.D0
           !enddo
           !RB(:) = abs(RA(:)-POS(:,jatom))*alat(:)
           !if (RB(1).lt.1d-2 .and. RB(2).lt.1d-2 .and. RB(3).lt.1d-2) then
           !   ! We found symmetry operation i, which does not change POS, i.e., leaves POS invariant
           !   NTHIS=NTHIS+1
           !   R0(:) = R0(:) + RA(:)
           !endif
           
           XA=TAU(1,I)+1.D0
           YA=TAU(2,I)+1.D0
           ZA=TAU(3,I)+1.D0                                                      
           DO J=1,3 
              XA=XA+dble(IZ(J,1,I))*POS(J,JATOM)                               
              YA=YA+dble(IZ(J,2,I))*POS(J,JATOM) 
              ZA=ZA+dble(IZ(J,3,I))*POS(J,JATOM)
           ENDDO
           XA = mod(XA,1.D0)
           YA = mod(YA,1.D0)
           ZA = mod(ZA,1.D0)
           if(xa.gt.0.999999)XA=XA-1.D0
           if(ya.gt.0.999999)YA=YA-1.D0
           if(Za.gt.0.999999)ZA=ZA-1.D0
           XB=ABS(XA-X00)*aa                                      
           YB=ABS(YA-Y00)*bb                 
           ZB=ABS(ZA-Z00)*cc                                      
           
           IF(XB.LT.1d-2.AND.YB.LT.1d-2.AND.ZB.LT.1d-2) THEN        
              NTHIS=NTHIS+1                                             
              X0=X0+XA
              Y0=Y0+YA
              Z0=Z0+ZA
           ENDIF
        ENDDO
        if(nthis.gt.1)then
           ! To increase precision, we average over all possible symmetry related positions, which should give exactly the same
           ! position, but due to numerics are slihtly different.
           TT=1.D0/dble(nthis)                         
           POS(1,JATOM)=X0*TT
           POS(2,JATOM)=Y0*TT
           POS(3,JATOM)=Z0*TT
        endif
!88      format(a,3f16.13,2i3)
     ENDDO
  ENDDO
  !---------------------------------------------------------------------  
  !
  !INDEX=0                                                           
  NCOUNT=0
  DO JATOM=1,NAT
     lfirst = 1+sum(mult(1:jatom-1))
     !INDEX1=INDEX+1
     DO M=1,MULT(JATOM)
        latom = lfirst + m -1
        !INDEX=INDEX+1
        !print *, 'lfirst-index1=', lfirst-index1
        !print *, 'latom-index=', latom-index
        
        found=.false.
        DO i=1,iord
           ! X <- Gamma * pos(:,lfirst)
           X(:) = matmul( transpose(iz(:,:,i)), pos(:,lfirst) )
           X(:)=X(:)+tau(:,i)
           ! X1 is very small if pos(:,latom) = Gamma * pos(:,lfirst)
           X1(:) = mod( abs(X(:)-POS(:,latom))+toler, 1.d0)-toler
           if ( maxval(abs(X1)).LT.TOLER2) then
              NCOUNT=NCOUNT+1
              TAUIJ(:,latom)=TAU(:,I)
              ROTIJ(:,:,latom)=IZ(:,:,I)
              found=.true.
              EXIT
           endif
           !....check positions for centered lattices
           if (lattic(1:1).eq.'B') then
              x1(:) = mod( x1(:)+0.5d0+toler, 1.d0)-toler
              if ( maxval(abs(X1)).LT.TOLER2) THEN
                 NCOUNT=NCOUNT+1                                             
                 TAUIJ(1:3,latom)=TAU(1:3,I)
                 ROTIJ(1:3,1:3,latom)=IZ(1:3,1:3,I)
                 found=.true.
                 EXIT
              END IF
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXY') then
              x1(1:2) = mod(x1(1:2)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN
                 NCOUNT=NCOUNT+1                                             
                 TAUIJ(1:3,latom)=TAU(1:3,I)
                 ROTIJ(1:3,1:3,latom)=IZ(1:3,1:3,I)
                 found=.true.
                 EXIT
              END IF
              x1(1:2)=mod(x1(1:2)+0.5d0,one)
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXZ') then
              x1(1)=mod(x1(1)+0.5d0+toler,one)-toler
              x1(3)=mod(x1(3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 TAUIJ(1:3,latom)=TAU(1:3,I)
                 ROTIJ(1:3,1:3,latom)=IZ(1:3,1:3,I)
                 found=.true.
                 EXIT
              END IF
              x1(1)=mod(x1(1)+0.5d0,one)
              x1(3)=mod(x1(3)+0.5d0,one)
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CYZ') then
              x1(2:3)=mod(x1(2:3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 TAUIJ(1:3,latom)=TAU(1:3,I)
                 ROTIJ(1:3,1:3,latom)=IZ(1:3,1:3,I)
                 found=.true.
                 EXIT
              END IF
           end if
        ENDDO
        if (.not. found) then
           !           Error: no symmetry operation found
           CALL OUTERR('ROTDEF','no symmetry operation found.')
           WRITE(ERRMSG,9000) jatom, latom          
           CALL OUTERR('ROTDEF',ERRMSG)
           WRITE(ERRMSG,9010) (POS(I1,lfirst),I1=1,3) 
           CALL OUTERR('ROTDEF',ERRMSG)
           WRITE(ERRMSG,9020) (POS(I1,latom),I1=1,3) 
           CALL OUTERR('ROTDEF',ERRMSG)
           STOP 'ROTDEF - Error'
        end if
     ENDDO
  ENDDO
  IF(NCOUNT.NE.sum(mult)) THEN
     CALL OUTERR('ROTDEF','ROTIJ not defined for all atoms of basis.')
     WRITE (ERRMSG,9030) NCOUNT
     CALL OUTERR('ROTDEF',ERRMSG)
     STOP 'ROTDEF - Error'
  ENDIF
  RETURN     
  !
9000 FORMAT ('for jatom, index',I2,i2)
9010 FORMAT ('atomposition of jatom',3F12.7)
9020 FORMAT ('atomposition of index',3F12.7)
9030 FORMAT ('NCOUNT=',I2)
END SUBROUTINE get_rotij_tauij
