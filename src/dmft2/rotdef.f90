SUBROUTINE ROTDEF (iz,tau,iord,nat,pos,ndif,rotij,tauij,mult,lattic)
  !                                                                       
  !     THE MATRICES  ROTIJ(3,3,INDEX)  TRANSFORM THE POSITION OF AN      
  !     ATOM TO ITS CORRESPONDING POSITION OF AN EQUIVALENT ATOM.
  !               
  use structure, only: aa,bb,cc
  use com_mpi,   ONLY: myrank, master
  IMPLICIT NONE
  INTEGER,INTENT(IN)     :: iord,nat,ndif
  INTEGER,INTENT(IN)     :: iz(3,3,iord),mult(nat)
  CHARACTER*4,INTENT(IN) :: lattic
  REAL(8),INTENT(IN)     :: tau(3,iord)
  REAL(8)                :: pos(3,ndif)
  REAL(8),INTENT(OUT)    :: rotij(3,3,ndif),tauij(3,ndif)
  ! locals
  CHARACTER*67    :: ERRMSG
  REAL(8)         :: x(3),x1(3),toler,toler2
  INTEGER         :: i,i1,j,m,index,index1,jatom,ncount
  LOGICAL         :: FOUND
  !
  ! Additions for averaging
  integer         :: jtimes,nthis
  real*8          :: X00, Y00, Z00, X0, Y0, Z0,XA,YA,ZA
  real*8          :: XB, YB, ZB, TT
  !
  toler = 1.D-7
  toler2=1.5d0*toler            
  !
  !---------------------------------------------------------------------  
  !  Patch up symmetry issues with positions in limited precision
  DO JTIMES=1,2
     !      write(88,*)'Pass ',Jtimes
     DO JATOM=1,ndif
        NTHIS=0
        X00=POS(1,JATOM)
        Y00=POS(2,JATOM)
        Z00=POS(3,JATOM)
        X0=0
        Y0=0
        Z0=0                                           
        DO I=1,IORD                                              
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
           TT=1.D0/dble(nthis)                         
           POS(1,JATOM)=X0*TT
           POS(2,JATOM)=Y0*TT
           POS(3,JATOM)=Z0*TT
           !       write(88,88)'OLD     ',X00,Y00,Z00
           !       write(88,88)'Average ',X0*TT,Y0*TT,Z0*TT,NTHIS,JATOM
        endif
        !88      format(a,3f16.13,2i3)
     ENDDO
  ENDDO
  !---------------------------------------------------------------------  
  !
  INDEX=0                                                           
  NCOUNT=0                                                          
  DO JATOM=1,NAT                                                 
     INDEX1=INDEX+1                                                 
     DO M=1,MULT(JATOM)                                          
        INDEX=INDEX+1
        FOUND = .FALSE.
        DO I=1,IORD                                              
           x(1:3)=0.d0
           x=MATMUL(TRANSPOSE(iz(1:3,1:3,i)),pos(1:3,index1))
           x(1:3)=x(1:3)+tau(1:3,i)
           x(1) = MOD(x(1) + 1.d0, 1.d0)
           x(2) = MOD(x(2) + 1.d0, 1.d0)
           x(3) = MOD(x(3) + 1.d0, 1.d0)
           x1(1:3)=MOD(ABS(X(1:3)-POS(1:3,INDEX))+toler,1.d0)-toler
           !           WRITE(*,*) 'JATOM,INDEX,I',JATOM,INDEX,I                   
           !           WRITE(*,*) ABS(X1(1:3)),toler
           IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
              NCOUNT=NCOUNT+1                                             
              TAUIJ(1:3,INDEX)=TAU(1:3,I)
              ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
              FOUND = .TRUE.
              EXIT ! GOTO 30                                                     
           END IF
           !....check positions for centered lattices
           if(lattic(1:1).eq.'B') then
              x1(1:3)=mod(x1(1:3)+0.5d0+toler,1.d0)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 TAUIJ(1:3,INDEX)=TAU(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
                 FOUND = .TRUE.
                 EXIT !GOTO 30                                                     
              END IF
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXY') then
              x1(1:2)=mod(x1(1:2)+0.5d0+toler,1.d0)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 TAUIJ(1:3,INDEX)=TAU(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
                 FOUND = .TRUE.
                 EXIT !GOTO 30                                                     
              END IF
              x1(1:2)=mod(x1(1:2)+0.5d0,1.d0)
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXZ') then
              x1(1)=mod(x1(1)+0.5d0+toler,1.d0)-toler
              x1(3)=mod(x1(3)+0.5d0+toler,1.d0)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 TAUIJ(1:3,INDEX)=TAU(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
                 FOUND = .TRUE.
                 EXIT ! GOTO 30                                                     
              END IF
              x1(1)=mod(x1(1)+0.5d0,1.d0)
              x1(3)=mod(x1(3)+0.5d0,1.d0)
           endif
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CYZ') then
              x1(2:3)=mod(x1(2:3)+0.5d0+toler,1.d0)-toler
              IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                 NCOUNT=NCOUNT+1                                             
                 TAUIJ(1:3,INDEX)=TAU(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
                 FOUND = .TRUE.
                 EXIT !GOTO 30                                                     
              END IF
           end if
        ENDDO
        IF (.NOT.FOUND) THEN
           !
           !           Error: no symmetry operation found
           !
           CALL OUTERR('ROTDEF','no symmetry operation found.')
           WRITE(ERRMSG,9000) jatom, index          
           CALL OUTERR('ROTDEF',ERRMSG)
           WRITE(ERRMSG,9010) (POS(I1,index1),I1=1,3) 
           CALL OUTERR('ROTDEF',ERRMSG)
           WRITE(ERRMSG,9020) (POS(I1,INDEX),I1=1,3) 
           CALL OUTERR('ROTDEF',ERRMSG)
           STOP 'ROTDEF - Error'
       ELSE
          if (myrank.EQ.master) WRITE(6,'(A,I3,1x,A,I3,1x,A,I3,1x)') 'Operation', I, 'transforms atom ', index, 'from first atom ', index1
        ENDIF
     ENDDO
  ENDDO
  IF(NCOUNT.NE.INDEX) THEN
     CALL OUTERR('ROTDEF','ROTIJ not defined for all atoms of basis.')
     WRITE (ERRMSG,9030) NCOUNT
     CALL OUTERR('ROTDEF',ERRMSG)
     STOP 'ROTDEF - Error'
  ENDIF

  RETURN     
  
9000 FORMAT ('for jatom, index',I2,i2)
9010 FORMAT ('atomposition of jatom',3F12.7)
9020 FORMAT ('atomposition of index',3F12.7)
9030 FORMAT ('NCOUNT=',I2)
END SUBROUTINE ROTDEF
