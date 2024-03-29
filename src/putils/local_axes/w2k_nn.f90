!program nn
!  REAL*8 :: dlimit, dstmax, dfac
!  CHARACTER*80 :: case
!  
!  dlimit=0.d0 
!  dstmax=-1   
!  dfac=2      
!  case='Sr2IrO4'
!  call w2k_nn(case,dfac,dlimit,dstmax)
!  
!end program nn


module variable_fields
  real*8, allocatable :: pos(:,:)
end module variable_fields

SUBROUTINE w2knn(case,dfac_,dlimit_,dstmax_)
  use variable_fields
  IMPLICIT REAL*8 (A-H,O-Z)
  !     Constant parameter definition
  PARAMETER (NNN=100000)                                              
  PARAMETER (mshell= 100)
  !
  CHARACTER*80, intent(in) :: case
  REAL*8, intent(in)       :: dfac_, dlimit_,dstmax_
  !
  !f2py real*8  optional,intent(in)        :: dfac_ = 2.0
  !f2py real*8  optional,intent(in)        :: dlimit_ = 0.0
  !f2py real*8  optional,intent(in)        :: dstmax_ = -1
  !
  !     NATO    NUMBER OF  NONEQUIVALENT ATOMS            
  !     NNN     NUMBER OF NEIGHBOURS                                      
  !     NDIF    NUMBER OF ALL ATOMS IN UNITCELL                           
  !     MSHELL  Max. Nr. of SHELLS
  !
  ! from struk
  real*8            BR2(3,3)
  real*8            BR1(3,3)
  CHARACTER*4       LATTIC                                            
  CHARACTER*10, allocatable :: name(:)
  ! 
  LOGICAL           ORTHO,error
  LOGICAL           SHOWALL
  CHARACTER*10      KNAME
  CHARACTER*11      STATUS,FORM                                     
  CHARACTER*79      TITLE                                           
  CHARACTER*100     FNAME,fname1                                           
  !-----------------------------------------------------------------------
  real*8  A(3),PIA(3),VI   
  DIMENSION DIF(3),XX(3),PP(3),P(3),DISTS(NNN),NR(NNN),PNN(3,NNN)   
  DIMENSION NNAT(NNN),help(3)                                            
  
  LOGICAL, allocatable :: ovlap(:),used(:)
  real*8,  allocatable :: rmt(:),zz(:),shdist(:,:)
  real*8,  allocatable :: zzorig(:),shellz(:,:,:)
  real*8,  allocatable :: V(:),POSN(:,:)
  real*8,  allocatable :: R0(:),R0N(:),DH(:)  
  real*8,  allocatable :: rmtt(:),rmttn(:),zzn(:),zzo(:) 
  real*8,  allocatable :: ROTLOC(:,:,:),ROTIJ(:,:,:),TAUIJ(:,:)  
  integer,  allocatable :: JRJ(:),JRJN(:)
  integer,  allocatable :: IATNR(:),MULT(:),icnt(:,:,:) 
  integer,  allocatable :: iz(:,:),ishell(:),ityp(:),imult(:)
  CHARACTER*10, allocatable :: namen(:)
  CHARACTER*100 :: f_20, f_66, f_67, f_21, f_68, f_06
  !-----------------------------------------------------------------------
  logical there
  dfac=dfac_
  dlimit=dlimit_
  dstmax=dstmax_
  
  ! 20,'Sr2IrO4.struct',   'old',    'formatted',0
  ! 66,'Sr2IrO4.outputnn', 'unknown','formatted',0
  ! 67,'Sr2IrO4.nnshells', 'unknown','formatted',0
  ! 68,'Sr2IrO4.rotlm_',   'unknown','formatted',0
  ! 21,'Sr2IrO4.struct_nn','unknown','formatted',0
  f_20 = TRIM(ADJUSTL(case))//'.struct'
  f_66 = TRIM(ADJUSTL(case))//'.outputnn_'
  f_67 = TRIM(ADJUSTL(case))//'.nnshells_'
  f_21 = TRIM(ADJUSTL(case))//'.struct_nn'
  f_68 = TRIM(ADJUSTL(case))//'.rotlm_'
  f_06 = TRIM(ADJUSTL(case))//'.info_'

  WRITE(6,*) 'structure file=', f_90
  
  open(20, FILE=f_20, STATUS='old', FORM='formatted')
  open(66, FILE=f_66, STATUS='unknown', FORM='formatted')
  open(67, FILE=f_67, STATUS='unknown', FORM='formatted')

  open(6, FILE=f_06, STATUS='unknown', FORM='formatted')
  
  fname1 = f_20
  
  !     This should be how much larger (in 1D) the DFT distances are relative to true
  DFTERR = 1.01D0        !This is about right for PBE -- could look in case.in0?
  !
  !     See if the user has calibrated the lattice parameter scaling difference
  !inquire(file='.latcalib',exist=there)
  !if(there)then
  !   open(unit=99,file='.latcalib')
  !   !       This should be how much larger (in 1D) the DFT distances are relative to true
  !   read(99,*)DFTERR
  !   close(unit=99)
  !endif
  SHOWALL=.TRUE.
  
  if(dfac .lt. 0.)then
     SHOWALL=.FALSE.
     dfac=abs(dfac)
  endif
  if(dlimit.lt.1.d-7) dlimit=1.d-5
  ishellmax=99
  
  RTB2=SQRT(3.)/2.                                                  
  !                                                                       
  !                                                                       
  !.....START READING FILE STRUCT, writing new struct (21)                   
  READ(20,1510) TITLE                                               
  write(66,1510) TITLE
  READ(20,1511) LATTIC,NAT,title                                          
  write(66,1511) LATTIC,NAT,title                                          
  !    allocate nato-arrays
  allocate (  name(nat) )
  allocate ( ovlap(nat) )
  allocate ( rmt(nat),zz(nat) )
  allocate ( zzorig(nat),V(nat),R0(nat),DH(nat),rmtt(nat) )
  allocate ( ROTLOC(3,3,nat),JRJ(nat),IATNR(nat),MULT(nat) ) 
  allocate ( POS(3,48*nat*2) )
  !     READ IN LATTICE CONSTANTS                                         
  read(20,'(6F10.6)') A(1),A(2),A(3),alpha,beta,gamma
  if(alpha.eq.0.d0) alpha=90.d0                                       
  if(beta .eq.0.d0) beta =90.d0                                       
  if(gamma.eq.0.d0) gamma=90.d0                                       
  write(66,'(6F10.6)') A(1),A(2),A(3),alpha,beta,gamma
  !17 FORMAT(6F10.6)                                                    
  if(dstmax .lt. 0)dstmax=max(20.d0,0.61d0*max(a(1),a(2),a(3)))  
  if(nat .eq. 1) then
     dstmax=max(20.d0,a(1),a(2),a(3))*1.1
  else
     dstmax=min(40.d0,dstmax)
  endif

  WRITE(6,*) 'A[1]=', A(1), 'A[2]=', A(2), 'A[3]=', A(3)

  WRITE(6,*) 'DSTMAX:',dstmax                  
  iix=5
  iiy=5
  iiz=5
  if(lattic(1:1).eq.'P'.or.lattic(1:1).eq.'H') then
     iix=max(1,nint(dstmax/a(1)))+1
     iiy=max(1,nint(dstmax/a(2)))+1
     iiz=max(1,nint(dstmax/a(3)))+1
  endif
  
  DO i=1,1000
     if(nat .gt. 1)then
        if( (iix-1)*a(1) .gt. dstmax*2.d0)then
           iix=max(1,iix-1)
        else if( (iiy-1)*a(2) .gt. dstmax*2.d0)then
           iiy=max(1,iiy-1)
        else if( (iiz-1)*a(3) .gt. dstmax*2.d0)then
           iiz=max(1,iiz-1)
           EXIT
        endif
     endif
  ENDDO
  
  if(lattic(1:3).eq.'CXY') then     !fix for orthorombic CXY with very different a,b
     iix=max(iix,iiy)
     iiy=iix
  endif
  WRITE(6,*) 'iix,iiy,iiz',iix,iiy,iiz,a(1)*iix,a(2)*iiy,a(3)*iiz

  !     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
  !     NONEQUIVALENT ATOMS                                               
  INDEX=0                                                           
  DO JATOM = 1,NAT        ! 50
     INDEX=INDEX+1                                                  
     READ(20,1012) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM)  
     write(66,1012) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM)  
     IF (MULT(JATOM).EQ.0) THEN                                     
        write(66,1020) JATOM,INDEX,MULT(JATOM)                       
        STOP ' NNN: MULT EQ 0'                                      
     ENDIF
     DO M=1,MULT(JATOM)-1 
        INDEX=INDEX+1                                            
        READ(20,1011) IATNR(JATOM),( POS(J,INDEX),J=1,3)         
        write(66,1011) IATNR(JATOM),( POS(J,INDEX),J=1,3)         
     ENDDO

     READ(20,1050) NAME(JATOM),JRJ(JATOM),R0(JATOM),RMTT(jatom),  zz(jatom)  
     zzorig(jatom)=zz(jatom)      
     if(name(jatom)(3:3).ne.' ') then
        write(66,*) 'NAMED ATOM: ',name(jatom), 'Z changed to ', 'IATNR+999 to determine equivalency'
        write(*,*) 'NAMED ATOM: ',name(jatom), 'Z changed to ', 'IATNR+999 to determine equivalency'
        zz(jatom)=999+jatom
     endif
     write(66,1049) NAME(JATOM),JRJ(JATOM),R0(JATOM),RMTT(jatom), zzorig(jatom)        
     if((jrj(jatom)/2)*2.eq.jrj(jatom)) then
        write(*,*) 'WARNING: JRJ of atom',jatom,' is even:',jrj(jatom)
        write(*,*) 'CHANGE it to ODD number !!!!'
        write(66,*) 'WARNING: JRJ of atom',jatom,' is even:', jrj(jatom)
        write(66,*) 'CHANGE it to ODD number !!!!'
     endif
     DH(JATOM)=LOG(RMTT(jatom)/R0(JATOM)) / (JRJ(JATOM)-1)                 
     RMT(JATOM)=R0(JATOM)*EXP( DH(JATOM)*(JRJ(JATOM)-1) )           
     READ(20,1051) ((ROTLOC(I1,I2,JATOM),I1=1,3),I2=1,3)               
     write(66,1051) ((ROTLOC(I1,I2,JATOM),I1=1,3),I2=1,3)               
  ENDDO !50   CONTINUE                                                          
  
  call reduce_alloc_r8_2d(3,48*nat,3,index)
  allocate ( shdist(index,mshell),shellz(index,mshell,48) )
  allocate ( namen(index),POSN(3,index),R0N(index),rmttn(index) )
  allocate ( zzn(index),zzo(index),ROTIJ(3,3,index),TAUIJ(3,index) )
  allocate ( JRJN(index),icnt(index,mshell,48),iz(index,mshell) )   
  allocate ( used(nat*index),ishell(nat*index),ityp(nat*index) )
  allocate ( imult(nat*index) )
  !                                                                       
  !                                                                       
  !.....SET UP LATTICE, AND IF REQUIRED ROTATION MATRICES                 
  CALL DIRLAT (BR1,BR2,ortho,lattic,nat,alpha,beta,gamma,A)
  !                                           
  write(66,*) ' '
  write(66,*) 'Bond-Valence Sums are calculated for current lattice parameters'
  write(66,*) 'and rescaled ones by 1 %. (You can put scaling into  .latcalib)' 
  pi=4.d0*atan(1.d0)
  cosgam=cos(gamma/180.d0*pi)
  singam=sin(gamma/180.d0*pi)           
  INDEX=0                                                           
  DO JATOM=1,NAT       ! 200
     do i=1,nat
        ovlap(i)=.true.
     enddo
     DO M=1,MULT(JATOM)  ! 190
        INDEX=INDEX+1                                                     
        DO J=1,3                                                      
           XX(J)=POS(J,INDEX)
        ENDDO
if(SHOWALL.or.(M.EQ.1)) write(66,'(/,A,I3,2X,A,I3,2X,A10,A,3F10.5)') ' ATOM:',JATOM,'EQUIV.',M,NAME(JATOM),' AT',XX(1),XX(2),XX(3)
        NC=0             
        DO I1=-iix,iix              ! 180
           DO I2=-iiy,iiy           ! 180
              DO I3=-iiz,iiz        ! 180
                 IF(ortho) THEN                                   
                    P(1)=I1*BR2(1,1)+I2*BR2(2,1)+I3*BR2(3,1)                    
                    P(2)=I1*BR2(1,2)+I2*BR2(2,2)+I3*BR2(3,2)                    
                    P(3)=I1*BR2(1,3)+I2*BR2(2,3)+I3*BR2(3,3)                    
                 ELSE                                                        
                    P(1)=I1                                                     
                    P(2)=I2                                                     
                    P(3)=I3                                                     
                    IF(LATTIC(1:3).eq.'CXZ') THEN
                       P(1)=I1*0.5d0+i3*0.5d0
                       P(2)=I2
                       P(3)=-I1*0.5d0+i3*0.5d0
                    END IF
                 ENDIF
                 K=0                                                           
                 DO JAT=1,NAT           ! 120
                    DO MM=1,MULT(JAT)   ! 110
                       K=K+1                                                             
                       DIST=0.                                                           
                       DO L=1,3                                                      
                          PP(L)=POS(L,K)+P(L)
                          DIF(L)=XX(L)-PP(L)
                       ENDDO
                       
                       IF (.not.ortho) THEN                                       
                          help(1)=dif(1)  
                          help(2)=dif(2)  
                          help(3)=dif(3)  
                          if(lattic(1:1).eq.'R') then
                             dif(1)=help(1)*BR2(1,1)+help(2)*BR2(2,1)+help(3)*BR2(3,1)             
                             dif(2)=help(1)*BR2(1,2)+help(2)*BR2(2,2)+help(3)*BR2(3,2)             
                             dif(3)=help(1)*BR2(1,3)+help(2)*BR2(2,3)+help(3)*BR2(3,3)           
                          elseif(lattic(1:3).eq.'CXZ') then
                             dif(1)=help(1)*singam            
                             dif(2)=(help(1)*cosgam*a(1)+help(2)*a(2))/a(2)             
                             dif(3)=help(3)           
                          else
                             dif(1)=(help(1)*BR2(1,1)*a(1)+help(2)*BR2(2,1)*a(2)+help(3)*BR2(3,1)*a(3))/a(1)
                             dif(2)=(help(1)*BR2(1,2)*a(1)+help(2)*BR2(2,2)*a(2)+help(3)*BR2(3,2)*a(3))/a(2)
                             dif(3)=(help(1)*BR2(1,3)*a(1)+help(2)*BR2(2,3)*a(2)+help(3)*BR2(3,3)*a(3))/a(3)
                          endif
                       ENDIF
                       DO L=1,3                                                      
                          DIST=DIST+DIF(L)*DIF(L)*A(L)*A(L)
                       ENDDO
                       DIST=SQRT(DIST)
                       
                       IF(DIST.GT.dstmax) CYCLE
                       IF(DIST.LT..001) CYCLE
                       NC=NC+1         
                       if(nc.gt.nnn) stop ' nnn too small'                                     
                       DISTS(NC)=DIST                                                    
                       NNAT(NC)=JAT                                                      
                       DO L=1,3                                                      
                          PNN(L,NC)=PP(L)
                       ENDDO
                    ENDDO !110 CONTINUE
                 ENDDO !120 CONTINUE                                                          
              ENDDO !  180 CONTINUE
           ENDDO    ! 180
        ENDDO       ! 180
        CALL ORD2(DISTS,NR,NC)                                            
        N1=1                                                              
        N2=NR(N1)                                                         
        N3=NNAT(N2)             
        SUMRAD=RMT(JATOM)+RMT(N3)                                         
        IF(M.EQ.1) THEN                                                   
           IF(SUMRAD.LE.DISTS(N1))THEN                                       
              WRITE(*,'(/,3X,A,I3,2X,A10,A,I3,2X,A10)') ' ATOM',JATOM,NAME(JATOM),' ATOM',N3,NAME(N3)
    WRITE(*, '(A,I3,A,F7.5,A,I3,A,F7.5,/,A,F8.5,2X,A,F8.5)') ' RMT(',JATOM,')=',RMT(JATOM),' AND RMT(',N3,')=',RMT(N3),' &
      & SUMS TO',SUMRAD,'LT.  NN-DIST=',DISTS(1)
              WRITE(66,'(A,I3,A,F7.5,A,I3,A,F7.5,/,A,F8.5,2X,A,F8.5)') ' RMT(',JATOM,')=',RMT(JATOM),' & 
         &AND RMT(',N3,')=',RMT(N3),' SUMS TO',SUMRAD,'LT.  NN-DIST=',DISTS(1)
           ELSE                                                              
              ovlap(n3)=.false.                        
              write(66,4) JATOM,RMT(JATOM),N3,RMT(N3), SUMRAD, DISTS(N1)
4             FORMAT(/,'   ERROR !!!!!!!!!!!!!!!', /,' RMT(',I3,')=',F7.5,' AND RMT(',I3,')=',F7.5,/,' &
              &SUMS TO',F8.5,' GT NNN-DIST=',F8.5)
              WRITE(*,4) JATOM,RMT(JATOM),N3,RMT(N3), SUMRAD, DISTS(N1)         
           ENDIF
        ENDIF
        !.....determination of "equal" atoms
        !
        olddis=0.d0
        ishell(index)=0
        !
        bva=0.0
        bvaerr=0.0

        !!!! XX(:) = POS(:,INDEX)
        !!!! N1=1,...NC
        !!!! DISTS(N1) -- distance
        !!!! NNAT(NC)  -- isort
        !!!! PNN(:3, N2=NR(N1) ) -- neigbor
        
        DO N1=1,NC
           N2=NR(N1)                                                         
           N3=NNAT(N2)                                                       
           !
           !     Include BVA
           i1=nint(zzorig(jatom))
           i2=nint(zzorig(N3))

           !! You can restore this setting this statement to .true.
           if (.false.) then
              call bvan(i1,i2,scale,dbva)
              if(scale .gt. 0)then
                 val=exp( (dbva - dists(n1)*0.529177) /scale)
                 bva=bva+val
                 valerr=exp( (dbva - dists(n1)*0.529177/DFTERR) /scale)
                 bvaerr=bvaerr+valerr
              endif
           endif
              
           !
           SUMRAD=RMT(JATOM)+RMT(N3)                                        
           if(dists(n1).lt.dfac*dists(1)) then
              if(SHOWALL.or.(M.EQ.1)) write(66,3) N3,NAME(N3),(PNN(L,N2),L=1,3),DISTS(N1),DISTS(N1)*0.529177  
              IF(ovlap(n3).and.SUMRAD.GE.DISTS(N1)) THEN
                 ovlap(n3)=.false.                        
                 write(66,4) JATOM,RMT(JATOM),N3,RMT(N3), SUMRAD, DISTS(N1)         
                 WRITE(*,4) JATOM,RMT(JATOM),N3,RMT(N3), SUMRAD, DISTS(N1)         
              end if
           end if
           !
           if((dists(n1)-olddis).gt.dlimit) then
              !.....new shell
              ishell(index)=ishell(index)+1
              iz(index,ishell(index))=1
              olddis=dists(n1)  
              if(ishell(index).ge.mshell) EXIT !goto 187
              icnt(index,ishell(index),iz(index,ishell(index)))=1
              shellz(index,ishell(index),iz(index,ishell(index)))=zz(n3)
              shdist(index,ishell(index))=olddis
           else
              !.....old shell
              do i=1,iz(index,ishell(index))
                 if(zz(n3).eq.shellz(index,ishell(index),i)) then
                    icnt(index,ishell(index),i)=icnt(index,ishell(index),i)+1
                    EXIT !goto 186
                 endif
              enddo
              iz(index,ishell(index))=iz(index,ishell(index))+1
              shellz(index,ishell(index),iz(index,ishell(index)))=zz(n3)
              icnt(index,ishell(index),iz(index,ishell(index)))=1          
           endif
           !186         continue
           !
        ENDDO !185 CONTINUE
        !187      CONTINUE
        if( (SHOWALL.or.(M.EQ.1)).and.(bva .gt. 0.001D0) )then
           write(66,'(A,i3,A,i2,A,A,A,2F8.2)') 'Atom ',JATOM,' equiv ',M,' ',NAME(JATOM),' Bond-Valence Sum ',bva,bvaerr
        endif
        !
        !.....limit shells to some maximum
        write(67,*) '---------------- ATOM',index
        do i=1,ishell(index)
           if(shdist(index,i).gt.dstmax-dstmax/3.d0) then
              ishell(index)=i
              !goto 190
              EXIT
           endif
           write(67,*) ' SHELL:',i,shdist(index,i)
           do j=1,iz(index,i)
              write(67,*) icnt(index,i,j),shellz(index,i,j)
           enddo
        enddo
        !
     ENDDO !  190 CONTINUE
  ENDDO !200 CONTINUE                                                          
  !                      
  write(66,552)
  INDEX=0                                                           
  DO  JATOM=1,NAT                
     DO  M=1,MULT(JATOM)                                            
        INDEX=INDEX+1
        zzo(index)=zz(jatom)                                                  
        write(66,553)index,((shdist(index,i),icnt(index,i,j),shellz(index,i,j),j=1,iz(index,i)),i=1,4)
     enddo
  enddo
  write(66,*)
  !
  inat=0
  do i=1,index
     ityp(i)=i
     imult(i)=0
     used(i)=.false.  
  enddo
  write(66,554)
  ityp1=0    
  ij=0          
  i=0                           
  do i0=1,nat            ! 500
     do i00=1,mult(i0)   ! 500
        i=i+1
        if(used(i)) CYCLE !goto 500
        write(66,*)
        do j=i,index     ! 501
           if(zz(i0).ne.zzo(j)) CYCLE !goto 501
           !     compare atom with index i with shells of other atoms (Fuhr)
           if (ishell(i).ne.ishell(j)) then
              write(66,559) i,j,ishell(i),ishell(j)
              CYCLE !goto 501
           endif
           !     compare atom with index i with all other atoms
           do i1=1,ishell(i)-1   ! 510
              if(abs(shdist(i,i1)-shdist(j,i1)).gt.dlimit.and.shdist(i,i1).lt.2.d000*max(a(1),a(2),a(3))) then
                 write(66,550) i,j,i1,shdist(i,i1),shdist(j,i1)
                 goto 501
              endif
              do i2=1,iz(i,i1)   ! 511
                 do j2=1,iz(j,i1)
                    if((icnt(i,i1,i2).eq.icnt(j,i1,j2)).and.(shellz(i,i1,i2).eq.shellz(j,i1,j2))) goto 511
                    if(shdist(i,i1).gt.2.d0000*max(a(1),a(2),a(3))) goto 511
                 enddo
                 write(66,551) i,j,i1,icnt(i,i1,i2),icnt(j,i1,j2-1),shellz(i,i1,i2),shellz(j,i1,j2-1),shdist(i,i1)
                 goto 501
511              CONTINUE
              enddo !511              continue
           enddo       !510           continue
           write(66,555) i,j
           if(i.eq.j) then
              ityp1=ityp1+1
              namen(ityp1)=name(i0)
              JRJN(ityp1)=jrj(i0)
              R0N(ityp1)=r0(i0)
              RMTTN(ityp1)=rmtt(i0)
              zzn(ityp1)=zzorig(i0)
           endif
           ityp(j)=ityp1
           if(inat.lt.ityp(j)) inat=ityp(j)
           imult(ityp1)=imult(ityp1)+1
           used(j)=.true.
           ij=ij+1
           posn(1,ij)=pos(1,j)
           posn(2,ij)=pos(2,j)
           posn(3,ij)=pos(3,j)
501        continue
        ENDDO
     ENDDO ! 500
  ENDDO ! 500

  write(66,*)      
  error=.false.
  INDEX=0                                                           
  DO  JATOM=1,NAT                
     write(66,556) jatom,mult(jatom),imult(jatom)  
     if(mult(jatom).ne.imult(jatom)) then
        error=.true.
        write(66,*)'WARNING: MULT not equal. The new multiplicity is',' different from the old one'  
        write(6,*)'WARNING: Mult not equal. PLEASE CHECK outputnn-file'
     end if
     DO  M=1,MULT(JATOM)                                            
        INDEX=INDEX+1                                                     
        if(jatom.ne.ityp(index)) then
           error=.true.
           write(66,557) index,jatom,ityp(index)
           write(66,*) 'WARNING: ITYP not equal. The new type is',' different from the old one'  
           write(6,*)'WARNING: ityp not equal. PLEASE CHECK outputnn-file'
        endif
     enddo
  enddo
  !

  open(68, FILE=f_68, STATUS='unknown', FORM='formatted')
  WRITE(68,*) 'BR1'
  DO JR=1,3
     WRITE(68, '(3F15.10)') BR1(1,JR), BR1(2,JR), BR1(3,JR) !--- Writting BR1 ----!
  ENDDO
  !WRITE(68,*) 'BR2'
  !DO JR=1,3
  !   WRITE(68, '(3F15.10)') BR2(1,JR), BR2(2,JR), BR2(3,JR) !--- Writting BR1 ----!
  !ENDDO
  
  close(20)
  close(66)
  close(67)
  close(68)
  if(.not.error) return !stop 'NN ENDS'
  
  write(66,*)
  write(66,*) 'NEW LIST of EQUIVALENT POSITIONS written to',' case.struct_nn'   
  !.....write new struct file
  !

  open(21, FILE=f_21, STATUS='unknown', FORM='formatted')
  write(21,1510) TITLE
  write(21,1511) LATTIC,inat,title                                  
  write(21,'(6F10.6)') A(1),A(2),A(3),alpha,beta,gamma
  index=0
  do jatom=1,inat
     INDEX=INDEX+1                                                     
     write(21,1012) -JATOM,( POSN(J,INDEX),J=1,3 ),iMULT(JATOM),8
     write(66,*)
     write(66,1011) -JATOM,( POSN(J,INDEX),J=1,3 ),NAMEN(JATOM),zzn(jatom) 
     DO M=1,iMULT(JATOM)-1                                     
        INDEX=INDEX+1                                            
        write(21,1011) -JATOM,( POSN(J,INDEX),J=1,3)         
        write(66,1011) -JATOM,( POSN(J,INDEX),J=1,3)         
     ENDDO
     if(jatom.lt.10) then
        write(21,1149) NAMEN(JATOM)(1:2),jatom,JRJN(JATOM),R0N(JATOM),RMTTN(jatom),zzn(jatom)       
     else if(jatom.lt.100) then
        write(21,1148) NAMEN(JATOM)(1:2),jatom,JRJN(JATOM),R0N(JATOM),RMTTN(jatom),zzn(jatom)       
     else
        write(21,1147) NAMEN(JATOM)(1:2),jatom,JRJN(JATOM),R0N(JATOM),RMTTN(jatom),zzn(jatom)       
     endif
     write(21,1512) 1.0,0.0,0.0
     write(21,1512) 0.0,1.0,0.0
     write(21,1512) 0.0,0.0,1.0
  enddo
  write(21,*) 0
  !                 
  WRITE(6,*) ' '
  WRITE(6,*) " NN created a new ",trim(fname1),"_nn file"   
                                                     
  close(21)

!3 FORMAT(' ATOM:',I3,2X,A10,'AT',3F8.4,' IS',F9.5,' A.U.',f10.5,' ANG')
3 FORMAT(' ATOM:',I3,2X,A10,'AT',3F14.10,' IS',F9.5,' A.U.',f10.5,' ANG')

550 format(' atom:',i4,' and ATOM:',i4,' differ in shell',i3,2f10.5)
551 format(' atom:',i4,' and ATOM:',i4,' differ in shell',i3,2i4,2f7.1,f10.5)
552 format(//,'SHELL STRUCTURE for all ATOMS:',/,'ATOM  | DISTANCE   #of NEIGHBORS   Z |')
553 format(i3,1x,8('|',f6.3,i3,f6.1))
554 format(/,'LISTING of INDEX of all EQUIVALENT ATOMS:')
555 format(' ATOM:',i4,' and ATOM:',i4,' are equivalent')
556 format(/,' ATOM KIND:',i4,'  OLD and NEW MULTIPLICITY:  ',2i4)
557 format(5x,'ATOM INDEX:',i4,'  OLD and NEW ATOM KIND:',2i4)
559 format(' atom:',i4,' and ATOM:',i4,' differ in number of shells',i4,'ne',i4)

1011 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8,5x,a10,f5.1)         
1012 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8,/,15X,I2,17X,I2)     
1020 FORMAT(///,3X,'ERROR IN EBND  : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
1049 FORMAT(A10,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f5.1)      
1149 FORMAT(A2,i1,7x,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f5.1)  
1148 FORMAT(A2,i2,6x,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f5.1)  
1147 FORMAT(A2,i3,5x,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f5.1)  
1050 FORMAT(A10,5X,I5,5X,F10.8,5X,F10.4,5x,f5.1)                           
1051 FORMAT(20X,3F10.7)                                                
1510 FORMAT(A79)                                                       
1511 FORMAT(A4,23X,I3,/,13X,A4)                                        
1512 FORMAT(20x,3f10.7)


  
END SUBROUTINE w2knn



SUBROUTINE REDUCE_ALLOC_R8_2D(N1OLD,N2OLD,N1NEW,N2NEW)
  !..REDUCES UP TO 4 DIMENSIONAL REAL*8 ARRAYS FROM OLD TO NEW ALLOCATIONS
  !
  use variable_fields,  a => pos 
  real*8,  allocatable ::  HELP(:,:)
  ALLOCATE ( HELP(N1NEW,n2new) )
  DO I=1,N1NEW
     DO j=1,N2NEW
        HELP(I,j)=A(I,j)
     ENDDO
  ENDDO
  DEALLOCATE ( A )
  ALLOCATE ( A(N1NEW,n2new) )
  DO I=1,N1NEW
     DO j=1,N2NEW
        A(I,j)=HELP(I,j)
     ENDDO
  ENDDO
  DEALLOCATE ( HELP )      
END SUBROUTINE REDUCE_ALLOC_R8_2D


SUBROUTINE DIRLAT (BR1,BR2,ortho,lattic,nat,alpha,beta,gamma,A)
  !                                                                       
  !     LATGEN GENERATES THE BRAVAIS MATRIX BR2(3,3), WHICH TRANSFORMS    
  !     A RECIPROCAL LATTICE VECTOR OF A SPECIAL COORDINATE SYSTEM (IN    
  !     UNITS OF 2 PI/A(I)) INTO CARTESIAN SYSTEM                         
  !     Convention:  R_i = (i,*)
  !                                                                       
  !use struk, only: br2,lattic
  IMPLICIT REAL*8 (A-H,O-Z)
  REAL*8, intent(out)    :: BR1(3,3), BR2(3,3)
  LOGICAL, intent(out)   :: ORTHO
  CHARACTER*4, intent(in):: lattic
  REAL*8, intent(in)     :: alpha,beta,gamma,A(3)
  REAL*8 :: PIA(3), SQRT3!, _ALPHA_(3)

  WRITE(6,*) 'A=', A
  pi = ACOS(-1.0D0)
  SQRT3=SQRT(3.0D+0)
  !_ALPHA_(1) = alpha*pi/180.d0
  !_ALPHA_(2) = beta*pi/180.d0
  !_ALPHA_(3) = gamma*pi/180.d0
  PIA(1)=2.D0*PI/A(1)
  PIA(2)=2.D0*PI/A(2)
  PIA(3)=2.D0*PI/A(3)
  
  gamma1=gamma*pi/180.d0
  beta1=beta*pi/180.d0
  alpha1=alpha*pi/180.d0
  cosg1=(cos(gamma1)-cos(alpha1)*cos(beta1))/sin(alpha1)/sin(beta1)
  gamma0=acos(cosg1)


  WRITE(6,*) 'PIA=', PIA, 'LATTIC=', LATTIC(1:1)

  IF(LATTIC(1:1).EQ.'H') GOTO 10                                    
  IF(LATTIC(1:1).EQ.'S') GOTO 20                                    
  IF(LATTIC(1:1).EQ.'P') GOTO 20                                    
  IF(LATTIC(1:1).EQ.'B') GOTO 30                                    
  IF(LATTIC(1:1).EQ.'F') GOTO 40                                    
  IF(LATTIC(1:1).EQ.'C') GOTO 50                                    
  IF(LATTIC(1:1).EQ.'R') GOTO 80                                    
  STOP 'LATTIC WRONG'                                               
  !                                                                       
  !.....HEXAGONAL CASE                                                    
10 CONTINUE
  BR2(1,1)=SQRT(3.d0)/2.d0                                              
  BR2(1,2)=-.5d0                                                      
  BR2(1,3)=0.0d0                                                      
  BR2(2,1)=0.0d0                                                     
  BR2(2,2)=1.0d0                                                      
  BR2(2,3)=0.0d0                                                      
  BR2(3,1)=0.0d0                                                      
  BR2(3,2)=0.0d0                                                      
  BR2(3,3)=1.d0                                                       
  ORTHO=.FALSE.
  !
  BR1(1,1)=2.D0/SQRT3*PIA(1)                                        
  BR1(1,2)=1.D0/SQRT3*PIA(1)                                        
  BR1(1,3)=0.0D0                                                    
  BR1(2,1)=0.0D0                                                    
  BR1(2,2)=PIA(2)                                                   
  BR1(2,3)=0.0D0                                                    
  BR1(3,1)=0.0D0                                                    
  BR1(3,2)=0.0D0                                                    
  BR1(3,3)=PIA(3)                                                   

  WRITE(6,*) 'BR1=', BR1

  GOTO 100                                                          
  !                                                                       
  !.....PRIMITIVE LATTICE CASE 
20 continue
  !
  BR2(1,1)=1.0d0*sin(gamma0)*sin(beta1)                
  BR2(1,2)=1.0d0*cos(gamma0)*sin(beta1)                 
  BR2(1,3)=1.0d0*cos(beta1)                                   
  BR2(2,1)=0.0d0                                                      
  BR2(2,2)=1.0d0*sin(alpha1)
  BR2(2,3)=1.0d0*cos(alpha1)
  BR2(3,1)=0.0d0                                                      
  BR2(3,2)=0.0d0                                                      
  BR2(3,3)=1.0d0                                                      
  ORTHO=.TRUE. 
  if(gamma.ne.90.d0) ortho=.false.                                      
  if(beta.ne.90.d0) ortho=.false.                                      
  if(alpha.ne.90.d0) ortho=.false.                                      
  !        write(*,*) alpha0,beta0,gamma0,ortho,br2
  !
  SINBC=SIN(alpha1)
  COSAB=COS(gamma1)
  COSAC=COS(beta1)
  COSBC=COS(alpha1)
  WURZEL=SQRT(SINBC**2-COSAC**2-COSAB**2+2*COSBC*COSAC*COSAB)
  BR1(1,1)= SINBC/WURZEL*PIA(1)
  BR1(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
  BR1(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
  BR1(2,1)= 0.0
  BR1(2,2)= PIA(2)/SINBC
  BR1(2,3)= -PIA(3)*COSBC/SINBC
  BR1(3,1)= 0.0
  BR1(3,2)= 0.0
  BR1(3,3)= PIA(3)
  GOTO 100     
  !                                                                       
  !.....BC CASE (DIRECT LATTICE)                                          
30 CONTINUE                                                          
  BR2(1,1)=-0.5d0                                                     
  BR2(1,2)=0.5d0                                                      
  BR2(1,3)=0.5d0                                                      
  BR2(2,1)=0.5d0                                                     
  BR2(2,2)=-0.5d0                                                     
  BR2(2,3)=0.5d0                                                      
  BR2(3,1)=0.5d0                                                      
  BR2(3,2)=0.5d0                                                      
  BR2(3,3)=-0.5d0                                                     
  ORTHO=.TRUE.
  !
  BR1(1,1)=PIA(1)                                                   
  BR1(1,2)=0.0D0                                                    
  BR1(1,3)=0.0D0                                                    
  BR1(2,1)=0.0D0                                                    
  BR1(2,2)=PIA(2)                                                   
  BR1(2,3)=0.0D0                                                    
  BR1(3,1)=0.0D0                                                    
  BR1(3,2)=0.0D0                                                    
  BR1(3,3)=PIA(3)                                                   
  GOTO 100                                                          
  !                                                                       
  !.....FC CASE (DIRECT LATTICE)                                          
40 CONTINUE                                                          
  BR2(1,1)=0.0d0                                                      
  BR2(1,2)=0.5d0                                                      
  BR2(1,3)=0.5d0                                                      
  BR2(2,1)=0.5d0                                                      
  BR2(2,2)=0.0d0                                                      
  BR2(2,3)=0.5d0                                                      
  BR2(3,1)=0.5d0                                                      
  BR2(3,2)=0.5d0                                                      
  BR2(3,3)=0.0d0                                                      
  ORTHO=.TRUE.
  !
  BR1(1,1)=PIA(1)                                                   
  BR1(1,2)=0.0D0                                                    
  BR1(1,3)=0.0D0                                                    
  BR1(2,1)=0.0D0                                                    
  BR1(2,2)=PIA(2)                                                   
  BR1(2,3)=0.0D0                                                    
  BR1(3,2)=0.0D0                                                    
  BR1(3,1)=0.0D0                                                    
  BR1(3,3)=PIA(3)                                                   
  GOTO 100                                                          
  !                                                                       
  !.....CXY  CASE (DIRECT LATTICE)                                          
50 CONTINUE                                                          
  IF(LATTIC(2:3).EQ.'XZ') GOTO 60                                    
  IF(LATTIC(2:3).EQ.'YZ') GOTO 70                                    
  BR2(1,1)=0.5d0                                                      
  BR2(1,2)=-0.5d0                                                     
  BR2(1,3)=0.0d0                                                      
  BR2(2,1)=0.5d0                                                      
  BR2(2,2)=0.5d0                                                      
  BR2(2,3)=0.0d0                                                      
  BR2(3,1)=0.0d0                                                      
  BR2(3,2)=0.0d0                                                      
  BR2(3,3)=1.0d0                                                      
  ORTHO=.TRUE.
  !
  BR1(1,1)=PIA(1)                                                   
  BR1(1,2)=0.0D0                                                    
  BR1(1,3)=0.0D0                                                    
  BR1(2,1)=0.0D0                                                    
  BR1(2,2)=PIA(2)                                                   
  BR1(2,3)=0.0D0                                                    
  BR1(3,1)=0.0D0                                                    
  BR1(3,2)=0.0D0                                                    
  BR1(3,3)=PIA(3)                                                   
  GOTO 100                                                          
  !                                                                       
  !.....CXZ  CASE (DIRECT LATTICE)                                          
60 CONTINUE 
  !.....CXZ ORTHOROMBIC CASE
  if(gamma.eq.90.d0) then
     BR2(1,1)=0.5d0                                                      
     BR2(1,2)=0.0d0                                                     
     BR2(1,3)=-0.5d0                                                      
     BR2(2,1)=0.0d0                                                      
     BR2(2,2)=1.0d0                                                      
     BR2(2,3)=0.0d0                                                      
     BR2(3,1)=0.5d0                                                     
     BR2(3,2)=0.0d0                                                      
     BR2(3,3)=0.5d0                                                      
     ORTHO=.TRUE.
     !
     BR1(1,1)=PIA(1)                                                   
     BR1(1,2)=0.0D0                                                    
     BR1(1,3)=0.0D0                                                    
     BR1(2,1)=0.0D0                                                    
     BR1(2,2)=PIA(2)                                                   
     BR1(2,3)=0.0D0                                                    
     BR1(3,1)=0.0D0                                                    
     BR1(3,2)=0.0D0                                                    
     BR1(3,3)=PIA(3)                                                   
     GOTO 100                                                          
  ELSE
     !.....CXZ MONOCLINIC CASE
     write(*,*) 'gamma not equal 90'
     SINAB=SIN(gamma1)
     COSAB=COS(gamma1)
     !
     BR2(1,1)=0.5d0*sinab                                                
     BR2(1,2)=0.5d0*cosab                                               
     BR2(1,3)=-0.5d0                                                      
     BR2(2,1)=0.0d0                                                      
     BR2(2,2)=1.0d0                                                      
     BR2(2,3)=0.0d0                                                      
     BR2(3,1)=0.5d0*sinab                                               
     BR2(3,2)=0.5d0*cosab                                                
     BR2(3,3)=0.5d0                                                      
     ORTHO=.FALSE.
     !
     BR1(1,1)= PIA(1)/SINAB 
     BR1(1,2)= -PIA(2)*COSAB/SINAB
     BR1(1,3)= 0.0                                                   
     BR1(2,1)= 0.0                                                      
     BR1(2,2)= PIA(2)                                                     
     BR1(2,3)= 0.0                                                     
     BR1(3,1)= 0.0                                                     
     BR1(3,2)= 0.0                                                     
     BR1(3,3)= PIA(3)                                                     
     GOTO 100
  ENDIF

  !                                                                       
  !.....CYZ  CASE (DIRECT LATTICE)                                          
70 CONTINUE                                                          
  BR2(1,1)=1.0d0                                                      
  BR2(1,2)=0.0d0                                                     
  BR2(1,3)=0.0d0                                                      
  BR2(2,1)=0.0d0                                                      
  BR2(2,2)=0.5d0                                                      
  BR2(2,3)=0.5d0                                                      
  BR2(3,1)=0.0d0                                                      
  BR2(3,2)=-0.5d0                                                     
  BR2(3,3)=0.5d0                                                      
  ORTHO=.TRUE.
  !!
  BR1(1,1)=PIA(1)                                                   
  BR1(1,2)=0.0D0                                                    
  BR1(1,3)=0.0D0                                                    
  BR1(2,1)=0.0D0                                                    
  BR1(2,2)=PIA(2)                                                   
  BR1(2,3)=0.0D0                                                    
  BR1(3,1)=0.0D0                                                    
  BR1(3,2)=0.0D0                                                    
  BR1(3,3)=PIA(3)                                                   
  GOTO 100
  
  !.....RHOMBOHEDRAL CASE
80 CONTINUE
  BR2(1,1)=1/2.d0/sqrt(3.d0)
  BR2(1,2)=-1/2.d0                                                     
  BR2(1,3)=1/3.d0                                                      
  BR2(2,1)=1/2.d0/SQRT(3.d0)                                          
  BR2(2,2)=1*0.5d0                                                
  BR2(2,3)=1/3.d0                                                      
  BR2(3,1)=-1/SQRT(3.d0)                                         
  BR2(3,2)=0.d0                                                
  BR2(3,3)=1/3.d0                                                      
  ORTHO=.FALSE.
  !
  BR1(1,1)=1.D0/SQRT(3.D0)*PIA(1)                                          
  BR1(1,2)=1.D0/SQRT(3.D0)*PIA(1)                                          
  BR1(1,3)=-2.d0/sqrt(3.d0)*PIA(1)                                         
  BR1(2,1)=-1.0d0*PIA(2)                                                  
  BR1(2,2)=1.0d0*PIA(2)                                                    
  BR1(2,3)=0.0d0*PIA(2)                                                    
  BR1(3,1)=1.0d0*PIA(3)                                                    
  BR1(3,2)=1.0d0*PIA(3)                                                    
  BR1(3,3)=1.0d0*PIA(3)                                                    
  GOTO 100                                                         
  !                                                                       
100 CONTINUE                                                          
  write(66,*) 'Bravais Matrix:'
  write(66,999) br2
999 format(3f15.5)
  !                                                                       
  RETURN                                                            
END SUBROUTINE DIRLAT


SUBROUTINE ORD2(A,NR,IMAX)                                        
  !     ORDERS ELEMENTS IN ARRAY A INCREASING IN SIZE                     
  !       REORDERS CORRESPONDING INDICES (NR)                             
  IMPLICIT REAL*8 (A-H,O-Z)
  LOGICAL CONT                                                      
  DIMENSION A(*),NR(*)                                              
  DO I=1,IMAX  ! 50
     NR(I)=I
  ENDDO
  if(imax.eq.1) return
  
  CONT=.TRUE.
  do while(CONT)
     CONT=.FALSE.
     DO I=2,IMAX
        IF(A(I).LT.A(I-1)) THEN
           !       INTERCHANGE I AND (I-1) ELEMENT IN ARRAYS A AND NR              
           HLP=A(I)                                                        
           A(I)=A(I-1)                                                     
           A(I-1)=HLP                                                      
           NHLP=NR(I)                                                      
           NR(I)=NR(I-1)                                                   
           NR(I-1)=NHLP                                                    
           CONT=.TRUE.                                                     
        ENDIF
     ENDDO
  enddo
  RETURN                                                
END SUBROUTINE

