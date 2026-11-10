SUBROUTINE readstruct
  ! Reads struct file and stores in module struct
  ! 23/2-01 GM
  USE struct
  USE symetr
  USE reallocate
  use parallel, only: abort_parallel  
  implicit real*8 (a-h,o-z)
  REAL*8              :: test,ninety,postest(3)
  CHARACTER*80    WSPACE
  LOGICAL there, response_scale_exists
  REAL*8 response_scale
  TEST=1.D-5
  NINETY=90.0D0
  
  READ(20,1000) title                                               
  READ(20,1010) lattic,nat,cform,irel    
  ALLOCATE (rmt(nat),iatnr(nat),mult(nat),isplit(nat))
  ALLOCATE (rotloc(3,3,nat))
  ALLOCATE (r0(nat),dx(nat),jri(nat))
  ALLOCATE (aname(nat))
  ALLOCATE (zz(nat))
  ALLOCATE (pos(3,48*nat) )
  !     READ IN LATTICE CONSTANTS                                         
  READ(20,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
  ! Diagnostic only: avoid 6-decimal lattice rounding in strain finite differences.
  ! RMT, fractional sites, axes and native potential algorithms are unchanged.
  inquire(file='lapw0.strain_scale',exist=response_scale_exists)
  if(response_scale_exists)then
     open(994,file='lapw0.strain_scale',status='old',action='read')
     read(994,*)response_scale
     close(994)
     if(.not.(response_scale.gt.0.d0))error stop 'invalid diagnostic strain scale'
     aa=aa*response_scale
     bb=bb*response_scale
     cc=cc*response_scale
     write(6,'(A,ES24.16)')'DIAGNOSTIC uniform lattice scale: ',response_scale
  endif
  a(1)=aa
  a(2)=bb
  a(3)=cc
  IF(ABS(alpha(1)).LT.test) alpha(1)=ninety
  IF(ABS(alpha(2)).LT.test) alpha(2)=ninety
  IF(ABS(alpha(3)).LT.test) alpha(3)=ninety
!  IF(irel.EQ.'RELA') rel=.TRUE.                                     
!  IF(irel.EQ.'NREL') rel=.FALSE.                                    
!
!        read crystal-structure (atompositions, symmetry-parameters,
!                                muffin-tin radius, ...)
!        'INDEX' counts all atoms in the unit cell,
!        'JATOM' counts only the notequivalent atoms
!
  index=0                                                           
  DO jatom = 1,nat                                               
     index=index+1
     !
     READ(20,1030,iostat=ios) iatnr(jatom),( pos(j,index),j=1,3 ), &
          mult(jatom),isplit(jatom) 
     IF(ios /= 0) THEN
        WRITE(6,*) iatnr(jatom),( pos(j,index),j=1,3 ), &
          mult(jatom),isplit(jatom) 
        WRITE(6,*) 'ERROR IN STRUCT FILE READ'
        call abort_parallel
        STOP
     ENDIF
     IF (mult(jatom) .EQ. 0) THEN
        !
        !           illegal number of equivalent atoms
        !
        WRITE (6,6000) jatom, index, mult(jatom)
        GOTO 930
     ENDIF
     DO m=1,mult(jatom)-1                                     
        index=index+1                                            
        !         if(isplit(jatom).eq.999.and.lxdos.ne.3) goto 956
        !
        READ(20,1031) iatnr(jatom),( pos(j,index),j=1,3)         
     ENDDO
     READ(20,1050) aname(jatom),jri(jatom),r0(jatom),rmt(jatom), &
          zz(jatom)
     dx(jatom)=LOG(rmt(jatom)/r0(jatom)) / (jri(jatom)-1)           
     rmt(jatom)=r0(jatom)*EXP( dx(jatom)*(jri(jatom)-1) )           
     READ(20,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)                
  ENDDO  
  READ(20,1151) iord
  if(iord.eq.0) then
    CALL OUTERR('LAPW0','IORD .EQ. 0')
    GOTO 999
  endif
  ALLOCATE(iz(3,3,iord),tau(3,iord),inum(iord))
  nsym=iord
  DO j=1,iord
     READ(20,1101) ( (iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3 ),inum(j)
  ENDDO
  ndif=index
  ALLOCATE (rotij(3,3,index),tauij(3,index))
  CALL doreallocate(pos, 3, index)
! Correct Rotation matrices
  call SymmRot(rotloc,NAT)
!
!     Are high-precision atoms available ?
        read(20,'(A)',err=1999,end=1999)WSPACE
        there=.false.
        if(WSPACE(1:17).eq.'Precise positions') there=.true.
        if(there)then
                jj=0
                do jatom=1,nat
                        do j=1,mult(Jatom)
                                jj=jj+1
!                                read(20,*)pos(1:3,jj)
                                read(20,*)postest(1:3)
                           do i=1,3
           if(abs(modulo(postest(i),1.d0)-modulo(pos(i,jj),1.d0)).gt.1.d-8) then
              print*, 'Precise positions overwritten',postest(i),pos(i,jj) 
              goto 1999 
           endif
                           enddo
                           pos(1:3,jj)=postest(1:3)
                        enddo
                enddo
        endif
1999    continue

  RETURN

930 INFO = 3
  CALL OUTERR('LAPW0','MULT .EQ. 0')
  GOTO 999
999 call abort_parallel
  STOP 'LAPW0 - Error'
  
1000 FORMAT(A80)                                                       
1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
1020 FORMAT(6F10.6)                                          
1030 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8,/,15X,I2,17X,I2)          
1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
1051 FORMAT(20X,3F10.8)                                             
1151 FORMAT(I4)
1101 FORMAT(3(3I2,F11.8/),I8)
6000 FORMAT(///,3X,'ERROR IN LAPW0 : MULT(JATOM)=0 ...', &
          /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  
END SUBROUTINE readstruct
