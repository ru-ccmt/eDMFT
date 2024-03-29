! @Copyright 2007 Kristjan Haule
! 

SUBROUTINE ReadFirstPart_of_in1(PRNTWF, SPRSWF, WFTAPE, efermi, NT, RKM, LNSMAX, NSLMAX, LMAX)
  !        read output-format parameters
  IMPLICIT NONE
  LOGICAL, intent(out) :: PRNTWF, SPRSWF, WFTAPE
  REAL*8, intent(out)  :: efermi, RKM
  INTEGER, intent(out) :: LNSMAX, NT
  INTEGER, intent(in)  :: NSLMAX, LMAX
  ! locals
  character   :: line*80
  CHARACTER*5 :: MODUS
  INTEGER     :: NT1, i
  READ(5,'(a)') line
  read(line,'(a5)') MODUS    ! beginning of case.in1
  IF (MODUS .EQ. 'WFPRI') THEN
     PRNTWF = .TRUE.
     SPRSWF = .FALSE.
     WFTAPE = .TRUE.
  ELSEIF (MODUS .EQ. 'WFFIL') THEN  ! default
     PRNTWF = .FALSE.
     SPRSWF = .FALSE.
     WFTAPE = .TRUE.
  ELSEIF (MODUS .EQ. 'ENFIL') THEN
     PRNTWF = .FALSE.
     SPRSWF = .TRUE.
     WFTAPE = .FALSE.
  ELSE
     PRNTWF = .FALSE.
     SPRSWF = .TRUE.
     WFTAPE = .FALSE.
  ENDIF
  ! From the first line of case.in1 extracts the fermi energy
  efermi=0.5d0
  do i=6,70
     if(line(i:i+1).eq.'EF') then
        read(line(i+3:80),*,err=1234) efermi
     endif
  enddo
1234 continue
  !
  !        read in R-MT times K-MAX (used to check convergence)
  !
  READ(5,*) RKM, NT1, LNSMAX   ! RKM==R-MT*K-MAX; NT1==MAX L IN WF; LNSMAX==V-NMT
  IF (LNSMAX .EQ. 0) LNSMAX = 2
  IF (LNSMAX .GE. NSLMAX) LNSMAX = NSLMAX - 1  ! can not be larger than nslmax=5
  IF (NT1 .GT. 0)    NT = NT1 + 1
  IF (NT  .GT. LMAX) NT = LMAX                 ! can not be greated that lmax=13
  ! NT is lmax+1 but updated to be compatible with case.in1 specified maximum l.
END SUBROUTINE ReadFirstPart_of_in1
