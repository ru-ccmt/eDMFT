subroutine read_inso(emm,xms,fl,jspin,kpot,ipr,irlotot)
  USE param, ONLY: labc, labc2, nato, lomax, lmax, nrad
  USE loabcr, ONLY:
  USE lolog, ONLY: lso
  USE rlolog, ONLY: loorext
  USE rpars, ONLY: nlr, extei, extde
  USE structure, ONLY: mult
  USE vns, ONLY:
  USE couples, ONLY: allocate_couplo
  USE MPI, ONLY: Qprint
  IMPLICIT NONE
  ! ipr  -- parameter for debug output;
  ! kpot -- determines how the potential for spin-polarized calculation is determined.
  ! fl   -- also eigenvectors are needed
  INTEGER, intent(out) :: ipr, kpot, irlotot   
  LOGICAL, intent(out) :: fl
  INTEGER, intent(in)  :: jspin
  REAL*8, intent(out)  :: emm(2), xms(3)
  ! locals
  character*5 ::  vect
  REAL*8      :: edum1, edum2
  INTEGER     :: iatoff(nato), jatom, nrelcase, noff, i, j, l
  logical     :: italk
  !real*8, parameter :: pi = 3.141592653589793238462643d0

  !   Read parameters from unit 5,i.e. case.inso
  read(5,'(A5)') vect
  read(5,*)labc,ipr,kpot
  labc2=(labc+1)**2

  if (vect.eq.'WFFIL') fl=.true.

  if (jspin.eq.2) then
     if (kpot.eq.1) then
        if (Qprint) write(8,*)' Averaged potential used when calculating dV/dr'
        if (Qprint) write(6,*)' Averaged potential used when calculating dV/dr'
     else
        if (Qprint) write(8,*)' Potential not averaged when calculating dV/dr'
        if (Qprint) write(6,*)' Potential not averaged when calculating dV/dr'
     endif
  endif
  if (fl) then
     if (Qprint) write(6,*)'LMAX=',labc,' S-O eigenvectors and eigenvalues'
  else
     if (Qprint) write(6,*)'LMAX=',labc,' S-O eigenvalues only'
  endif

  !  Read energy ranges for which H_so is calculated and
  !  orientation for polarization axis for s-p-systems.
  read(5,*) emm(1),emm(2)         ! emin, emax
  read(5,*) xms(1),xms(2),xms(3)  ! quantization direction

  if (Qprint) write(6,'(A,f19.4,A,f19.4)') ' Emin=',emm(1),' Emax=',emm(2)
  if (emm(2) .lt. emm(1)) stop 'EMAX < EMIN'

  !  read in parameters for basis extension by relativistic basis function.
  nlr(1:nato)=0
  read(5,*,end=644)nrelcase
  do i=1,nrelcase 
     read(5,*) jatom,edum1,edum2
     extei(jatom,1)=edum1
     extde(jatom,1)=edum2
     nlr(jatom)=1 
  enddo
644 continue

  ! how many atoms have SO switched off ?
  italk=.true.
  read(5,*,end=343) noff,(iatoff(i),i=1,noff)
  goto 345
343 continue
  italk=.false.
345 continue
  if (italk) then
     do i=1,noff
        lso(iatoff(i)) = .false.
        if (Qprint) write(6,*)'SOC on atom ',iatoff(i),'SWITCHED OFF !!!'
     enddo
  endif

  ! sort entries for rlo extension by angular momentum.
  do jatom = 1, nato
     do l = 0, lomax
        loorext(l,jatom) = .false.
        do j = 1, nlr(jatom)
           loorext(1,jatom) = .true.    ! It seems only p-orbitals have these...
           irlotot = (2*1+1)*mult(nato) ! these are p_1/2 states....
        enddo
     enddo
  enddo

  return
!549 FORMAT(A5)
END subroutine read_inso


