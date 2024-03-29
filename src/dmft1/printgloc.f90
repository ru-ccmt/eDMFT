! @Copyright 2007 Kristjan Haule
! 

module PrintG
contains
  SUBROUTINE PrintGloc(fh_dos, fh_gc, fh_dt, Glc, gloc, gtot, Deltac, omega, csize, csizes, nl, ll, legend, iatom, ncix, nomega, natom, norbitals, maxsize, Ry2eV)
    !-- In the old version we used:  fh_dos = 100, fh_gc = 120, fh_dt = 140
    !---Currently we use          :  fh_dos = 500, fh_gc = 180, fh_dl = 280
    USE com_mpi, ONLY: myrank, master 
    USE param, ONLY: cmp_partial_dos
    IMPLICIT NONE
    INTEGER, intent(in)    :: fh_dos, fh_gc, fh_dt
    COMPLEX*16, intent(in) :: Glc(maxsize,ncix,nomega), gtot(nomega), Deltac(maxsize,ncix,nomega)
    COMPLEX*16, allocatable, intent(in) :: gloc(:,:) ! gloc(norbitals,nomega)
    REAL*8, intent(in)     :: omega(nomega)
    INTEGER, intent(in)    :: csize(ncix), csizes(ncix), iatom(natom)
    INTEGER, intent(in)    :: nl(natom), ll(natom,4)
    CHARACTER*30, intent(in):: legend(maxsize,ncix)
    INTEGER, intent(in) :: nomega, norbitals, natom, ncix, maxsize
    REAL*8, intent(in)  :: Ry2eV
    ! local
    INTEGER :: L, iom, lcase, i, icix, itape, jtape, icase, iorb
    REAL*8     :: pi

    if (myrank.NE.master) RETURN

    pi=ACOS(-1.0D0)

    ! Header for correlated
    do icix=1,ncix
       ! Header
       itape = fh_gc+icix
       jtape = fh_dt+icix
       write(itape,'(A7,14x)',advance='no') '# omega'
       write(jtape,'(A7,14x)',advance='no') '# omega'
       do i=1,csize(icix)
          write(itape,'(A28)',advance='no') legend(i,icix)       !end relativistic DOS
          write(jtape,'(A28)',advance='no') legend(i,icix)       !end relativistic DOS
       enddo
       write(itape,*)
       write(jtape,*)
    enddo

    do iom=1,nomega
       do icix=1,ncix
          ! Header
          itape = fh_gc+icix
          jtape = fh_dt+icix
          write(itape,'(f19.12,1x)',advance='no') omega(iom)*Ry2eV
          write(jtape,'(f19.12,1x)',advance='no') omega(iom)*Ry2eV
          do i=1,csizes(icix)
             write(jtape,'(2f19.12)',advance='no') Deltac(i,icix,iom)*Ry2eV
          enddo
          do i=1,csize(icix)
             write(itape,'(2f19.12)',advance='no') Glc(i,icix,iom)/Ry2eV
          enddo
          write(itape,*)
          write(jtape,*)
       enddo
    enddo


    ! Header for non-correlated
    itape = fh_dos
    write(itape,'(A14,6x)',advance='no') '# omega'
    write(itape,'(A5,9x)',advance='no') 'total'
    do icase=1,natom
       do lcase=1,nl(icase)
          L=ll(icase,lcase)
          write(itape,'(A2,I2,1x,A2,I2,5x)',advance='no') 'a=', iatom(icase), 'L=', L 
       enddo
    enddo
    write(itape,*)
    do iom=1,nomega
       write(itape,'(f14.8,1x)',advance='no') omega(iom)*Ry2eV
       write(itape,'(f14.8)',advance='no') -aimag(gtot(iom))/pi/Ry2eV

       if (allocated(gloc)) then
          do iorb=1,norbitals
             write(itape,'(f14.8)',advance='no') -aimag(gloc(iorb,iom))/pi/Ry2eV
          enddo
       end if

       write(itape,*)
    enddo

    do icix=1,ncix
       close(fh_gc+icix)
       close(fh_dt+icix)
    enddo
    close(fh_dos)
    return
  END SUBROUTINE PrintGloc

end module PrintG
