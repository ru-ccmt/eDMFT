! @Copyright 2007 Kristjan Haule
! 

module p_project
  IMPLICIT NONE
  REAL*8, allocatable ::  P_rfi(:,:)
  REAL*8, allocatable :: rix_mat(:,:)
  INTEGER :: n_al_ucase, max_lcase
  INTEGER, allocatable :: al_ucase(:,:), al_ucase_(:,:), l_al_ucase(:), j_al_ucase(:), nl_case(:)
  INTEGER, parameter :: kNri=5  ! Number of points for the wave function in the interstitials is 2**kNri+1
  REAL*8, allocatable :: rix(:), dri(:)
  !COMPLEX*16, allocatable :: al_interstitial(:,:,:,:)
  ! locals
  INTEGER :: max_nrad, nrf
  REAL*8, allocatable ::  P_rfk(:,:,:)
  REAL*8, allocatable ::  phi_jl(:, :)
           
CONTAINS
    
  SUBROUTINE p_allocate(P_filename)
    use param, ONLY: nrad, lmax2, nloat, lomax
    use com_mpi, ONLY : myrank, master, Bcast_Size, Bcast_Projector
    use structure, ONLY: rmt, jri, r0, mult
    use com, ONLY: nat
    USE dmfts, ONLY: iso, natom, iatom, isort, cix, ll, nl !, iso, projector, nl, ll, cix, crotloc, ncix, maxdim, maxsize
    !USE w_atpar, ONLY: w_ilo, w_loor
    USE lo,    ONLY: elo_store, elo, rlo, loor, lapw, nlo, nlov, nlon, ilo
    USE atspdt,ONLY: e_store, el
    IMPLICIT NONE
    CHARACTER*100, intent(in) :: P_filename
    ! locals
    INTEGER :: icase, iucase, il, lcase, fhp, ncase, ii, ir, nr, nr0, jatom, lc, l, i, Nri, lmax, nrf_
    CHARACTER*2 :: cm
    INTEGER, allocatable :: tmp_ind(:,:)
    REAL*8,  allocatable :: rx(:), ypp(:)
    !INTEGER, allocatable :: nrj(:)
    REAL*8 :: r, yval, ypval, yppval
    INTEGER:: csort(nat), maxucase, jlo
    fhp=1001

    ! determines maxucase and csort
    csort=0
    maxucase=1
    do icase=1,natom
       if (csort(isort(iatom(icase))).EQ.0) then
          csort(isort(iatom(icase))) = maxucase
          maxucase = maxucase + 1
       endif
    enddo
    maxucase = maxucase-1
    
    
    allocate( nl_case(maxucase) )
    
    max_lcase=0
    DO icase=1,natom
       iucase = csort(isort(iatom(icase)))  ! some atoms are equivalent, creates unique index to inequivalent correlated atoms
       il=0
       do lcase=1,nl(icase)
          if (cix(icase,lcase).GT.0) il=il+1 ! only correlated states will count
       enddo
       nl_case(iucase)=il                    ! we have at most nl_case correlated l's. This is usually just one, but we will keep general possibility
       if (il.gt.max_lcase) max_lcase=il
    ENDDO
    
    allocate( al_ucase_(maxucase,max_lcase) )  ! contains unique index for correlated orbital. There are less unique correlated orbitals than there is norbital, because
    n_al_ucase=0                               ! the equivalent atoms are counted multiple times in norbital index
    do iucase=1,maxucase                       ! unique correlated atoms
       do il=1,nl_case(iucase)                 ! correlated l's on thie unique atom
          n_al_ucase = n_al_ucase+1            ! unique index for (atom,il)
          al_ucase_(iucase,il) = n_al_ucase    ! save this unique index for later
       enddo
    enddo

    lmax=0
    allocate( l_al_ucase(n_al_ucase), j_al_ucase(n_al_ucase) )       ! gives l given unique index (atom,il)
    allocate( al_ucase(natom,lmax2+1) )     ! given (icase,lcase) it gives unique index of unique correlated orbital
    al_ucase(:,:)=0
    do icase=1,natom
       iucase = csort(isort(iatom(icase)))   ! some atoms are equivalent, creates unique index to inequivalent correlated atoms
       il=0
       do lcase=1,nl(icase)
          if (cix(icase,lcase).GT.0) then
             il=il+1
             ii = al_ucase_(iucase,il)            ! unique index created above
             al_ucase(icase,lcase) = ii           ! unique index saved to al_ucase
             l_al_ucase(ii) = ll(icase,lcase)     ! l can be accessed from unique index
             j_al_ucase(ii) = isort(iatom(icase)) ! jatom
             if ( ll(icase,lcase) .gt. lmax ) lmax = ll(icase,lcase)
          endif
       enddo
    enddo
    
    
    if (myrank.eq.master) then
       WRITE(6,'(A,I3)') 'max_lcase=', max_lcase
       do iucase=1,maxucase
          WRITE(6,'(A,I3,A,I3)') 'nl_case(', iucase, ')=', nl_case(iucase)
       enddo
       do icase=1,natom
          do lcase=1,nl(icase)
             if (al_ucase(icase,lcase).GT.0) WRITE(6,'(A,I3,A,I3,A,I3)') 'al_ucase(', icase, ',', lcase, ')=', al_ucase(icase,lcase)
          enddo
       enddo
       WRITE(6,'(A,I3)') 'n_al_ucase=', n_al_ucase
       do ii=1,n_al_ucase
          WRITE(6,'(A,I3,A,I3)') 'l_al_ucase(', ii,')=', l_al_ucase(ii)
          WRITE(6,'(A,I3,A,I3)') 'j_al_ucase(', ii,')=', j_al_ucase(ii)
       enddo
    end if

    Nri = 2**kNri+1

    !allocate( nrj(n_al_ucase) )
    allocate( P_rfi(Nri,n_al_ucase) )
    allocate( dri(n_al_ucase) )

    if (myrank.eq.master) then
       open(fhp, FILE=TRIM(P_filename), status='old', ERR=900)
       READ(fhp,*) cm, ncase, max_nrad
       
       if (ncase.NE.n_al_ucase) WRITE(6,*) 'ERROR: number of correlated orbitals from projectorw.dat and n_al_ucase do not agree!'
       
       allocate( P_rfk(max_nrad,2,n_al_ucase) )
       
       do ii=1,n_al_ucase
          l = l_al_ucase(ii)

          READ(fhp,*) cm, nr , nr0, jatom, lc
          
          if (nr0 .NE. jri(jatom)) WRITE(6,*) 'ERROR: number of radial points from structure and projectorw.dat do not agree',nr0, jri(jatom)
          if (lc .NE.l) WRITE(6,*) 'ERROR: The orbital momentum l from projectorw.dat and case.indmfl do not agree', lc, l

          !nrj(ii)=nr
          
          allocate(rx(nr))
          do ir=1,nr
             READ(fhp,*) rx(ir), P_rfk(ir,1,ii), P_rfk(ir,2,ii)
          enddo
          if ( abs(rx(1)-r0(jatom)).GT.1e-10 ) WRITE(6,*) 'ERROR: first point in logarithmic mesh from structure and projectorw.dat do not agree', rx(1), r0(jatom)
          if ( abs(rx(nr0)-rmt(jatom)).GT.1e-6) WRITE(6,*) 'ERROR: the muffin-thin radius from structure and projectorw.dat do not agree', rx(nr0), rmt(jatom)
          
          dri(ii) = (rx(nr)-rx(nr0))/(Nri-1.)

          if (nr.gt.nr0) then

             allocate( ypp(nr) )
             CALL spline_cubic_set (nr, rx, P_rfk(:,1,ii), 0, 0.0D0, 0, 0.0D0, ypp )
             r=rx(nr0)
             do i=1,Nri
                CALL spline_cubic_val (nr, rx, P_rfk(:,1,ii), ypp, r, yval, ypval, yppval )
                P_rfi(i,ii)=yval
                r=r+dri(ii)
             enddo
             deallocate( ypp )
          else
             P_rfi(:,ii)=0.0
          end if

          !open(199, file='p_inters.dat',status='unknown')
          !r=rx(nr0)
          !do i=1,Nri
          !   WRITE(199,*) r, P_rfi(i,ii)
          !   r=r+dri(ii)
          !enddo
          !close(199)
          
          deallocate(rx)
       enddo
       close(fhp)
    endif

    
    CALL Bcast_Size(max_nrad)

    if (myrank.ne.master) allocate( P_rfk(max_nrad,2,n_al_ucase) )
    
    CALL Bcast_Projector(P_rfk, P_rfi, dri, max_nrad, n_al_ucase, Nri)
    
    allocate( rix(Nri) )

    nrf=2
    do icase=1,natom
       jatom = isort(iatom(icase))
       
       ! This computes on the fly the following quantities: el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo
       CALL get_ilo(e_store,elo_store,jatom,mult,nat,lmax2,lomax,nloat, el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo)
       do lcase=1,nl(icase)
          if (cix(icase,lcase).GT.0) then
             l = ll(icase,lcase)
             nrf_=2
             DO jlo=1,ilo(l)
                if (loor(jlo,l)) nrf_ = 2+jlo !nrf_=nrf_+1
             ENDDO
             if (nrf_.gt.nrf) nrf=nrf_
          endif
       enddo
    enddo


    if (myrank.eq.master) then
       WRITE(6,'(A,I3)') 'nrf=', nrf
    endif
    
    allocate( rix_mat(nrf,n_al_ucase) )
    
    RETURN
900 CONTINUE
    print *, 'Could not find file for projector wave function ', TRIM(P_filename), '!'
    WRITE(6,*) 'Could not find file for projector wave function ', TRIM(P_filename), '!'
    STOP
  END SUBROUTINE p_allocate
  
  SUBROUTINE p_cmp_rixmat()
    use structure, ONLY: jri, r0, dx, mult
    use param, ONLY: nrad, lmax2, lomax, nloat
    use com, ONLY: rel, nat
    use com_mpi, ONLY : myrank, master
    USE w_atpar, ONLY: w_RRAD1, w_RRAD2, w_RADE1, w_RADE2, w_a1lo, w_b1lo!, w_loor, w_ilo
    !
    USE lo,    ONLY: elo_store, elo, rlo, loor, lapw, nlo, nlov, nlon, ilo
    USE atspdt,ONLY: e_store, el
    IMPLICIT NONE
    ! External function
    REAL*8  :: Rint13n
    ! locals
    INTEGER :: ii, l, jatom, irf, nr0, jlo, i
    REAL*8  :: r0_, dx_, r
    
    rix_mat(:,:)=0

    do ii=1,n_al_ucase
       l     = l_al_ucase(ii)
       jatom = j_al_ucase(ii)
       nr0 = jri(jatom)
       r0_ = r0(jatom)
       dx_ = dx(jatom)

       !rfk(:nr0,1,1) = w_RRAD1(:nr0,l,jatom)
       !rfk(:nr0,2,1) = w_RRAD2(:nr0,l,jatom)
       rix_mat(1,ii) = Rint13n(nr0, w_RRAD1(:nr0,l,jatom), w_RRAD2(:nr0,l,jatom), P_rfk(:nr0,1,ii), P_rfk(:nr0,2,ii), r0_, dx_, rel)
       
       !rfk(:nr0,1,2) = w_RADE1(:nr0,l,jatom)
       !rfk(:nr0,2,2) = w_RADE2(:nr0,l,jatom)
       rix_mat(2,ii) = Rint13n(nr0, w_RADE1(:nr0,l,jatom), w_RADE2(:nr0,l,jatom), P_rfk(:nr0,1,ii), P_rfk(:nr0,2,ii), r0_, dx_, rel)

       !rix_mat(2,ii) = Rint13n(nr0, w_RADE1(:nr0,l,jatom), w_RADE2(:nr0,l,jatom),w_RRAD1(:nr0,l,jatom), w_RRAD2(:nr0,l,jatom), r0_, dx_, rel)

       ! This computes on the fly the following quantities: el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo
       CALL get_ilo(e_store,elo_store,jatom,mult,nat,lmax2,lomax,nloat, el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo)

       DO jlo=1,ilo(l)
          if (loor(jlo,l)) then
             !rfk(:nr0,1,2+jlo) = w_a1lo(:nr0,jlo,l,jatom)
             !rfk(:nr0,2,2+jlo) = w_b1lo(:nr0,jlo,l,jatom)
             rix_mat(2+jlo,ii) = Rint13n(nr0, w_a1lo(:nr0,jlo,l,jatom),  w_b1lo(:nr0,jlo,l,jatom), P_rfk(:nr0,1,ii), P_rfk(:nr0,2,ii), r0_, dx_, rel)
          endif
       ENDDO
       
       !do irf=1,nrf
       !   rix_mat(irf,ii) = Rint13n(nr0, w_rfk(:nr0,1,irf,ii), w_rfk(:nr0,2,irf,ii), P_rfk(:nr0,1,ii), P_rfk(:nr0,2,ii), r0_, dx_, rel)
       !enddo
    enddo

    if (myrank.eq.master) then
       WRITE(6,*)
       WRITE(6,'(A)') 'Overlap of projector function (from projectorw.dat) and the solution of Dirac equation'
       do ii=1,n_al_ucase
          do irf=1,nrf
             WRITE(6,'(A,I3,A,I2,A,F18.15)') '    rix_mat(',ii, ',irf=',irf,')=', rix_mat(irf,ii)
          enddo
       enddo
    endif

    !open(99,file="DEBUG.brisi", status='unknown')
    !do ii=1,n_al_ucase
    !   l     = l_al_ucase(ii)
    !   jatom = j_al_ucase(ii)
    !   nr0 = jri(jatom)
    !   r0_ = r0(jatom)
    !   dx_ = dx(jatom)
    !   do i=1,nr0
    !      r = r0_*exp(dx_*(i-1))
    !      WRITE(99,'(7f15.10)') r, P_rfk(i,1,ii), P_rfk(i,2,ii), w_RRAD1(i,l,jatom), w_RRAD2(i,l,jatom), w_RADE1(i,l,jatom), w_RADE2(i,l,jatom)
    !   enddo
    !enddo
          
  end SUBROUTINE p_cmp_rixmat
  
  SUBROUTINE p_deallocate()
    IMPLICIT NONE
    deallocate( al_ucase, l_al_ucase, j_al_ucase )
    deallocate( P_rfk, P_rfi )
    deallocate( dri )
    deallocate( nl_case )
    deallocate( al_ucase_ )
    deallocate( rix_mat )
    !deallocate( nrj )
    deallocate( rix )
  END SUBROUTINE p_deallocate
end module p_project
