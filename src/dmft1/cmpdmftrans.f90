! @Copyright 2007 Kristjan Haule
! 

module DMFTProjector
  contains
    SUBROUTINE CmpDMFTrans(DMFTrans, DMFTU, alml, ri_mat, projector, cfX, L, iorb1, cix_orb, nindo, iso, nemin, nbands, nind, maxdim2, icix, idt, icase, lcase, lmax2, nume, nrf, nrfmax, norbitals)
      ! The DMFT transformation is created, which transforms between Kohn-Sham basis and DMFT basis.
      ! It is created by :
      ! DMFTransform(i,j,L,M) = <Y_{L}(r_1)|psi_{ki}(r_1)><psi_{kj}(r_2)|Y_{M}(r_2)>*delta(r_1-r_2)
      ! Where <|> stands for the integral over space angles of r_1 and r_2 only, 
      ! while delta(r_1-r_2) is delta function in radial distance only.
      !
      ! The transformation has the following important properties:
      !
      ! \sum_i DMFTrans(i,i,L1,L2) = delta(L1,L2)
      !
      ! \sum_L DMFTrans(i,j,L,L) = delta(i,j) . Be careful! For this to be valid, one needs to sum over all L to high cutoff. 
      !                                         Something we never really need to do.
      !
      ! sum_{ij} DMFTrans(i,j,L1,L2) * DMFTrans(j,i,L2p,L1p) = delta(L1,L1p)*delta(L2,L2p)
      ! which makes shure that a quantity G_{LL'} when transformed to G_{ij} and back to G_{LL'} takes the same form
      !
      USE p_project,ONLY: rix_mat, al_ucase, al_interstitial
      IMPLICIT NONE
      COMPLEX*16, allocatable, intent(inout) :: DMFTrans(:,:,:,:,:) !DMFTrans(nbands,nbands,maxdim2,maxdim2,norbitals)
      COMPLEX*16, intent(inout) :: DMFTU(nbands,maxdim2,norbitals)
      COMPLEX*16, intent(in)  :: alml(0:lmax2,2*lmax2+1,nume,nrf,2)
      REAL*8,     intent(in)  :: ri_mat(nrf,nrf,0:lmax2,2)
      COMPLEX*16, intent(in)  :: cfX(maxdim2,maxdim2,norbitals,norbitals)
      INTEGER,    intent(in)  :: L, iorb1, cix_orb(norbitals), nindo(norbitals), iso, nemin, nbands, nind, icix, idt
      INTEGER,    intent(in)  :: lmax2, nume, nrf, maxdim2, nrfmax, norbitals
      INTEGER,    intent(in)  :: projector, icase, lcase
      ! locals
      INTEGER :: N1, num1, num2, is1, m1, lms1, ind1, is2, m2, lms2, ind2, irf1, irf2, it, nind1, nind2, iorb2
      COMPLEX*16 :: csum, cc, corc, sm1, sm2, ff
      REAL*8     :: fct
      COMPLEX*16, allocatable :: Rx(:,:), URx(:,:), tmp(:,:)!, cft(:,:)
      INTEGER, allocatable    :: non_zero(:,:,:), all_non_zero(:)
      REAL*8 :: small
      !COMPLEX*16 :: tuu, tud, tdu, tdd
      !COMPLEX*16, allocatable :: Uu(:,:), Ud(:,:)
      INTEGER :: iind


      small = 1e-6   ! cutoff for ri_mat
      N1=2*L+1

      !ALLOCATE( cft(maxdim2,maxdim2) )
      ALLOCATE( tmp(nbands,maxdim2) )

      if (ALLOCATED(DMFTrans)) then
         allocate( non_zero(2,nrf**2,iso), all_non_zero(iso) )
         ! for optimization purposes only ! Many terms in ri_mat are zero, because u*dot{u}=0
         do is1=1,iso
            it=1
            do irf1=1,nrfmax
               do irf2=1,nrfmax
                  if (abs(ri_mat(irf1,irf2,L,is1)).gt.small) then
                     non_zero(1,it,is1) = irf1
                     non_zero(2,it,is1) = irf2
                     it = it+1
                  endif
               enddo
            enddo
            all_non_zero(is1)=it-1
         enddo
         
         !$OMP PARALLEL DO SHARED(DMFTrans) PRIVATE(Rx,num1,num2,is1,m1,lms1,ind1,is2,m2,lms2,ind2,csum,fct,it,irf1,irf2) SCHEDULE(STATIC)
         do num1=1,nbands      !nemin,nemax     ! over bands
            allocate( Rx(nind,nind) )
            do num2=1,nbands   !nemin,nemax     ! over bands
               if (icix.eq.0 .and. num1.ne.num2) CYCLE  ! for non correlated bands, we need only band1=band2 transformation
               do is1=1,iso            ! over spin-1
                  do m1=-L,L           ! over m-1
                     lms1=L+1+m1
                     ind1=L+1+idt*m1+N1*(is1-1)
                     is2=is1
                     do m2=-L,L     ! over m-2
                        lms2=L+1+m2
                        ind2=L+1+idt*m2+N1*(is2-1)                                   
                        csum=0
                        fct = idt**(m1+m2)
                        do it=1,all_non_zero(is1) ! This optimization is because ri_mat is zero for most of components
                           irf1 = non_zero(1,it,is1)
                           irf2 = non_zero(2,it,is1)
                           if (idt.gt.0) then
                              csum = csum + alml(L,lms1,num1+nemin-1,irf1,is1)*conjg(alml(L,lms2,num2+nemin-1,irf2,is2))*(ri_mat(irf1,irf2,L,is1)*fct)
                           else
                              csum = csum + conjg(alml(L,lms1,num1+nemin-1,irf1,is1))*alml(L,lms2,num2+nemin-1,irf2,is2)*(ri_mat(irf1,irf2,L,is1)*fct)
                           endif
                        enddo
                        Rx(ind1,ind2) = csum
                     enddo
                  enddo
               enddo
               DMFTrans(num1,num2,:nind,:nind,iorb1) = Rx 
            enddo
            deallocate( Rx )
         enddo
         !$OMP END PARALLEL DO

         deallocate( non_zero, all_non_zero )
      endif

      if (icix.gt.0) then
         allocate( URx(nbands,nind) )

         if (abs(projector).eq.1) then
            ! For correlated orbitals, we will also create transformation matrix U, such that P(LL'ij)=U^\dagger*U

            !$OMP PARALLEL DO SHARED(URx) PRIVATE(num1,is1,m1,lms1,ind1,cc,irf1) SCHEDULE(STATIC)
            do num1=1,nbands      !nemin,nemax     ! over bands
               do is1=1,iso            ! over spin-1
                  do m1=-L,L           ! over m-1
                     lms1=L+1+m1
                     ind1=L+1+idt*m1+N1*(is1-1)
                     cc=0
                     do irf1=1,nrfmax
                        cc = cc + alml(L,lms1,num1+nemin-1,irf1,is1)*ri_mat(irf1,1,L,is1)*(idt)**m1
                     enddo
                     if (idt.gt.0) then
                        URx(num1,ind1) = cc
                     else
                        URx(num1,ind1) = conjg(cc)
                     endif
                  enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         elseif (abs(projector).eq.2 .or. abs(projector).eq.4) then
            ! For correlated orbitals, we will also create transformation matrix U, such that P(LL'ij)=U^\dagger*U

            !$OMP PARALLEL DO SHARED(URx) PRIVATE(num1,is1,m1,lms1,ind1,cc,sm1,sm2,irf1,irf2,ff) SCHEDULE(STATIC)
            do num1=1,nbands      !nemin,nemax     ! over bands
               do is1=1,iso            ! over spin-1
                  do m1=-L,L           ! over m-1
                     lms1=L+1+m1
                     ind1=L+1+idt*m1+N1*(is1-1)
                     cc=0
                     sm1=0
                     sm2=0
                     do irf1=1,nrfmax
                        do irf2=1,nrfmax
                           ff = alml(L,lms1,num1+nemin-1,irf1,is1)*conjg(alml(L,lms1,num1+nemin-1,irf2,is1))
                           sm1 = sm1 + ff * ri_mat(irf1,irf2,L,is1)
                           sm2 = sm2 + ff * ri_mat(irf1,1,L,is1) * ri_mat(1,irf2,L,is1)
                        enddo
                        cc = cc + alml(L,lms1,num1+nemin-1,irf1,is1)*ri_mat(irf1,1,L,is1)*(idt)**m1
                     enddo
                     cc = cc * sqrt(abs(sm1/sm2))
                     if (idt.gt.0) then
                        URx(num1,ind1) = cc
                     else
                        URx(num1,ind1) = conjg(cc)
                     endif
                  enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         elseif (abs(projector).eq.3) then
            ! For correlated orbitals, we will also create transformation matrix U, such that P(LL'ij)=U^\dagger*U

            !$OMP PARALLEL DO SHARED(URx) PRIVATE(num1,is1,m1,lms1,ind1,cc,corc,irf1,irf2) SCHEDULE(STATIC)
            do num1=1,nbands      !nemin,nemax     ! over bands
               do is1=1,iso            ! over spin-1
                  do m1=-L,L           ! over m-1
                     lms1=L+1+m1
                     ind1=L+1+idt*m1+N1*(is1-1)
                     irf1 = 1
                     cc = alml(L,lms1,num1+nemin-1,irf1,is1)*sqrt(ri_mat(irf1,irf1,L,is1))*(idt)**m1
                     ! Adding correction due to \cdot{u} and local orbitals
                     ! The results is  A(i,L,irf=1)*sqrt(o(L)) * sqrt( 1 + \sum_{irf>1} |A(i,L,irf)|^2*o(L,irf))
                     if ( abs(cc)>1e-6 ) then
                        corc = 0
                        do irf1=1,nrfmax
                           do irf2=1,nrfmax
                              if (irf1.eq.1 .and. irf2.eq.1) cycle
                              corc = corc + alml(L,lms1,num1+nemin-1,irf1,is1)*conjg(alml(L,lms1,num1+nemin-1,irf2,is1))*ri_mat(irf1,irf2,L,is1)
                           enddo
                        enddo
                        cc = cc*sqrt(1 + abs(corc)/abs(cc)**2)
                     endif
                     if (idt.gt.0) then
                        URx(num1,ind1) = cc
                     else
                        URx(num1,ind1) = conjg(cc)
                     endif
                  enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         elseif (abs(projector).eq.5 .or. abs(projector).eq.6) then
            iind = al_ucase(icase,lcase)
            !$OMP PARALLEL DO SHARED(URx,rix_mat) PRIVATE(num1,is1,m1,lms1,ind1,cc,irf1) SCHEDULE(STATIC)
            do num1=1,nbands  
               do is1=1,iso   
                  do m1=-L,L  
                     lms1=L+1+m1
                     ind1=L+1+idt*m1+N1*(is1-1)
                     cc=0
                     do irf1=1,nrfmax
                        cc = cc + alml(L,lms1,num1+nemin-1,irf1,is1)*rix_mat(irf1,iind)*(idt)**m1
                     enddo
!!! New
                     cc = cc + al_interstitial(lms1,num1,is1,lcase)*(idt)**m1
!!! New
                     if (idt.gt.0) then
                        URx(num1,ind1) = cc
                     else
                        URx(num1,ind1) = conjg(cc)
                     endif
                  enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         else
            print *, 'Only projector=[1,2,3,4,5,6] is allowed!'
            stop
         endif

         nind1 = (2*L+1)*iso
         !if (iso.eq.2) then
         !   no1 = nind1/iso
         !   ALLOCATE( Uu(nbands,no1), Ud(nbands,no1))
         !   tuu = Rspin(1,1,iorb1)
         !   tud = Rspin(1,2,iorb1)
         !   tdu = Rspin(2,1,iorb1)
         !   tdd = Rspin(2,2,iorb1)
         !   Uu(:,:) = URx(:,1:no1)
         !   Ud(:,:) = URx(:,no1+1:2*no1)
         !   URx(:,1:no1)       = Uu*tuu + Ud*tdu
         !   URx(:,no1+1:2*no1) = Uu*tud + Ud*tdd
         !   DEALLOCATE( Uu, Ud )
         !endif

         DO iorb2=1,norbitals
            if ( cix_orb(iorb2).NE. icix ) CYCLE
            nind2 = nindo(iorb2)
            call zgemm('N','C', nbands, nind2, nind1, (1.d0,0.d0), URx,nbands, cfX(:,:,iorb2,iorb1),maxdim2, (0.d0,0.d0), tmp,nbands)
            DMFTU(:,:nind2,iorb2) = DMFTU(:,:nind2,iorb2) + conjg(tmp(:,:nind2))
         ENDDO

         deallocate(URx)
      endif

      !deallocate( cft )
      deallocate( tmp )
    end SUBROUTINE CmpDMFTrans
end module DMFTProjector
