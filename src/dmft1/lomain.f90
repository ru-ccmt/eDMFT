SUBROUTINE LoMAIN(nemin,nemax,lfirst,latom,n,jatom,isym,LC,iso,crotloc_x_BR1)
  USE param
  USE structure, ONLY: mult, rotij, tauij, pos, tau
  USE sym2
  USE abc
  USE case
  USE matpar,   only: atpar, alo, nlo, nlon, ilo
  IMPLICIT NONE
  INTEGER, intent(in) :: nemin, nemax, lfirst, latom, n, jatom, isym, LC, iso
  REAL*8,  intent(in) :: crotloc_x_BR1(3,3)
  !IMPLICIT REAL*8 (A-H,O-Z)
  COMPLEX*16 :: YL((LMAX2+1)*(LMAX2+1))
  COMPLEX*16 :: PHSHEL,PH_SPIN(2)
  COMPLEX*16,ALLOCATABLE   :: phs(:)
  DATA           IMAG/(0.0D0,1.0D0)/         
  COMPLEX*16 :: IMAG, cfac
  REAL*8     :: PI, TWOPI, ARG123, ARG2, ARGT, ARGT2
  real*8     :: BKROT(3), BKROT2(3), BKRLOC(3)
  INTEGER    :: I, L, jlo, jneq, M1, is, NUM, irf
!-----------------------------------------------------------------------------     
  PI=ACOS(-1.0D0)
  TWOPI=2.D0*PI
  ALLOCATE(phs(nume))
  i=n-(nlo+nlon)
  DO L=0,LoMAX
     do jlo=1,ilo(l)
        do jneq=1,mult(jatom)
           DO M1=-l,+l
              i=i+1  
              if (.not.(L.eq.LC)) CYCLE
              ! BKROT = R_a.(k+K) transforms to the reducible k-point
              BKROT = matmul(TMAT(:,:,isym), BK3(:,I) )
              ! BKROT2 = R_n.R_a.(k+K), transformation from the first atom to an equivalent atom               
              BKROT2 = matmul(rotij(:,:,latom), BKROT)
              !---- BR1 transforms integer reciprocal lattice vectors, as given in the VECTORLIST of LAPW1, into cartesian system ----!
              ! BKROT3 = R_n.R_a.(k+K), but in cartesian coordinate system
              !BKROT3 = matmul(BR1, BKROT2)
              !!---- BKRLOC = crotloc.R_n.R_a.(k+K),  rotates according to the user specified local coordinate system.
              !BKRLOC = matmul(crotloc(:,:,icase), BKROT3)
              BKRLOC = matmul(crotloc_x_BR1, BKROT2)
              !---- YLM = Y_{L}(Rotloc.R_g.(k+K))
              CALL YLM (BKRLOC,LMAX2,YL)
              ! (R_n.R_a.(k+K)) *  R(first) * 2pi
              ! ARG123 = (R_g.(k+K)) *  R(latom) * 2pi
              ARG123 = dot_product(BKROT2, POS(:,lfirst))*TWOPI 
              ! ARGT = (k+K)*tau(isym) * 2pi
              ARGT = dot_product(BK3(:3,I), TAU(:3,isym))*TWOPI
              ! ARGT2 = (R_a.(k+K)).tau_n * 2pi
              ARGT2 = dot_product(BKROT, tauij(:3,latom))*TWOPI
              ! ARG2 = (R_a.(k+K)) *  shft * 2pi
              ARG2  = dot_product(BKROT, shft(:3,latom))*TWOPI
              ! PHSEHL = e^{I*2pi*( (R_a.(k+K))*tau_n + (K+k)*tau(isym) + (R_n.R_a.(k+K)*R(first)))}
              PHSHEL=EXP(IMAG*(ARG123+ARG2+ARGT+ARGT2))
              DO is=1,iso
                 PH_SPIN(IS)=EXP((2*is-3)*IMAG*PHASE(ISYM)/2)
                 PHS(nemin:nemax) = A(I,nemin:nemax,is)*PHSHEL
                 DO irf=1,nrf
                    DO num=nemin,nemax
                       cfac = ALo(l,jlo,irf,is)*PHS(num)*PH_SPIN(is)
                       alm(1:(2*l+1),num,irf,is)=alm(1:(2*l+1),num,irf,is) + cfac * dconjg(YL((l*l+1):(l+1)**2))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        enddo
     enddo
  ENDDO
  return                   
END SUBROUTINE LoMAIN
