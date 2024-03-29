SUBROUTINE mmatcv(l,iat,ivmax,nemin,nemax,nnj,mm,is,lso,spin,nk)
  use comi, nat=>natti
  use mxyz
  use moments
  use ablm, alfa=>alm,beta=>blm,gama=>clm
  use xrpar
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
  !......COMPUTING : < PSI'_core| GRAD | PSI >...............................
  !................. | PSI > := (ALM(NB) UL + BLM(NB) UPL) YLM .........
  !................. | PSI_CORE> := (UL_CORE) YLM .........
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  !-------------------------------------------------------------------
  IMPLICIT REAL*8 (A-H,O-Z)
  INCLUDE 'param.inc'
  integer(4) :: par,z,zmin,zmax
  integer(4), intent(in) :: nnj,mm 
  real(8) :: CB,CCBB  !,CB1,CB2
  real(8) :: lll,mmm,la,ma,m2,ms,msup,msdn,mmm1,lll1
  integer(4), intent(in) :: iat,nemin
  COMPLEX*16  CZERO,IMAG
  COMPLEX*16 :: P_1(3)
  
  logical, intent(in) :: lso, spin
  integer (4) :: denom1, denom2, lm, l1mm1, l1m, l1m1
  integer (4) :: lm1mm1, lm1m, lm1m1       
  real (8) :: flm1, flm2, flm3, flm4, flm5, flm6
  complex*16 :: xpy1_1, xmy1_1, z1_1
  complex*16 :: xpy2_1, xmy2_1, z2_1
  
  !ad 
  DATA    CZERO/(0.0D0,0.0D0)/IMAG/(0.0D+0,1.0D+0)/ 
  DATA    ONE/1.0D+0/TWO/2.0D+0/THREE/3.0D+0/
  
  ms = 0.5
  msup = 0.5
  msdn = -0.5
  
  !........INLINE FUNCTIONS FOR YLM OTHOGONALITY ....................
  !      
  !       F0(LC)   = (TWO*LC+ONE)*(TWO*LC+THREE)
  !       F1(LC,M) = - SQRT(( LC + M + ONE ) * ( LC + M + TWO ) / F0(LC))
  !       F3(LC,M) =   SQRT(( LC - M + ONE ) * ( LC - M + TWO ) / F0(LC))
  !       F5(LC,M) =   SQRT(( LC - M + ONE ) * ( LC + M + ONE ) / F0(LC))
  !..................................................................
  
  do iv=nemin,nemax
     MX_(iv)=CZERO
     MY_(iv)=CZERO
     MZ_(iv)=CZERO
     pxpy_1(iv)=CZERO
     pxmy_1(iv)=CZERO
     pz_1(iv)=CZERO
     mxcv_1=CZERO
     mycv_1=CZERO
     mzcv_1=CZERO	 
  enddo
    
  !ad..........................................................M=-L,L 

  CB=0            
  DO M=-L,L
     
     lll1=real(nnj)              
     mmm1=real(mm)
     lll=lll1/2
     mmm=mmm1/2
     la=real(l)
     ma=real(m)      
     
     if(lso.and.spin) then
        
        if(is.eq.1)then
           m2=msup
           CB = cleb(la,ms,lll,ma,msup,mmm)
           !            write(*,*)'j1=',la,'j2=',ms,'j3=',lll, & 
           !             &     'm1=',ma,'m2=',msup,'m3=',mmm,'CBup=',CB
        else
           m2= msdn
           CB = cleb(la,ms,lll,ma,msdn,mmm)
        endif
     else
        write(6,*) 'XMCD calculation requires spinpolarized AND spin-orbit setup'
        stop
     endif
     
     do iv=nemin,nemax		
        denom1=dble((2*l+1)*(2*l+3))
        denom2=dble((2*l-1)*(2*l+1))
        !       lm -> l,m                
        lm = l*l+l+m+1 
        !       l1mm1 -> l+1,m-1                
        l1mm1=lm+2*l+1
        !       l1m -> l+1,m                
        l1m = l1mm1 + 1
        !       l1m1 -> l+1,m+1                
        l1m1 = l1m + 1
        !       lm1mm1 -> l-1,m-1                
        lm1mm1=lm-2*l-1
        !       lm1m -> l-1,m                
        lm1m = lm1mm1 + 1
        !       lm1m1 -> l-1,m+1                
        lm1m1 = lm1m + 1
        !       calculate the F_lm coefficients
        flm1=-sqrt(dble((l+m+1)*(l+m+2))/denom1)
        flm2=sqrt(dble((l-m-1)*(l-m))/denom2)
        flm3=sqrt(dble((l-m+1)*(l-m+2))/denom1)
        flm4=-sqrt(dble((l+m-1)*(l+m))/denom2)
        flm5=sqrt(dble((l-m+1)*(l+m+1))/denom1)
        flm6=sqrt(dble((l-m)*(l+m))/denom2)
        !       calculate the matrix elements <core/p/valence>
        xpy1_1=CZERO
        xmy1_1=CZERO
        z1_1=CZERO
        if (l.gt.0)then
           if (l+m-2.ge.0) then 
              xpy1_1=alfa(iv,lm1mm1)*cmplx(iucl1ul(iat,is),0.0d0,8)+ beta(iv,lm1mm1)*cmplx(iucl1udl(iat,is),0.0d0,8) + gama(iv,lm1mm1)*cmplx(iucl1ulol(iat,is),0.0d0,8)
           endif
           if (l-m-2.ge.0)then
              xmy1_1=alfa(iv,lm1m1)*cmplx(iucl1ul(iat,is),0.0d0,8)+ beta(iv,lm1m1)*cmplx(iucl1udl(iat,is),0.0d0,8)+ gama(iv,lm1m1)*cmplx(iucl1ulol(iat,is),0.0d0,8)
           endif
           if ((l+m-1.ge.0).and.(l-m-1.ge.0))then
              z1_1=alfa(iv,lm1m)*cmplx(iucl1ul(iat,is),0.0d0,8)+ beta(iv,lm1m)*cmplx(iucl1udl(iat,is),0.0d0,8)+ gama(iv,lm1m)*cmplx(iucl1ulol(iat,is),0.0d0,8) 
           endif
        endif
        xpy2_1=alfa(iv,l1mm1)*cmplx(iuclul1(iat,is),0.0d0,8)+ beta(iv,l1mm1)*cmplx(iucludl1(iat,is),0.0d0,8) + gama(iv,l1mm1)*cmplx(iuclulol1(iat,is),0.0d0,8)
        xmy2_1=alfa(iv,l1m1)*cmplx(iuclul1(iat,is),0.0d0,8)+ beta(iv,l1m1)*cmplx(iucludl1(iat,is),0.0d0,8) + gama(iv,l1m1)*cmplx(iuclulol1(iat,is),0.0d0,8)
        z2_1=alfa(iv,l1m)*cmplx(iuclul1(iat,is),0.0d0,8)+ beta(iv,l1m)*cmplx(iucludl1(iat,is),0.0d0,8)+ gama(iv,l1m)*cmplx(iuclulol1(iat,is),0.0d0,8)
        
        !       we make use of the fact that:
        !       F^1-(l-1,m-1)=f^4_{lm}
        !       F^2-(l-1,m-1)=f^3_{lm}
        !       F^3-(l-1,m-1)=f^2_{lm}
        !       F^4-(l-1,m-1)=f^1_{lm}
        !       F^5-(l-1,m-1)=f^6_{lm}
        !       F^6-(l-1,m-1)=f^5_{lm}
        
        pxpy_1(iv) = pxpy_1(iv) + CB*(flm4 * xpy1_1 + flm3 * xpy2_1)   ! grad_x+i*grad_y
        pxmy_1(iv) = pxmy_1(iv) + CB*(flm2 * xmy1_1 + flm1 * xmy2_1)   ! grad_x-i*grad_y
        pz_1(iv)   = pz_1(iv) + CB*(flm6 * z1_1   + flm5 * z2_1)       ! grad_z
        mxcv_1(iv)=5.0d-1*(pxpy_1(iv)+pxmy_1(iv))                      ! grad_x
        mycv_1(iv)=5.0d-1*imag*(pxpy_1(iv)-pxmy_1(iv)) !THE REMOVED A leading minus...  ! -grad_y
        mzcv_1(iv)=pz_1(iv)                                            ! grad_z
     enddo
  enddo
  
  !ad................................................end.......M=-L,L
  !ad    Transformation to Cartesian coordinates
  !ad
  do iv=nemin,nemax
     p_1(1)=mxcv_1(iv)
     p_1(2)=mycv_1(iv)
     p_1(3)=mzcv_1(iv)
     MX_(iv)=p_1(1)
     MY_(iv)=p_1(2)
     MZ_(iv)=p_1(3)
     p_1(1)=CZERO           
     p_1(2)=CZERO           
     p_1(3)=CZERO
  enddo
  RETURN
end SUBROUTINE mmatcv

real*8 function cleb(j1,j2,j3,m1,m2,m3)
  Implicit None
  real(8), intent (in) :: j1,j2,j3,m1,m2,m3
  integer(4) :: jj1, jj2, jj3, mm1, mm2, mm3
  integer(4) :: ka1, ka2, ka3, i
  Real (8) :: ThreeJSymbol
  External ThreeJSymbol
  jj1 = nInt(2*j1)
  jj2 = nint(2*j2)
  jj3 = nint(2*j3)
  mm1 = nInt(2*m1)
  mm2 = nint(2*m2)
  mm3 = nint(2*m3)
  !-------check arguments-------------------
  ka1=0
  ka2=0
  ka3=0
  do i=-abs(jj1),jj1,2
     if(mm1.eq.i) ka1=ka1+1
  enddo
  do i=-abs(jj2),jj2,2
     if(mm2.eq.i) ka2=ka2+1
  enddo
  do i=-abs(jj3),jj3,2
     if(mm3.eq.i) ka3=ka3+1
  enddo
  
  if((ka1.ne.1).or.(ka2.ne.1).or.(ka3.ne.1)) then
     Write (6, '("Error(clebgor): non-physical arguments!")')
     stop
  endif
  !-------end check------------------
  If ((jj1 .Lt. 0) .Or. (jj2 .Lt. 0) .Or. (jj3 .Lt. 0) .Or. (Abs(m1).Gt. j1) .Or. (Abs(m2) .Gt. j2) .Or. (Abs(m3) .Gt. j3)) Then
     Write (6, '("Error(clebgor): non-physical arguments!")')
     Stop
  End If
  If ((jj1 .Eq. 0) .And. (jj2 .Eq. 0) .And. (jj3 .Eq. 0)) Then
     cleb = 1.d0
     Return
  End If
  If ((jj1 .Gt. 100) .Or. (jj2 .Gt. 100) .Or. (jj3 .Gt. 100)) Then
     Write (6,*)
     Write (6, '("Error(clebgor): angular momenta out of range : ",3F5.2)') j1, j2, j3
     Write (6,*)
     Stop
  End If
  If ((mm1+mm2-mm3 .Ne. 0) .Or. (j2+j3 .Lt. j1) .Or. (j1+j3 .Lt. j2).Or. (j1+j2 .Lt. j3)) Then
     cleb = 0.d0
     Return
  End If
  cleb = Sqrt (dble(2*j3+1)) * ThreeJSymbol (j1, j2, j3, m1, m2,-m3)
  If (Mod(nint(j1-j2+m3), 2) .Ne. 0) cleb = - cleb
  Return
end function cleb


real*8 function ThreeJSymbol(L1, L2, L3, M1, M2, M3)
  ! !INPUT/OUTPUT PARAMETERS:
  !   L1,L2,L3,M1,M2,M3   :   integers specifying the 3j-symbol (L1   L2   L3)
  !                                                              M1   M2   M3
  ! !DESCRIPTION:
  !   Calculates Wigner 3j symbols using the Racah formula.
  ! !REVISION HISTORY:
  !   Updated November 2004 (Kevin Jorissen)
  !EOP
  implicit none
  real(8), intent(in) :: l1,l2,l3,m1,m2,m3
  integer(4) :: jj1, jj2, jj3, mm1, mm2, mm3
  integer k,kmin,kmax
  real*8 tjs
  real*8, external :: fact
  jj1 = nInt(2*l1)
  jj2 = nint(2*l2)
  jj3 = nint(2*l3)
  mm1 = nInt(2*m1)
  mm2 = nint(2*m2)
  mm3 = nint(2*m3)
  if ((abs(m1).gt.l1).or.(abs(m2).gt.l2).or.(abs(m3).gt.l3)) then
     ThreeJSymbol = dble(0)
     return
  else
     tjs = dble(0)
     kmin = max0(0, -nint(l3-l2+m1), -nint(l3-l1-m2))
     kmax = min0(nint(l1+l2-l3), nint(l1-m1), nint(l2+m2))
     do k = kmin, kmax
        tjs = tjs + ((-1)**k) / fact(k) / fact(nint(l1+l2-l3-k)) / fact(nint(l1-m1-k)) / fact(nint(l2+m2-k)) / fact(nint(l3-l2+m1+k)) / fact(nint(l3-l1-m2+k))
     enddo
     tjs = tjs * (-1)**(nint(l1-l2-m3))
     tjs = tjs * dsqrt (fact(nint(l1+l2-l3)) * fact(nint(l1-l2+l3)) * fact(nint(-l1+l2+l3)) / fact(nint(l1+l2+l3+1)) * fact(nint(l1+m1)) * fact(nint(l1-m1)) * fact(nint(l2+m2)) * fact(nint(l2-m2))* fact(nint(l3+m3)) * fact(nint(l3-m3)) )
     ThreeJSymbol = tjs
  endif
  return
end function ThreeJSymbol


REAL*8 FUNCTION FACT(N)
  ! !INPUT/OUTPUT PARAMETERS:
  !   n   :  argument
  ! !DESCRIPTION:
  !   Calculates a factorial.
  ! !REVISION HISTORY:
  !     Updated November 2004 (Kevin Jorissen)
  !EOP
  implicit none
  !  INPUT
  integer,intent(in) :: n
  !  LOCALS
  INTEGER N1
  REAL*8 FA
  
  FA = dble(1)
  IF (N.GT.0) THEN
     DO N1 = 1, N
        FA = FA * DBLE(N1)
     ENDDO
  ENDIF
  FACT = FA
  RETURN
END FUNCTION FACT
