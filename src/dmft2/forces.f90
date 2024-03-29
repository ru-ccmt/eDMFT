! @Copyright 2007 Kristjan Haule
! 

module Forces
  logical  :: Qforce                   ! whether we are computing forces at all
  logical, allocatable  :: forcea(:,:) ! focea(0:3,nat) : which atoms can have nonzero forces because they are at non-symmetric positions
CONTAINS
  SUBROUTINE init_forces(nat)
    IMPLICIT NONE
    INTEGER, intent(in) :: nat
    allocate( forcea(0:3,nat) )
    forcea(:,:)=.False.
  END SUBROUTINE init_forces
  SUBROUTINE fini_forces()
    deallocate( forcea )
  END SUBROUTINE fini_forces

  Function Force1_DMFT(nemin,nemax,DM_nemaxx,lmx,alm,blm,clm,aalm,bblm,cclm,weight,el,elo,wEpsw,ri_mat,ilo,lmax2,nume,nloat,lomax,nbands_dft,nbands) result(fsph)
    IMPLICIT NONE
    REAL*8 :: fsph(3)
    INTEGER, intent(in)    :: nemin, nemax, lmx, lmax2, nume, nloat,lomax, DM_nemaxx, nbands_dft, nbands
    COMPLEX*16, intent(in) :: alm((lmax2+1)*(lmax2+1),nume), blm((lmax2+1)*(lmax2+1),nume), clm((lomax+1)*(lomax+1),nume,nloat)
    COMPLEX*16, intent(in) :: aalm((lmax2+1)*(lmax2+1),nbands_dft,3), bblm((lmax2+1)*(lmax2+1),nbands_dft,3), cclm((lomax+1)*(lomax+1),nbands_dft,nloat,3)
    COMPLEX*16, intent(in) :: wEpsw(nbands_dft,nbands_dft)
    REAL*8, intent(in)     :: weight(nume)
    REAL*8, intent(in)     :: el(0:lmax2), elo(0:lomax,1:nloat), ri_mat(2+nloat,2+nloat,0:lmax2)
    INTEGER, intent(in)    :: ilo(0:lmax2)
    ! local variables
    REAL*8     :: tforce(3)
    COMPLEX*16 :: zene, afac, bfac, cfac(nloat), csum
    INTEGER    :: num, nim, l, m, ly, jlo, jlop, i, j
    COMPLEX*16, allocatable :: aofac(:), bofac(:), cofac(:,:), tmp(:)
    REAL*8, allocatable     :: sqt(:)
    ! Implementation of Eq.~305 in the notes: (K-K')<chi_K'|(-nabla^2+V_KS)*delta_{ij} w_i - wepsw_{ij} |chi_K> A*_{iK'} A_{iK}

    allocate( aofac(nbands_dft), bofac(nbands_dft), cofac(nbands_dft,nloat), tmp(nbands_dft) )
    allocate( sqt(nbands_dft) )

    fsph(:)=0.d0
    DO l=0,lmx
       DO m=-l,l
          ly=l*(l+1)+m+1
          aofac=0.d0 ; bofac=0.d0 ; cofac=0.d0

          !$OMP PARALLEL DO PRIVATE(nim,jlo,csum,jlop) SHARED(aofac,bofac,cofac) SCHEDULE(STATIC)
          DO num=nemin,DM_nemaxx
             nim=num-nemin+1
             aofac(nim) = alm(ly,num)
             bofac(nim) = blm(ly,num) * ri_mat(2,2,l)
             DO jlo=1,ilo(l)
                aofac(nim) = aofac(nim) + clm(ly,num,jlo) * ri_mat(1,2+jlo,l)
                bofac(nim) = bofac(nim) + clm(ly,num,jlo) * ri_mat(2,2+jlo,l) 
                csum=0.d0
                DO jlop=1,ilo(l)
                   csum = csum + clm(ly,num,jlop) * ri_mat(2+jlo,2+jlop,l)
                ENDDO
                cofac(nim,jlo) =  alm(ly,num) * ri_mat(1,2+jlo,l) + blm(ly,num) * ri_mat(2,2+jlo,l) + csum
             ENDDO
          ENDDO
          !$OMP END PARALLEL DO
          
          tmp(:) = matmul(conjg(wEpsw), aofac )
          aofac(:) = tmp(:)
          tmp(:) = matmul(conjg(wEpsw), bofac )
          bofac(:) = tmp(:)
          do jlo=1,ilo(l)
             tmp(:) = matmul(conjg(wEpsw), cofac(:,jlo) )
             cofac(:,jlo) = tmp(:)
          enddo
          
          tforce(:)=0.d0
          !$OMP PARALLEL DO PRIVATE(nim,afac,bfac,jlo,zene,csum,jlop,cfac) SCHEDULE(STATIC) REDUCTION(+:tforce)
          DO num=nemin,DM_nemaxx
             nim=num-nemin+1
             afac = 2.d0*alm(ly,num)*el(l) + blm(ly,num)
             bfac = alm(ly,num) + 2.d0*blm(ly,num)*el(l) * ri_mat(2,2,l)
             DO jlo=1,ilo(l)
                zene = elo(l,jlo)+el(l)
                afac = afac + clm(ly,num,jlo) * zene * ri_mat(1,2+jlo,l)
                bfac = bfac + clm(ly,num,jlo) * ( ri_mat(1,2+jlo,l) + zene * ri_mat(2,2+jlo,l) )
                csum=0.d0
                DO jlop=1,ilo(l)
                   csum = csum + clm(ly,num,jlop) * ri_mat(2+jlo,2+jlop,l) * (elo(l,jlo) + elo(l,jlop))
                ENDDO
                cfac(jlo) = alm(ly,num) * zene * ri_mat(1,2+jlo,l) + blm(ly,num) * ri_mat(1,2+jlo,l) + blm(ly,num) * zene * ri_mat(2,2+jlo,l) + csum
             ENDDO
             tforce(:) = tforce(:) + aimag( aalm(ly,nim,:) * dconjg( afac * weight(num) - 2.d0 * aofac(nim) ) )
             tforce(:) = tforce(:) + aimag( bblm(ly,nim,:) * dconjg( bfac * weight(num) - 2.d0 * bofac(nim) ) )
             do jlo=1,ilo(l)
                tforce(:) = tforce(:) + aimag( cclm(ly,nim,jlo,:) * dconjg( cfac(jlo) * weight(num) - 2.d0 * cofac(nim,jlo) ) )
             enddo
          ENDDO
          !$OMP END PARALLEL DO
          fsph(:) = fsph(:) + tforce(:)
       ENDDO
    ENDDO
    deallocate( aofac, bofac, cofac, tmp )
    deallocate( sqt )
    !print *, 'Force1=', fsph(:)
  END Function Force1_DMFT

  Function Force1(nemin,nemax,DM_nemaxx,lmx,alm,blm,clm,aalm,bblm,cclm,weight,el,elo,E,ri_mat,ilo,lmax2,nume,nloat,lomax,nbands_dft) result(fsph)
    IMPLICIT NONE
    REAL*8 :: fsph(3)
    INTEGER, intent(in)    :: nemin, nemax, lmx, lmax2, nume, nloat,lomax, DM_nemaxx, nbands_dft
    COMPLEX*16, intent(in) :: alm((lmax2+1)*(lmax2+1),nume), blm((lmax2+1)*(lmax2+1),nume), clm((lomax+1)*(lomax+1),nume,nloat)
    COMPLEX*16, intent(in) :: aalm((lmax2+1)*(lmax2+1),nbands_dft,3), bblm((lmax2+1)*(lmax2+1),nbands_dft,3), cclm((lomax+1)*(lomax+1),nbands_dft,nloat,3)
    REAL*8, intent(in)     :: weight(nume)
    REAL*8, intent(in)     :: el(0:lmax2), elo(0:lomax,1:nloat), E(nume), ri_mat(2+nloat,2+nloat,0:lmax2)
    INTEGER, intent(in)    :: ilo(0:lmax2)
    ! local variables
    REAL*8     :: tforce(3)
    COMPLEX*16 :: zene, afac, bfac, cfac, cfacf(nloat)
    INTEGER    :: num, nim, l, m, ly, jlo, jlop
    ! Implementation of Eq.~68 in the notes: (K-K')<chi_K'|-nabla^2+V_KS-epsi|chi_K> A*_{iK'} A_{iK}
    fsph(:)=0.d0
    DO num=nemin,DM_nemaxx !nemax
       nim=num-nemin+1
       tforce(:)=0.d0
       DO l=0,lmx
          DO m=-l,l
             ly=l*(l+1)+m+1
             afac = 2.d0*alm(ly,num)*(el(l)-E(num)) + blm(ly,num)
             bfac = alm(ly,num) + 2.d0*blm(ly,num)*(el(l)-E(num)) * ri_mat(2,2,l)
             cfacf = (0.0d0,0.0d0)
             DO jlo=1,ilo(l)
                zene = elo(l,jlo)+el(l)-2.d0*E(num)
                afac = afac + clm(ly,num,jlo) * zene * ri_mat(1,2+jlo,l)
                bfac = bfac + clm(ly,num,jlo) * ( ri_mat(1,2+jlo,l) + zene*ri_mat(2,2+jlo,l) )
                cfacf(jlo) =  ( blm(ly,num) + alm(ly,num) * zene ) * ri_mat(1,2+jlo,l) + blm(ly,num) * zene * ri_mat(2,2+jlo,l)
             ENDDO
             tforce(:) = tforce(:) + dimag( aalm(ly,nim,:)*dconjg(afac) + bblm(ly,nim,:)*dconjg(bfac))
             DO jlo=1,ilo(l)
                DO jlop=1,ilo(l)
                   cfac = 2.d0*clm(ly,num,jlo)*(elo(l,jlo)-E(num))*ri_mat(2+jlo,2+jlop,l)
                   tforce(:) = tforce(:) + dimag( cclm(ly,nim,jlo,:)*dconjg(cfacf(jlop)+cfac) )
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       fsph(:) = fsph(:) + weight(num)*tforce(:)
    ENDDO
  END Function Force1
  
  Function Force1_corrected(nemin,nemax,DM_nemaxx,lmx,alm,blm,clm,aalm,bblm,cclm,weight,el,elo,E,ri_mat,ilo,lmax2,nume,nloat,lomax,nbands_dft) result(fsph)
    IMPLICIT NONE
    REAL*8 :: fsph(3)
    INTEGER, intent(in)    :: nemin, nemax, lmx, lmax2, nume, nloat,lomax, DM_nemaxx, nbands_dft
    COMPLEX*16, intent(in) :: alm((lmax2+1)*(lmax2+1),nume), blm((lmax2+1)*(lmax2+1),nume), clm((lomax+1)*(lomax+1),nume,nloat)
    COMPLEX*16, intent(in) :: aalm((lmax2+1)*(lmax2+1),nbands_dft,3), bblm((lmax2+1)*(lmax2+1),nbands_dft,3), cclm((lomax+1)*(lomax+1),nbands_dft,nloat,3)
    REAL*8, intent(in)     :: weight(nume)
    REAL*8, intent(in)     :: el(0:lmax2), elo(0:lomax,1:nloat), E(nume), ri_mat(2+nloat,2+nloat,0:lmax2)
    INTEGER, intent(in)    :: ilo(0:lmax2)
    ! local variables
    REAL*8     :: tforce(3)
    COMPLEX*16 :: zene, afac, bfac, cfac, cfacf(nloat)
    INTEGER    :: num, nim, l, m, ly, jlo, jlop
    ! Implementation of Eq.~68 in the notes: (K-K')<chi_K'|-nabla^2+V_KS-epsi|chi_K> A*_{iK'} A_{iK}
    fsph(:)=0.d0
    DO num=nemin,DM_nemaxx !nemax
       nim=num-nemin+1
       tforce(:)=0.d0
       DO l=0,lmx
          DO m=-l,l
             ly=l*(l+1)+m+1
             afac = 2.d0*alm(ly,num)*(el(l)-E(num)) + blm(ly,num)
             bfac = alm(ly,num) + 2.d0*blm(ly,num)*(el(l)-E(num)) * ri_mat(2,2,l)
             cfacf = (0.0d0,0.0d0)
             DO jlo=1,ilo(l)
                zene = elo(l,jlo)+el(l)-2.d0*E(num)
                afac = afac + clm(ly,num,jlo) * zene * ri_mat(1,2+jlo,l)
                bfac = bfac + clm(ly,num,jlo) * ( ri_mat(1,2+jlo,l) + zene*ri_mat(2,2+jlo,l) )
                cfacf(jlo) =  ( blm(ly,num) + alm(ly,num) * zene ) * ri_mat(1,2+jlo,l) + blm(ly,num) * zene * ri_mat(2,2+jlo,l)
             ENDDO
             tforce(:) = tforce(:) + dimag( aalm(ly,nim,:)*dconjg(afac) + bblm(ly,nim,:)*dconjg(bfac))
             DO jlo=1,ilo(l)
                DO jlop=1,ilo(l)
                   cfac = 2.d0*clm(ly,num,jlop)*(elo(l,jlop)-E(num))*ri_mat(2+jlo,2+jlop,l) ! This line is corrected from Wien2K. I believe it was a bug in Wien2k.
                   tforce(:) = tforce(:) + dimag( cclm(ly,nim,jlo,:)*dconjg(cfac) )         ! This line is corrected from Wien2K. I believe it was a bug in Wien2k.
                ENDDO
                tforce(:) = tforce(:) + dimag( cclm(ly,nim,jlo,:)*dconjg(cfacf(jlo)) )      ! This line is corrected from Wien2K. I believe it was a bug in Wien2k.
             ENDDO
          ENDDO
       ENDDO
       fsph(:) = fsph(:) + weight(num)*tforce(:)
    ENDDO
  END Function Force1_Corrected

  Function Force2(nemin,nemax,DM_nemaxx,Nhmx,alm,blm,clm,aalm,bblm,cclm,weight,Vts,lv,lpv,mv,mpv,ilo,lmax2,lomax,nume,nloat,ngau,nbands_dft) result(fnsp)
    IMPLICIT NONE
    REAL*8 :: fnsp(3)
    INTEGER, intent(in)    :: nemin, nemax, Nhmx, DM_nemaxx, nbands_dft, lomax
    COMPLEX*16, intent(in) :: alm((lmax2+1)*(lmax2+1),nume), blm((lmax2+1)*(lmax2+1),nume), clm((lomax+1)*(lomax+1),nume,nloat)
    COMPLEX*16, intent(in) :: aalm((lmax2+1)*(lmax2+1),nbands_dft,3), bblm((lmax2+1)*(lmax2+1),nbands_dft,3), cclm((lomax+1)*(lomax+1),nbands_dft,nloat,3)
    COMPLEX*16, intent(in) :: Vts(3,3,ngau)
    INTEGER, intent(in)    :: lv(ngau),lpv(ngau),mv(ngau),mpv(ngau)
    INTEGER, intent(in)    :: ilo(0:lmax2),lmax2,nume,nloat,ngau
    REAL*8, intent(in)     :: weight(nume)
    ! local variables
    REAL*8     :: tforce(3)
    COMPLEX*16 :: afac(3)
    INTEGER    :: ih, iq, jlo, jlop, lpy, ly, num, nim
    ! Implementation of Eq.~69 in the notes: (K-K') A*_{iK'} A_{iK} * <chi_K'|V_{non-sph}|chi_K>
    !   it also goes under Eq.~305 non-symmetric part
    ! tsp[:,:]=[[ tuu,   tud,  tuu12],
    !           [ tdu,   tdd,  tdu12],
    !           [ tuu21, tud21,tuu22]] 
    fnsp(:) = 0.d0
    DO num=nemin,DM_nemaxx! nemax
       nim=num-nemin+1
       tforce(:)=0.d0
       DO ih=1,Nhmx
          lpy= lpv(ih)*(lpv(ih)+1)+mpv(ih)+1
          ly = lv(ih) *( lv(ih)+1)+ mv(ih)+1
          do iq=1,3
             afac(iq) = dconjg(alm(lpy,num)) * Vts(1,iq,ih) + dconjg(blm(lpy,num))* Vts(2,iq,ih)
             DO jlo=1,ilo(lpv(ih))
                afac(iq) = afac(iq) + dconjg(clm(lpy,num,jlo))* Vts(3,iq,ih)
             ENDDO
          enddo
          tforce(:) = tforce(:) + 2.d0*aimag( afac(1) * aalm(ly,nim,:) + afac(2)*bblm(ly,nim,:))
          DO jlop=1,ilo(lv(ih))
             tforce(:) = tforce(:) + 2.d0*aimag( afac(3) * cclm(ly,nim,jlop,:) )
          ENDDO
       ENDDO
       fnsp(:) = fnsp(:) + weight(num)*tforce(:)
    ENDDO
    !print *, 'Force2=', fnsp(:)
  END Function  Force2
  
  Function Force3(nemin,nemax,DM_nemaxx,lmx,alm,blm,clm,aalm,bblm,cclm,weight,ul_Rmt,dul_Rmt,Rmt,ilo,lmax2,lomax,nume,nloat,nbands_dft) result(fsph2)
    IMPLICIT NONE
    REAL*8 :: fsph2(3)
    INTEGER, intent(in)    :: nemin, nemax, lmx, lmax2, nume, nloat, DM_nemaxx, nbands_dft, lomax
    COMPLEX*16, intent(in) :: alm((lmax2+1)*(lmax2+1),nume), blm((lmax2+1)*(lmax2+1),nume), clm((lomax+1)*(lomax+1),nume,nloat)
    COMPLEX*16, intent(in) :: aalm((lmax2+1)*(lmax2+1),nbands_dft,3), bblm((lmax2+1)*(lmax2+1),nbands_dft,3), cclm((lomax+1)*(lomax+1),nbands_dft,nloat,3)
    REAL*8, intent(in)     :: ul_Rmt(nloat+2,0:lmax2), dul_Rmt(nloat+2,0:lmax2), Rmt
    INTEGER, intent(in)    :: ilo(0:lmax2)
    REAL*8, intent(in)     :: weight(nume)
    ! local variables
    REAL*8     :: tforce(3)
    COMPLEX*16 :: kinfac1,kinfac2,kinfac3(3),kinfac4(3)
    INTEGER    :: jlo, l, ly, m, num, nim
    ! Implementation of Eq.~70 in the notes: (K-K')A*_{iK'} A_{iK} * int dS (chi_{K'} \nabla chi_{K})
    !  It also goes under the name Eq.~307
    fsph2(:)=0.d0
    DO num=nemin,DM_nemaxx ! nemax
       nim=num-nemin+1
       tforce(:)=0.d0
       DO l=0,lmx
          DO m=-l,l
             ly=l*(l+1)+m+1
             kinfac1    = alm(ly,num)    *  ul_Rmt(1,l)  + blm(ly,num)    *  ul_Rmt(2,l)
             kinfac2    = alm(ly,num)    * dul_Rmt(1,l)  + blm(ly,num)    * dul_Rmt(2,l)
             kinfac3(:) = aalm(ly,nim,:) * dul_Rmt(1,l)  + bblm(ly,nim,:) * dul_Rmt(2,l)
             kinfac4(:) = aalm(ly,nim,:) *  ul_Rmt(1,l)  + bblm(ly,nim,:) *  ul_Rmt(2,l)
             DO jlo=1,ilo(l)
                kinfac1    = kinfac1    + clm(ly,num,jlo)    *  ul_Rmt(2+jlo,l)
                kinfac2    = kinfac2    + clm(ly,num,jlo)    * dul_Rmt(2+jlo,l)
                kinfac3(:) = kinfac3(:) + cclm(ly,nim,jlo,:) * dul_Rmt(2+jlo,l)
                kinfac4(:) = kinfac4(:) + cclm(ly,nim,jlo,:) *  ul_Rmt(2+jlo,l)
             ENDDO
             tforce(:) = tforce(:) + Rmt**2 * aimag(dconjg(kinfac1)*kinfac3(:)+dconjg(kinfac2)*kinfac4(:))
          ENDDO
       ENDDO
       fsph2(:) = fsph2(:) + weight(num)*tforce(:)
    ENDDO
    !print *, 'Force3=', fsph2(:)
  END Function Force3


  Function Force5_mine(nemin,nemax,lmx,alm,blm,clm,weight,ul_Rmt,dul_Rmt,Rmt,ilo,lmax2,lomax,nume,nloat) result(f71)
    IMPLICIT NONE
    REAL*8     :: f71(3)
    INTEGER, intent(in) :: nemin, nemax, lmx, lmax2, nume, nloat, lomax
    COMPLEX*16, intent(in) :: alm((lmax2+1)*(lmax2+1),nume), blm((lmax2+1)*(lmax2+1),nume), clm((lomax+1)*(lomax+1),nume,nloat)
    REAL*8, intent(in)     :: ul_Rmt(nloat+2,0:lmax2), dul_Rmt(nloat+2,0:lmax2), Rmt
    INTEGER, intent(in)    :: ilo(0:lmax2)
    REAL*8, intent(in)     :: weight(nume)
    !
    COMPLEX*16 :: YTY(3), tforce(3)
    REAL*8     :: cs(2,0:lmx), ds(2,0:lmx)
    COMPLEX*16 :: kinfac(2,(lmx+1)*(lmx+1) )
    COMPLEX*16, parameter, dimension(3) :: vmpm = [ ( 1.d0,0.d0), (0.d0,-1.d0), (0.d0,0.d0)]
    COMPLEX*16, parameter, dimension(3) :: vmmm = [ (-1.d0,0.d0), (0.d0,-1.d0), (0.d0,0.d0)]
    COMPLEX*16, parameter, dimension(3) :: vmem = [ ( 0.d0,0.d0), (0.d0, 0.d0), (0.d0,1.d0)]
    REAL*8     :: VFA, VFF
    INTEGER    :: l,m,lp,mp, ly, lyp, jlo, num, typ
    !*******************Statement Functions***************************
    VFA(L,M) = SQRT((L+M+1.)*(L+M+2.)/((2.*L+3.)*(2.*L+1.)))
    VFF(L,M) = SQRT((L+M+1.)*(L-M+1.)/((2.*L+3.)*(2.*L+1.)))
    !**************************************************************
    do l=0,lmx
       cs(1,l) = 0.5d0          ! Here T = e_r
       ds(1,l) = 0.5d0
       cs(2,l) = 0.5d0*l*(l+2.) ! Here T = \nabla . \nabla e_r
       ds(2,l) = 0.5d0*(l-1.)*(l+1.)
    enddo
    
    DO num=nemin,nemax
       DO l=0,lmx
          DO m=-l,l
             ly=l*(l+1)+m+1
             ! kinfac1 <-  a_{i,lm} S * du_{lm}(r)/dr|_S
             ! kinfac2 <-  a_{i,lm} u_{lm}(S)
             kinfac(1,ly) = (alm(ly,num) * dul_Rmt(1,l)  + blm(ly,num) * dul_Rmt(2,l)) * Rmt
             kinfac(2,ly) =  alm(ly,num) *  ul_Rmt(1,l)  + blm(ly,num) *  ul_Rmt(2,l)
             DO jlo=1,ilo(l)
                kinfac(1,ly) = kinfac(1,ly) + clm(ly,num,jlo) * dul_Rmt(2+jlo,l) * Rmt
                kinfac(2,ly) = kinfac(2,ly) + clm(ly,num,jlo) *  ul_Rmt(2+jlo,l)
             ENDDO
          ENDDO
       ENDDO

       ! Using analytical expressions for Int[Y*_{l'm'} T Y_{lm} dOmega]
       ! typ==1 <-  (a*_{i,l'm'} * S*du_{l'm'}/dr) * ( a_{i,lm} S*du_{lm}/dr ) * (nabla Y*_{l'm'})(nabla Y_{lm}) e_r 
       ! typ==2 <-  (a*_{i,l'm'} * u_{l'm'}(S) ) * ( a_{i,lm} u_{lm}(S) ) * (nabla Y*_{l'm'})(nabla Y_{lm}) e_r 
       tforce(:)=0.d0
       do l=0,lmx
          do m=-l,l
             ly=l*(l+1)+m+1
             do typ=1,2
                lp=l+1
                mp=m+1
                if (lp.le.lmx) then
                   lyp=lp*(lp+1)+mp+1
                   YTY(:) = cs(typ,l)*vfa(l,m)*vmpm(:)
                   tforce(:) = tforce(:) + dconjg(kinfac(typ, lyp)) * kinfac(typ,ly) * YTY
                endif
                mp=m-1
                if (lp.le.lmx .and. mp.ge.-lmx) then
                   lyp=lp*(lp+1)+mp+1
                   YTY(:) = cs(typ,l)*vfa(l,-m)*vmmm(:)
                   tforce(:) = tforce(:) + dconjg(kinfac(typ, lyp)) * kinfac(typ,ly) * YTY
                endif
                mp=m
                if (lp.le.lmx) then
                   lyp=lp*(lp+1)+mp+1
                   YTY(:) = 2*cs(typ,l)*vff(l,m)*vmem(:)
                   tforce(:) = tforce(:) + dconjg(kinfac(typ, lyp)) * kinfac(typ,ly) * YTY
                endif
                
                lp=l-1
                mp=m+1
                if (lp.ge.0 .and. mp.le.lp) then
                   lyp=lp*(lp+1)+mp+1
                   YTY(:) = -ds(typ,l)*vfa(lp,-mp)*vmpm(:)
                   tforce(:) = tforce(:) + dconjg(kinfac(typ, lyp)) * kinfac(typ,ly) * YTY
                endif
                mp=m-1
                if (lp.ge.0 .and. mp.ge.-lmx) then
                   lyp=lp*(lp+1)+mp+1
                   YTY(:) = -ds(typ,l)*vfa(lp,mp)*vmmm(:)
                   tforce(:) = tforce(:) + dconjg(kinfac(typ, lyp)) * kinfac(typ,ly) * YTY
                endif
                mp=m
                if (lp.ge.0 .and. mp.le.lp) then
                   lyp=lp*(lp+1)+mp+1
                   YTY(:) = 2*ds(typ,l)*vff(lp,mp)*vmem(:)
                   tforce(:) = tforce(:) + dconjg(kinfac(typ, lyp)) * kinfac(typ,ly) * YTY
                endif
             end do
          end do
       end do
       f71(:) = f71(:) - weight(num)*dreal(tforce(:))
    ENDDO
  End Function Force5_mine
  
  Function Force4_mine(Clm_new,Vlm,R,lmmax,lm,jri,dx,nrad,lmax2,ncom) result(fvdrho)
    ! Implementation of Eq.~72 in the notes: Int( V_{KS}(r) * nabla rho(r) )
    ! ************ Calculates the following three dimensional integral ****************
    !     Int[  V_{ks} * (nabla)rho ]
    !
    !     Both V_{ks} and rho are given in terms of spherical harmonics expansion
    !     The integral over angle is performed analytically, while the radial integral
    !     is performed numerically.
    !
    !     This implements Eq.72 which also appears as Eq.~308 in the notes.
    ! 
    !    For derivation of the formula, see: Computer Physics Communications 94, 31-48 (1996).
    ! *********************************************************************************
    IMPLICIT NONE
    REAL*8 :: fvdrho(3)
    REAL*8, intent(in)  :: Clm_new(nrad,0:ncom-1)  ! density rho(r)
    REAL*8, intent(in)  :: Vlm(nrad,0:ncom-1)     ! potential V(r)
    REAL*8, intent(in)  :: R(nrad), dx
    INTEGER, intent(in) :: lm(2,0:ncom-1), lmax2, nrad, ncom, lmmax, jri
    ! functions
    Interface
       real*8 FUNCTION int_radial(funct,RX,DELTA,nr)
         INTEGER, intent(in):: nr
         REAL*8, intent(in) :: RX(nr), funct(nr), DELTA
       end FUNCTION int_radial
    end Interface
    ! locals
    REAL*8 :: intfeld(2,0:lmax2,0:lmax2,-1:1,0:lmax2,0:lmax2,-1:1)
    REAL*8  :: toint(NRAD),derv(NRAD)
    INTEGER :: l1, m1, l2, m2, ik1, ik2, lm1, lm2, llmax, L, M, typ, lp, mp
    REAL*8 :: VFA, VFF, forvrx, forvry, forvrz, overal
    REAL*8, allocatable :: cs(:,:), ds(:,:)
    REAL*8 :: dlts(0:1), yTy(-1:1,-1:1)
    !*******************Statement Functions***************************
    VFA(L,M) = SQRT((L+M+1.)*(L+M+2.)/((2.*L+3.)*(2.*L+1.)))
    VFF(L,M) = SQRT((L+M+1.)*(L-M+1.)/((2.*L+3.)*(2.*L+1.)))
    !**************************************************************
    overal=-1  ! W2k uses non-standard convention for Ylm's, namely Ylm's from classical mechanics instead of more standard
               ! quantum convention, which has extra (-1)^m . As Kurki-Sunio uses standard quantum definition of Ylm's, there is an extra minus sign in xy components.
    
    llmax=iabs(lm(1,lmmax-1))
    allocate( cs(2,0:llmax), ds(2,0:llmax) )
    do l=0,llmax
       cs(1,l) = 0.5
       ds(1,l) = 0.5
       cs(2,l) = -0.5*l
       ds(2,l) = 0.5*(l+1)
    enddo
    !***************** Integrals of V*nabla*rho ********  
    forvrz=0.D0
    forvry=0.D0
    forvrx=0.D0
    intfeld(1:2,0:lmax2,0:lmax2,-1:1,0:lmax2,0:lmax2,-1:1)=0.D0   !1: Int( V_{l1,m1}(r) * drho_{l2,m2}(r)/dr r**2 dr)
                                                                  !2: Int( V_{l1,m1}(r) * rho_{l2,m2}(r) / r r**2 dr)
  
    DO lm1=0,lmmax-2 
       l1=abs(lm(1,lm1)) 
       m1=lm(2,lm1)
       ik1 = sign(1, lm(1,lm1) )
       DO lm2=lm1+1,lmmax-1
          l2 = abs(lm(1,lm2))
          m2 = lm(2,lm2)
          ik2 = sign(1, lm(1,lm2) )
          if ( (l2 .eq. (l1+1)) .and. ( (m2.eq.m1) .or. (m2.eq.(m1-1)) .or. (m2.eq.(m1+1)) )) then
             ! Computes intfeld2(l1,m1,l2,m2) = Integrate( V_{lm1}(r)*rho_{lm2}(r)/r * r**2 dr)
             !    Note that rho(r)*r**2 = Clm(r)
             !    intfeld2(l1,m1,l2,m2) = Int( V_{lm}(r) * (r**2*rho(r)_{l+1,m+(0,1,-1)})/r )
             toint(:jri) = Vlm(:jri,lm1) * Clm_new(:jri,lm2) / R(:jri)
             intfeld(2,l1,m1,ik1,l2,m2,ik2) = int_radial(toint(:jri),r(:jri),dx,jri)
             !    Here we just replace lm1 <-> lm2 so that the loop is only over l2>l1
             toint(:jri) = Vlm(:jri,lm2) * Clm_new(:jri,lm1) / R(:jri)
             intfeld(2,l2,m2,ik2,l1,m1,ik1) = int_radial(toint(:jri),r(:jri),dx,jri)
             !*************************************************************************
             ! Computes intfeld1(l1,m1,l2,m2) = Integrate( V_{lm1}(r)* d/dr rho_{lm2}(r) * r**2 dr)
             !    Note that rho(r)*r**2 = Clm(r) hence  d/dr ( Clm/r**2) = ( dClm/dr - 2 Clm/r )/r**2 and
             !    intfeld1(l1,m1,l2,m2) = Integrate( V_{lm1}(r)* ( dClm/dr - 2 Clm/r ) dr)
             call dfrad(derv(:jri),Clm_new(:jri,lm2),r(:jri),jri)
             toint(:jri) = Vlm(:jri,lm1) * ( derv(:jri) - 2.D0*Clm_new(:jri,lm2)/R(:jri) )
             intfeld(1,l1,m1,ik1,l2,m2,ik2) = int_radial(toint(:jri),r(:jri),dx,jri)
             !    Here we just replace lm1 <-> lm2 so that the loop is only over l2>l1 
             call dfrad(derv(:jri),Clm_new(:jri,lm1),r(:jri),jri)
             toint(:jri) = Vlm(:jri,lm2) * ( derv(:jri) - 2.D0*Clm_new(:jri,lm1)/R(:jri) )
             intfeld(1,l2,m2,ik2,l1,m1,ik1) = int_radial(toint(:jri),r(:jri),dx,jri)
             !*************************************************************************
          endif
       ENDDO
    ENDDO
    DO l=0,llmax
       DO m=0,l
          ! XY - components first
          mp = m+1
          dlts(:)=1.d0
          if (m.eq.0) then
             dlts(0) = sqrt(2.d0) ! because cos(0*phi)=1
             dlts(1) = 0.d0       ! because sin(0*phi)=0
          endif
          dlts(:) = dlts(:)*overal
          
          lp = l+1
          if (lp.le.llmax) then
             do typ=1,2
                yTy(+1,+1) = -cs(typ,l)*vfa(l,m)*dlts(0)   ! <y_{l'm',+}|T|y_{lm+}>
                yTy(-1,-1) = -cs(typ,l)*vfa(l,m)*dlts(1)   ! <y_{l'm',-}|T|y_{lm-}>
                yTy(+1,-1) =  cs(typ,l)*vfa(l,m)*dlts(1)   ! <y_{l'm',+}|T|y_{lm-}>
                yTy(-1,+1) = -cs(typ,l)*vfa(l,m)*dlts(0)   ! <y_{l'm',-}|T|y_{lm+}>
                forvrx = forvrx + intfeld(typ,lp,mp,+1,l,m,+1)*yTy(+1,+1)
                forvrx = forvrx + intfeld(typ,lp,mp,-1,l,m,-1)*yTy(-1,-1)
                forvry = forvry + intfeld(typ,lp,mp,+1,l,m,-1)*yTy(+1,-1)
                forvry = forvry + intfeld(typ,lp,mp,-1,l,m,+1)*yTy(-1,+1)
             enddo
          endif
          lp=l-1
          if (lp.ge.0 .and. mp.le.lp) then
             do typ=1,2
                yTy(+1,+1) =  ds(typ,l)*vfa(lp,-mp)*dlts(0)
                yTy(-1,-1) =  ds(typ,l)*vfa(lp,-mp)*dlts(1)
                yTy(+1,-1) = -ds(typ,l)*vfa(lp,-mp)*dlts(1)
                yTy(-1,+1) =  ds(typ,l)*vfa(lp,-mp)*dlts(0)
                forvrx = forvrx + intfeld(typ,lp,mp,+1,l,m,+1)*yTy(+1,+1)
                forvrx = forvrx + intfeld(typ,lp,mp,-1,l,m,-1)*yTy(-1,-1)
                forvry = forvry + intfeld(typ,lp,mp,+1,l,m,-1)*yTy(+1,-1)
                forvry = forvry + intfeld(typ,lp,mp,-1,l,m,+1)*yTy(-1,+1)
             enddo
          endif

          mp = m-1
          dlts(:)=1.d0
          if (mp.eq.0) then
             dlts(0) = sqrt(2.d0)
             dlts(1) = 0.d0
          endif
          dlts(:) = dlts(:)*overal
                
          lp = l+1
          if (mp.ge.0 .and. lp.le.llmax) then
             do typ=1,2
                yTy(+1,+1) = cs(typ,l)*vfa(l,-m)*dlts(0)
                yTy(-1,-1) = cs(typ,l)*vfa(l,-m)*dlts(1)
                yTy(+1,-1) = cs(typ,l)*vfa(l,-m)*dlts(0)
                yTy(-1,+1) =-cs(typ,l)*vfa(l,-m)*dlts(1)
                forvrx = forvrx + intfeld(typ,lp,mp,+1,l,m,+1)*yTy(+1,+1)
                forvrx = forvrx + intfeld(typ,lp,mp,-1,l,m,-1)*yTy(-1,-1)
                forvry = forvry + intfeld(typ,lp,mp,+1,l,m,-1)*yTy(+1,-1)
                forvry = forvry + intfeld(typ,lp,mp,-1,l,m,+1)*yTy(-1,+1)
             enddo
          endif
          lp=l-1
          if (mp.ge.0 .and. lp.ge.0) then
             do typ=1,2
                yTy(+1,+1) = -ds(typ,l)*vfa(lp,mp)*dlts(0)
                yTy(-1,-1) = -ds(typ,l)*vfa(lp,mp)*dlts(1)
                yTy(+1,-1) = -ds(typ,l)*vfa(lp,mp)*dlts(0)
                yTy(-1,+1) =  ds(typ,l)*vfa(lp,mp)*dlts(1)
                forvrx = forvrx + intfeld(typ,lp,mp,+1,l,m,+1)*yTy(+1,+1)
                forvrx = forvrx + intfeld(typ,lp,mp,-1,l,m,-1)*yTy(-1,-1)
                forvry = forvry + intfeld(typ,lp,mp,+1,l,m,-1)*yTy(+1,-1)
                forvry = forvry + intfeld(typ,lp,mp,-1,l,m,+1)*yTy(-1,+1)
             enddo
          endif
          ! Z - components
          mp = m
          dlts(:)=1.d0
          if (m.eq.0) then
             dlts(0)=1.d0
             dlts(1)=0.d0
          endif
          lp=l+1
          if (lp.le.llmax) then
             do typ=1,2
                yTy(+1,+1) = 2*cs(typ,l)*vff(l,m)*dlts(0)
                yTy(-1,-1) = 2*cs(typ,l)*vff(l,m)*dlts(1)
                forvrz = forvrz + intfeld(typ,lp,mp,+1,l,m,+1)*yTy(+1,+1)
                forvrz = forvrz + intfeld(typ,lp,mp,-1,l,m,-1)*yTy(-1,-1)
             enddo
          endif
          lp=l-1
          if (lp.ge.0 .and. mp.le.lp) then
             do typ=1,2
                yTy(+1,+1) = 2*ds(typ,l)*vff(lp,mp)*dlts(0)
                yTy(-1,-1) = 2*ds(typ,l)*vff(lp,mp)*dlts(1)
                forvrz = forvrz + intfeld(typ,lp,mp,+1,l,m,+1)*yTy(+1,+1)
                forvrz = forvrz + intfeld(typ,lp,mp,-1,l,m,-1)*yTy(-1,-1)
             enddo
          endif
       ENDDO
    ENDDO
    deallocate( cs, ds )
    
    fvdrho(1)=forvrx
    fvdrho(2)=forvry
    fvdrho(3)=forvrz
    return
  end function Force4_Mine
  
  Function Force4(Clm_new,Vlm,R,lmmax,lm,jri,dx,nrad,lmax2,ncom) result(fvdrho)
    ! Implementation of Eq.~72 in the notes: Int( V_{KS}(r) * nabla rho(r) )
    ! ************ Calculates the following three dimensional integral ****************
    !     Int[  V_{ks} * (nabla)rho ]
    !
    !     Both V_{ks} and rho are given in terms of spherical harmonics expansion
    !     The integral over angle is performed analytically, while the radial integral
    !     is performed numerically.
    !
    !    For derivation of the formula, see: Computer Physics Communications 94, 31-48 (1996).
    ! *********************************************************************************
    IMPLICIT NONE
    REAL*8 :: fvdrho(3)
    REAL*8, intent(in)  :: Clm_new(nrad,0:ncom-1)  ! density rho(r)
    REAL*8, intent(in)  :: Vlm(nrad,0:ncom-1)     ! potential V(r)
    REAL*8, intent(in)  :: R(nrad), dx
    INTEGER, intent(in) :: lm(2,0:ncom-1), lmax2, nrad, ncom, lmmax, jri
    ! functions
    Interface
       real*8 FUNCTION int_radial(funct,RX,DELTA,nr)
         INTEGER, intent(in):: nr
         REAL*8, intent(in) :: RX(nr), funct(nr), DELTA
       end FUNCTION int_radial
    end Interface
    ! locals
    REAL*8 :: intfeld1(0:lmax2,0:lmax2,-1:1,0:lmax2,0:lmax2,-1:1)
    REAL*8 :: intfeld2(0:lmax2,0:lmax2,-1:1,0:lmax2,0:lmax2,-1:1) 
    REAL*8  :: toint(NRAD),derv(NRAD)
    INTEGER :: l1, m1, l2, m2, mz, ik1, ik2, lm1, lm2, llmax, L, M
    REAL*8 :: VFA, VFF, fakt, faktmo, faktpo, forvrx, forvry, forvrz, sqhalf
    !*******************Statement Functions***************************
    VFA(L,M) = SQRT((L+M+1.)*(L+M+2.)/((2.*L+3.)*(2.*L+1.)))
    VFF(L,M) = SQRT((L+M+1.)*(L-M+1.)/((2.*L+3.)*(2.*L+1.)))
    !**************************************************************
    SQHALF=SQRT(0.5D0) 
    !***************** Integrals of V*nabla*rho ********  
    forvrz=0.D0
    forvry=0.D0
    forvrx=0.D0
    intfeld1(0:lmax2,0:lmax2,-1:1,0:lmax2,0:lmax2,-1:1)=0.D0   ! Int( V_{l1,m1}(r) * drho_{l2,m2}(r)/dr r**2 dr)
    intfeld2(0:lmax2,0:lmax2,-1:1,0:lmax2,0:lmax2,-1:1)=0.D0   ! Int( V_{l1,m1}(r) * rho_{l2,m2}(r) / r r**2 dr)
    
    llmax=iabs(lm(1,lmmax-1))
    DO lm1=0,lmmax-2 
       l1=abs(lm(1,lm1)) 
       m1=lm(2,lm1)
       ik1 = sign(1, lm(1,lm1) )
       DO lm2=lm1+1,lmmax-1
          l2 = abs(lm(1,lm2))
          m2 = lm(2,lm2)
          ik2 = sign(1, lm(1,lm2) )
          if ( (l2 .eq. (l1+1)) .and. ( (m2.eq.m1) .or. (m2.eq.(m1-1)) .or. (m2.eq.(m1+1)) )) then
             ! Computes intfeld2(l1,m1,l2,m2) = Integrate( V_{lm1}(r)*rho_{lm2}(r)/r * r**2 dr)
             !    Note that rho(r)*r**2 = Clm(r)
             !    intfeld2(l1,m1,l2,m2) = Int( V_{lm}(r) * (r**2*rho(r)_{l+1,m+(0,1,-1)})/r )
             toint(:jri) = Vlm(:jri,lm1) * Clm_new(:jri,lm2) / R(:jri)
             intfeld2(l1,m1,ik1,l2,m2,ik2) = int_radial(toint(:jri),r(:jri),dx,jri)
             !    Here we just replace lm1 <-> lm2 so that the loop is only over l2>l1
             toint(:jri) = Vlm(:jri,lm2) * Clm_new(:jri,lm1) / R(:jri)
             intfeld2(l2,m2,ik2,l1,m1,ik1) = int_radial(toint(:jri),r(:jri),dx,jri)
             !*************************************************************************
             ! Computes intfeld1(l1,m1,l2,m2) = Integrate( V_{lm1}(r)* d/dr rho_{lm2}(r) * r**2 dr)
             !    Note that rho(r)*r**2 = Clm(r) hence  d/dr ( Clm/r**2) = ( dClm/dr - 2 Clm/r )/r**2 and
             !    intfeld1(l1,m1,l2,m2) = Integrate( V_{lm1}(r)* ( dClm/dr - 2 Clm/r ) dr)
             call dfrad(derv(:jri),Clm_new(:jri,lm2),r(:jri),jri)
             toint(:jri) = Vlm(:jri,lm1) * ( derv(:jri) - 2.D0*Clm_new(:jri,lm2)/R(:jri) )
             intfeld1(l1,m1,ik1,l2,m2,ik2) = int_radial(toint(:jri),r(:jri),dx,jri)
             !    Here we just replace lm1 <-> lm2 so that the loop is only over l2>l1 
             call dfrad(derv(:jri),Clm_new(:jri,lm1),r(:jri),jri)
             toint(:jri) = Vlm(:jri,lm2) * ( derv(:jri) - 2.D0*Clm_new(:jri,lm1)/R(:jri) )
             intfeld1(l2,m2,ik2,l1,m1,ik1) = int_radial(toint(:jri),r(:jri),dx,jri)
             !*************************************************************************
          endif
       ENDDO
    ENDDO
    DO l=0,llmax-1
       DO mz=-l,l
          m=abs(mz)
          fakt=1.d0
          if (m.eq.0) fakt=2.d0
          forvrz = forvrz + 0.5d0*fakt*vff(l,mz)*( (2+l)*( intfeld2(l,  m,+1,l+1,m,+1) + intfeld2(l,  m,-1,l+1,m,-1)) &
                                                      -l*( intfeld2(l+1,m,+1,l,  m,+1) + intfeld2(l+1,m,-1,l,  m,-1)) &
                                                        +( intfeld1(l,  m,+1,l+1,m,+1) + intfeld1(l,  m,-1,l+1,m,-1)) &
                                                        +( intfeld1(l+1,m,+1,l,  m,+1) + intfeld1(l+1,m,-1,l,  m,-1)))
          
          if (mz.gt.0)    faktpo = -1.0d0
          if (mz.lt.(-1)) faktpo =  1.0d0
          if (mz.eq.0)    faktpo = -sqrt(2.0d0)
          if (mz.eq.(-1)) faktpo =  sqrt(2.0d0)
  
          if (mz.gt.1) faktmo = -1.0d0
          if (mz.lt.0) faktmo =  1.0d0
          if (mz.eq.0) faktmo =  sqrt(2.0d0)
          if (mz.eq.1) faktmo = -sqrt(2.0d0)
          !
          forvrx = forvrx - vfa(l, mz)*0.25D0*faktpo*(  (2+l)*( intfeld2(l,  abs(mz),  +1,l+1,abs(mz+1),+1) + intfeld2(l,  abs(mz),  -1,l+1,abs(mz+1),-1)) &
                                                           -l*( intfeld2(l+1,abs(mz+1),+1,l,  abs(mz),  +1) + intfeld2(l+1,abs(mz+1),-1,l,  abs(mz),  -1)) &
                                                             +( intfeld1(l,  abs(mz),  +1,l+1,abs(mz+1),+1) + intfeld1(l,  abs(mz),  -1,l+1,abs(mz+1),-1)) &
                                                             +( intfeld1(l+1,abs(mz+1),+1,l,  abs(mz),  +1) + intfeld1(l+1,abs(mz+1),-1,l,  abs(mz),  -1)))
          forvrx = forvrx + vfa(l,-mz)*0.25D0*faktmo*(  (2+l)*( intfeld2(l,  abs(mz),  +1,l+1,abs(mz-1),+1) + intfeld2(l,  abs(mz),  -1,l+1,abs(mz-1),-1)) &
                                                           -l*( intfeld2(l+1,abs(mz-1),+1,l,  abs(mz),  +1) + intfeld2(l+1,abs(mz-1),-1,l,  abs(mz),  -1)) &
                                                             +( intfeld1(l,  abs(mz),  +1,l+1,abs(mz-1),+1) + intfeld1(l,  abs(mz),  -1,l+1,abs(mz-1),-1)) &
                                                             +( intfeld1(l+1,abs(mz-1),+1,l,  abs(mz),  +1) + intfeld1(l+1,abs(mz-1),-1,l,  abs(mz),  -1)))
          faktmo = 1.0D0
          if (mz.eq.1 .or. mz.eq.0) faktmo = sqrt(2.0d0)
          faktpo = 1.0d0
          if (mz.eq.0 .or. mz.eq.(-1)) faktpo = sqrt(2.0d0)
          forvry = forvry + vfa(l, mz)*0.25D0*faktpo*( (2+l)*( intfeld2(l,  abs(mz),  +1,l+1,abs(mz+1),-1) - intfeld2(l,  abs(mz),  -1,l+1,abs(mz+1),+1)) &
                                                          +l*( intfeld2(l+1,abs(mz+1),+1,l,  abs(mz),  -1) - intfeld2(l+1,abs(mz+1),-1,l,  abs(mz),  +1)) &
                                                            +( intfeld1(l,  abs(mz),  +1,l+1,abs(mz+1),-1) - intfeld1(l,  abs(mz),  -1,l+1,abs(mz+1),+1)) &
                                                            -( intfeld1(l+1,abs(mz+1),+1,l,  abs(mz),  -1) - intfeld1(l+1,abs(mz+1),-1,l,  abs(mz),  +1)))
          forvry = forvry + vfa(l,-mz)*0.25D0*faktmo*( (2+l)*( intfeld2(l,  abs(mz),  +1,l+1,abs(mz-1),-1) - intfeld2(l,  abs(mz),  -1,l+1,abs(mz-1),+1)) &
                                                          +l*( intfeld2(l+1,abs(mz-1),+1,l,  abs(mz),  -1) - intfeld2(l+1,abs(mz-1),-1,l,  abs(mz),  +1)) &
                                                            +( intfeld1(l,  abs(mz),  +1,l+1,abs(mz-1),-1) - intfeld1(l,  abs(mz),  -1,l+1,abs(mz-1),+1)) &
                                                            -( intfeld1(l+1,abs(mz-1),+1,l,  abs(mz),  -1) - intfeld1(l+1,abs(mz-1),-1,l,  abs(mz),  +1)) )
       ENDDO
    ENDDO
    fvdrho(1)=forvrx
    fvdrho(2)=forvry
    fvdrho(3)=forvrz
    return
  end function Force4
  
  ! former-fsumai1
  SUBROUTINE Force_surface(ekink,kzz,pos,rotij,mult,Rmt,BR1,rotloc,iz,tau,iord,nwave,nat,ndif,Qcomplex)
    IMPLICIT NONE
    complex*16, intent(in) :: ekink(nwave)    ! 
    INTEGER,    intent(in) :: kzz(3,nwave)    ! all plane waves for interstital
    REAL*8,     intent(in) :: pos(3,ndif)     ! positions of all atoms in the unit cell
    REAL*8,     intent(in) :: rotij(3,3,ndif) ! transforms the position of an atom to its corresponding position of an equivalent atom
    INTEGER,    intent(in) :: mult(nat)       ! number of atom of each sort
    REAL*8,     intent(in) :: Rmt(nat), BR1(3,3)
    REAL*8,     intent(in) :: rotloc(3,3,nat) ! local axis transformation from struc file
    INTEGER,    intent(in) :: iz(3,3,iord)    ! matrix of transformations for group operations
    REAL*8,     intent(in) :: tau(3,iord)     ! atoms shifts for group operations
    INTEGER,    intent(in) :: iord            ! Number of group operations    
    INTEGER,    intent(in) :: nwave, nat, ndif
    LOGICAL,    intent(in) :: Qcomplex        ! do we require complex calculation (so or absence of inversin)
    ! locals
    REAL*8,PARAMETER :: pi=3.1415926535897932d0
    COMPLEX*16 :: factor, fstar(3), fphase(iord)
    INTEGER    :: lfirst, jatom, j, jj, ik
    INTEGER    :: Nstar, G_star(3,iord)
    REAL*8     :: fj(0:2), Gi(3), G_rot(3), p_this(3)
    REAL*8     :: Gnorm,frmt,pha,fsur_norm, rotloc_x_BR1(3,3), rotloc_x_BR1_x_rotij(3,3)
    REAL*8     :: fsur(3)
    lfirst=1
    DO jatom=1,nat
       ! skip G=(0,0,0): no contribution to ekin
       p_this(1:3) = pos(1:3,lfirst)*2.d0*pi
       fsur(:)=0.d0
       rotloc_x_BR1 = matmul(ROTLOC(:,:,jatom), BR1)
       rotloc_x_BR1_x_rotij = matmul(rotloc_x_BR1, rotij(:,:,lfirst) )
       ! SHOULD USE $OMP
       DO j=2,nwave
          Gi    = matmul(BR1,kzz(1:3,j))    ! plane wave transformed to cartesian coordinates
          Gnorm = dsqrt(sum(Gi(:)**2))
          CALL SPHBES(2, Rmt(jatom)*Gnorm, fj)   ! spherical bessel of j_l(|G|R_{MT})
          CALL STERN(kzz(1:3,j), Nstar, G_star, fphase, iz, tau, iord, Qcomplex)  ! Generates the star of G : G_star -> G, fphase=sum_{start-member=j}e^{2*pi*i*G*tau_j}, Nstar-> number of star members
          fstar(:)=0.d0
          DO jj=1,Nstar
             !     Calculate phase factor
             pha = dot_product(p_this, G_star(:,jj))      ! e^{2*pi*i*G*R_i}, where R_i is current atom
             !     swaped sin and cos because there's an i in the equation, i.e., i*e^{i*pha}
             factor= dcmplx(-sin(pha),cos(pha))*fphase(jj) ! factor = e^{i*G*R_i}*\sum_{start-members=j} e^{i*G*tau_j}
             !G_rot = matmul(rotij(:,:,lfirst),G_star(:,jj))
             !Gi    = matmul(BR1, G_rot)             ! Since we use rotloc here, we transform to local coordinate axis, hence even plane wave term is given in local axis
             !G_rot = matmul(ROTLOC(:,:,jatom),Gi)   ! Here we use rotloc and not crotloc, because the wave functions are compatible with rotloc in the structure file
             G_rot = matmul(rotloc_x_BR1_x_rotij, G_star(:,jj)) 
             fstar(:) =fstar(:) + factor*G_rot(:)         ! fstar = G*factor
          ENDDO
          ! Here we add a term : G/|G|*j1(|G|S)* e^{i*G*R_i}*{1/n_star}\sum_{start-members=j} i*e^{i*G*tau_j}
          fsur(:)=fsur(:) + dble(ekink(j)*fstar(:))/dble(Nstar)*fj(1)/Gnorm
       ENDDO
       frmt=4.d0*pi*Rmt(jatom)**2
       fsur(:) = fsur(:)*frmt     ! Now we have: 4*pi*S**2 * G/|G|*j1(|G|S)* e^{i*G*R_i}*{1/n_star}\sum_{start-members=j} i*e^{i*G*tau_j}
       fsur_norm = sqrt(sum(fsur**2))
       WRITE(6, '(2i3,a7,4e15.7)') jatom,1,'+ SUR', fsur_norm, (fsur(ik),ik=1,3)
       WRITE(21,79) jatom,jatom, fsur_norm*1000, (fsur(ik)*1000,ik=1,3)
       lfirst=lfirst+mult(jatom)
    ENDDO
79  FORMAT(':FSU',i3.3,':',1x,i3,'.ATOM',4F17.9) 
  END SUBROUTINE Force_surface

  Function Force_surface_extra(rhok,ekink,Vlm,rhoalm,kzz,pos,rotij,Rmt,BR1,rotloc,iz,tau,lm,iord,nwave,lmmax,ncom,jatom) result( fsur2 )
    IMPLICIT NONE
    REAL*8     :: fsur2(9)
    complex*16, intent(in) :: ekink(nwave)
    COMPLEX*16, intent(in) :: rhok(nwave)     ! interstitial charge
    REAL*8,     intent(in) :: Vlm(0:lmmax-1)  ! potential V(R_MT)
    REAL*8,     intent(in) :: rhoalm(0:lmmax-1)! charge at the MT-boundary from APW's
    INTEGER,    intent(in) :: kzz(3,nwave)    ! all plane waves for interstital
    REAL*8,     intent(in) :: pos(3)          ! positions of the first atom of this type
    REAL*8,     intent(in) :: rotij(3,3)      ! transforms the position of an atom to its corresponding position of an equivalent atom
    REAL*8,     intent(in) :: Rmt, BR1(3,3)   ! MT-sphere, transformation to cartesian coordinates
    REAL*8,     intent(in) :: rotloc(3,3)     ! local axis transformation from struc file
    INTEGER,    intent(in) :: iz(3,3,iord)    ! matrix of transformations for group operations
    REAL*8,     intent(in) :: tau(3,iord)     ! atoms shifts for group operations
    INTEGER,    intent(in) :: lm(2,0:ncom-1)
    INTEGER,    intent(in) :: iord            ! Number of group operations    
    INTEGER,    intent(in) :: nwave, lmmax, ncom, jatom
    ! locals
    REAL*8,PARAMETER :: pi=3.1415926535897932d0
    REAL*8, allocatable :: V_lm(:,:,:), Vs(:,:,:,:), jl(:)
    COMPLEX*16, allocatable :: Ylmc(:)
    INTEGER, allocatable    :: lm_ext(:,:)
    COMPLEX*16, allocatable :: rho_part2(:), rho_part3(:), rholm(:)
    COMPLEX*16 :: rho_part, i_to_l
    COMPLEX*16 :: factor, fphase(iord), imag, fstar(3), fkdc(3)
    INTEGER    :: lfirst, i, j, jj, ik, lm1, l, m, s, l_max, lp, mp, m1m, ind_p, ind_m, lmmax_ext
    INTEGER    :: Nstar, G_star(3,iord)
    REAL*8     :: Gi(3), G_rot(3), p_this(3), rotloc_x_BR1(3,3), rotloc_x_BR1_x_rotij(3,3)
    REAL*8     :: Gnorm,frmt,pha,fsur_norm, overal, cs, ds, ylm_pm, nrm, nrm0
    REAL*8     :: fri(3), fra(3), fsur(3)
    REAL*8     :: dlts(0:1), yTy(-1:1,-1:1)
    REAL*8     :: GR
    REAL*8     :: VFA, VFF
    !*******************Statement Functions***************************
    VFA(L,M) = SQRT((L+M+1.)*(L+M+2.)/((2.*L+3.)*(2.*L+1.)))
    VFF(L,M) = SQRT((L+M+1.)*(L-M+1.)/((2.*L+3.)*(2.*L+1.)))
    !*****************************************************************
    imag=dcmplx(0.d0,1.d0)

    !write(6,'(A,F7.4,A)') 'Density at the MT-boundary r=', Rmt, ' from augmented PW:'
    !DO lm1=0,lmmax-1
    !   write(6,'(I3,1x,I3,1x,F20.10)') lm(1,lm1), lm(2,lm1), rhoalm(lm1)
    !END DO
    
    l_max = maxval(abs(lm(1,0:lmmax-1)))

    ! lm will be extended so that the surface integral \oint dS rho*Vks is computed exactly
    lmmax_ext = lmmax + 2*(l_max+1)+1   ! We will add one l shell more
    allocate( lm_ext(2,0:lmmax_ext-1) ) ! lm with one extra shell is lm_ext
    lm_ext(:,0:lmmax-1) = lm(:,0:lmmax-1)
    lm_ext(1,lmmax) = l_max+1
    lm_ext(2,lmmax) = 0
    i=lmmax
    do m=1,l_max+1
       i=i+1
       lm_ext(1,i) = l_max+1
       lm_ext(2,i) = m
       i=i+1
       lm_ext(1,i) = -(l_max+1)
       lm_ext(2,i) = m
    enddo

    fsur(:)=0.d0
    fkdc(:)=0.d0
    !********** Starting density in the interstitials ************!
    rotloc_x_BR1 = matmul(ROTLOC, BR1)
    rotloc_x_BR1_x_rotij = matmul(rotloc_x_BR1, rotij )
    allocate( jl(0:l_max+2) )
    allocate( Ylmc((l_max+2)*(l_max+2)) )
    allocate( rho_part2(0:l_max+1), rho_part3(0:l_max+1) )
    allocate( rholm(0:lmmax_ext-1) )
!!! BRISI
    !allocate( drho_part2(0:l_max), drho_part3(0:l_max) )
    !allocate( drholm(0:lmmax-1) )
    !drholm(:)=0.d0
!!! BRISI
    rholm(:)=0.d0
    p_this(:) = pos(:)*2.d0*pi
    nrm0 = sqrt(2.0d0)    
    DO j=1,nwave
       Gi    = matmul(BR1, kzz(1:3,j) )    ! plane wave transformed to cartesian coordinates
       Gnorm = dsqrt(sum(Gi(:)**2))
       GR = Rmt*Gnorm
       CALL SPHBES(l_max+2, GR, jl)   ! spherical bessel of j_l(|G|R_{MT})
       CALL STERN( kzz(:,j), Nstar, G_star, fphase, iz, tau, iord, .true.)  ! Generates the star of G : G_star -> G, fphase=sum_{start-member=j}e^{2*pi*i*G*tau_j}, Nstar-> number of star members
       rho_part = 4.d0*pi*rhok(j)/dble(Nstar)
       do l=0,l_max+1
          i_to_l  = dcmplx(cos(0.5*pi*l),sin(0.5*pi*l))
          rho_part2(l) = rho_part * i_to_l * jl(l)
          !if (l.le.l_max .and. j.gt.1) then
          !   drho_part2(l) = rho_part * i_to_l * Gnorm * ( l/GR*jl(l) - jl(l+1) )
          !else
          !   drho_part2(l)=0.d0
          !endif
       enddo
       fstar(:)=0.d0
       !$OMP PARALLEL DO PRIVATE(pha,factor,G_rot,Ylmc,rho_part3,lm1,l,m,nrm,ylm_pm) SCHEDULE(STATIC) REDUCTION(+:fstar,rholm)
       DO jj=1,Nstar
          !     Calculate phase factor
          pha    = dot_product(p_this, G_star(:,jj))      ! e^{2*pi*i*G*R_i}, where R_i is current atom writen in the global coordinate system
          factor = dcmplx(cos(pha),sin(pha))*fphase(jj)   ! factor = e^{i*G*R_i}*\sum_{start-members=j} e^{i*G*tau_j}
          G_rot  = matmul(rotloc_x_BR1_x_rotij, G_star(:,jj)) ! Since we use rotloc here, we transform to local coordinate axis. This is needed as Vks is given in local coordinate system.
          
          fstar(:) = fstar(:) + (factor*imag) * G_rot(:)
          
          call YLM(G_rot, l_max+1, Ylmc, (l_max+2)*(l_max+2))
          rho_part3(:) = rho_part2(:) * factor
          !drho_part3(:) = drho_part2(:) * factor
          DO lm1=0,lmmax_ext-1
             l = abs(lm_ext(1,lm1)) 
             m = lm_ext(2,lm1)
             nrm = nrm0*(-1)**m   ! Here we need an extra (-1)**m because w2k uses uses non-standard convention for Ylm's, namely Ylm's from classical mechanics instead of more standard Q.M., which has extra (-1)^m
             if (m.eq.0) nrm=1.d0
             if (lm_ext(1,lm1).ge.0) then  ! s=+1
                ylm_pm =  dble(Ylmc(l**2+l+m+1))*nrm
             else                          ! s=-1
                ylm_pm = aimag(Ylmc(l**2+l+m+1))*nrm
             endif
             rholm(lm1) = rholm(lm1) + rho_part3(l) * ylm_pm
             !if (lm1.le.lmmax-1) then
             !   drholm(lm1) = drholm(lm1) + dble(drho_part3(l)*ylm_pm)
             !endif
          enddo
       ENDDO
       !$OMP END PARALLEL DO
       if (j.gt.1) then
          fsur(:) = fsur(:) + dble(ekink(j)*fstar(:))/dble(Nstar)*jl(1)/Gnorm ! G=(0,0,0) gives no contribution
          fkdc(:) = fkdc(:) + dble(rho_part * fstar(:)) * jl(1) * Gnorm
       endif
    ENDDO
    deallocate( rho_part2, rho_part3 )
    deallocate( jl, Ylmc )
!!! BRISI
    !deallocate( drho_part2, drho_part3 )
!!! BRISI
    fsur(:) = fsur(:)*(4.d0*pi*Rmt**2)
    fkdc(:) = fkdc(:)*(0.5*Rmt**2)
    
    !write(6,'(A,F7.4,A)') 'Density at the MT-boundary r=', Rmt, ' from interstitial PW:'
    !DO lm1=0,lmmax_ext-1
    !   write(6,'(I3,1x,I3,1x,F20.10)') lm_ext(1,lm1), lm_ext(2,lm1), dble(rholm(lm1))
    !END DO
    !********** Finished density in the interstitials ************!

    
!!! BRISI    
    !write(6,'(A,F7.4,A)') 'Radial Derivative of the Density at the MT-boundary r=', Rmt, ' from interstitial PW:'
    !DO lm1=0,lmmax-1
    !   write(6,'(I3,1x,I3,1x,F20.10)') lm(1,lm1), lm(2,lm1), drholm(lm1)
    !END DO
    !deallocate( drholm )
!!! BRISI
    
    overal=-1.d0 ! W2k uses non-standard convention for Ylm's, namely Ylm's from classical mechanics instead of more standard
                 ! quantum convention, which has extra (-1)^m . As Kurki-Sunio uses standard quantum definition of Ylm's, there is an extra minus sign in xy components.
    cs=0.5d0     ! matrix element for <y_lm|e_r|y_lm>
    ds=0.5d0     ! matrix element for <y_lm|e_r|y_lm>
    
    !WRITE(6,'(A,I3,A)') 'Potential at the MT boundary for jatom=', jatom, ':'
    !do lm1=0,lmmax-1
    !   WRITE(6,'(I3,1x,I3,1x,F20.10)') lm(1,lm1), lm(2,lm1), Vlm(lm1)
    !enddo
    
    allocate( V_lm(0:l_max,0:l_max,-1:1) )
    V_lm(:,:,:)=0.d0
    DO lm1=0,lmmax-1
       l = abs(lm(1,lm1)) 
       m = lm(2,lm1)
       s = sign(1, lm(1,lm1) )
       V_lm(l,m,s) = Vlm(lm1)
    ENDDO
    
    !WRITE(6,'(A,I3,A)') 'Potential at the MT boundary for jatom=', jatom, '='
    !do l=0,l_max
    !   do m=0,l
    !      do s=-1,1,2
    !         if (abs(V_lm(l,m,s)).gt.1e-10) WRITE(6,'(I3,1x,I3,1x,I3,1x,F20.10,1x)') l, m, s, V_lm(l,m,s)
    !      enddo
    !   enddo
    !enddo
    !WRITE(6,*)

    !  Vs([x,y,z],[-1,1],m,l)
    allocate( Vs(1:3,-1:1,0:(l_max+1),0:(l_max+1)) ) ! Since Vlm vanishes for l>l_max, Vs will vanish for l>l_max+1
    Vs (:,:,:,:) = 0.d0                          ! This is Vs_{lms} = sum_{l'm's'} V_{l'm's'} * I^1_{l'm's',lms}
    
    DO l=0,l_max+1
       DO m=0,l
          ! XY - components first
          mp = m+1
          dlts(:)=1.d0
          if (m.eq.0) then
             dlts(0) = sqrt(2.d0) ! because cos(0*phi)=1
             dlts(1) = 0.d0       ! because sin(0*phi)=0
          endif
          dlts(:) = dlts(:)*overal

          lp = l+1
          if (lp.le.l_max .and. mp.le.lp) then
             yTy(+1,+1) = -cs*vfa(l,m)*dlts(0)   ! <y_{l'm',+}|T|y_{lm+}>
             yTy(-1,-1) = -cs*vfa(l,m)*dlts(1)   ! <y_{l'm',-}|T|y_{lm-}>
             yTy(+1,-1) =  cs*vfa(l,m)*dlts(1)   ! <y_{l'm',+}|T|y_{lm-}>
             yTy(-1,+1) = -cs*vfa(l,m)*dlts(0)   ! <y_{l'm',-}|T|y_{lm+}>
             Vs(1,+1,m,l) = Vs(1,+1,m,l) + V_lm(lp,mp,+1)*yTy(+1,+1) ! intfeld(typ,lp,mp,+1,l,m,+1)*yTy(+1,+1)
             Vs(1,-1,m,l) = Vs(1,-1,m,l) + V_lm(lp,mp,-1)*yTy(-1,-1) ! intfeld(typ,lp,mp,-1,l,m,-1)*yTy(-1,-1)
             Vs(2,-1,m,l) = Vs(2,-1,m,l) + V_lm(lp,mp,+1)*yTy(+1,-1) ! intfeld(typ,lp,mp,+1,l,m,-1)*yTy(+1,-1)
             Vs(2,+1,m,l) = Vs(2,+1,m,l) + V_lm(lp,mp,-1)*yTy(-1,+1) ! intfeld(typ,lp,mp,-1,l,m,+1)*yTy(-1,+1)
          endif

          lp=l-1
          if (lp.ge.0 .and. lp.le.l_max .and. mp.le.lp) then
             yTy(+1,+1) =  ds*vfa(lp,-mp)*dlts(0)
             yTy(-1,-1) =  ds*vfa(lp,-mp)*dlts(1)
             yTy(+1,-1) = -ds*vfa(lp,-mp)*dlts(1)
             yTy(-1,+1) =  ds*vfa(lp,-mp)*dlts(0)
             Vs(1,+1,m,l) = Vs(1,+1,m,l) + V_lm(lp,mp,+1)*yTy(+1,+1) ! intfeld(typ,lp,mp,+1,l,m,+1)*yTy(+1,+1)
             Vs(1,-1,m,l) = Vs(1,-1,m,l) + V_lm(lp,mp,-1)*yTy(-1,-1) ! intfeld(typ,lp,mp,-1,l,m,-1)*yTy(-1,-1)
             Vs(2,-1,m,l) = Vs(2,-1,m,l) + V_lm(lp,mp,+1)*yTy(+1,-1) ! intfeld(typ,lp,mp,+1,l,m,-1)*yTy(+1,-1)
             Vs(2,+1,m,l) = Vs(2,+1,m,l) + V_lm(lp,mp,-1)*yTy(-1,+1) ! intfeld(typ,lp,mp,-1,l,m,+1)*yTy(-1,+1)
          endif

          mp = m-1
          dlts(:)=1.d0
          if (mp.eq.0) then
             dlts(0) = sqrt(2.d0)
             dlts(1) = 0.d0
          endif
          dlts(:) = dlts(:)*overal

          lp = l+1
          if (mp.ge.0 .and. lp.le.l_max .and. mp.le.lp) then
             yTy(+1,+1) = cs*vfa(l,-m)*dlts(0)
             yTy(-1,-1) = cs*vfa(l,-m)*dlts(1)
             yTy(+1,-1) = cs*vfa(l,-m)*dlts(0)
             yTy(-1,+1) =-cs*vfa(l,-m)*dlts(1)
             Vs(1,+1,m,l) = Vs(1,+1,m,l) + V_lm(lp,mp,+1)*yTy(+1,+1) ! intfeld(typ,lp,mp,+1,l,m,+1)*yTy(+1,+1)
             Vs(1,-1,m,l) = Vs(1,-1,m,l) + V_lm(lp,mp,-1)*yTy(-1,-1) ! intfeld(typ,lp,mp,-1,l,m,-1)*yTy(-1,-1)
             Vs(2,-1,m,l) = Vs(2,-1,m,l) + V_lm(lp,mp,+1)*yTy(+1,-1) ! intfeld(typ,lp,mp,+1,l,m,-1)*yTy(+1,-1)
             Vs(2,+1,m,l) = Vs(2,+1,m,l) + V_lm(lp,mp,-1)*yTy(-1,+1) ! intfeld(typ,lp,mp,-1,l,m,+1)*yTy(-1,+1)
          endif
          
          lp=l-1
          if (mp.ge.0 .and. lp.ge.0 .and. lp.le.l_max .and. mp.le.lp) then
             yTy(+1,+1) = -ds*vfa(lp,mp)*dlts(0)
             yTy(-1,-1) = -ds*vfa(lp,mp)*dlts(1)
             yTy(+1,-1) = -ds*vfa(lp,mp)*dlts(0)
             yTy(-1,+1) =  ds*vfa(lp,mp)*dlts(1)
             Vs(1,+1,m,l) = Vs(1,+1,m,l) + V_lm(lp,mp,+1)*yTy(+1,+1) ! intfeld(typ,lp,mp,+1,l,m,+1)*yTy(+1,+1)
             Vs(1,-1,m,l) = Vs(1,-1,m,l) + V_lm(lp,mp,-1)*yTy(-1,-1) ! intfeld(typ,lp,mp,-1,l,m,-1)*yTy(-1,-1)
             Vs(2,-1,m,l) = Vs(2,-1,m,l) + V_lm(lp,mp,+1)*yTy(+1,-1) ! intfeld(typ,lp,mp,+1,l,m,-1)*yTy(+1,-1)
             Vs(2,+1,m,l) = Vs(2,+1,m,l) + V_lm(lp,mp,-1)*yTy(-1,+1) ! intfeld(typ,lp,mp,-1,l,m,+1)*yTy(-1,+1)
          endif
          ! Z - components
          mp = m
          dlts(:)=1.d0
          if (m.eq.0) then
             dlts(0)=1.d0
             dlts(1)=0.d0
          endif
          lp=l+1
          if (lp.le.l_max .and. mp.le.lp) then
             yTy(+1,+1) = 2*cs*vff(l,m)*dlts(0)
             yTy(-1,-1) = 2*cs*vff(l,m)*dlts(1)
             Vs(3,+1,m,l) = Vs(3,+1,m,l) + V_lm(lp,mp,+1)*yTy(+1,+1) ! intfeld(typ,lp,mp,+1,l,m,+1)*yTy(+1,+1)
             Vs(3,-1,m,l) = Vs(3,-1,m,l) + V_lm(lp,mp,-1)*yTy(-1,-1) ! intfeld(typ,lp,mp,-1,l,m,-1)*yTy(-1,-1)
          endif
          lp=l-1
          if (lp.ge.0 .and. mp.le.lp .and. lp.le.l_max) then
             yTy(+1,+1) = 2*ds*vff(lp,mp)*dlts(0)
             yTy(-1,-1) = 2*ds*vff(lp,mp)*dlts(1)
             Vs(3,+1,m,l) = Vs(3,+1,m,l) + V_lm(lp,mp,+1)*yTy(+1,+1) ! intfeld(typ,lp,mp,+1,l,m,+1)*yTy(+1,+1)
             Vs(3,-1,m,l) = Vs(3,-1,m,l) + V_lm(lp,mp,-1)*yTy(-1,-1) ! intfeld(typ,lp,mp,-1,l,m,-1)*yTy(-1,-1)
          endif
       ENDDO
    ENDDO
    deallocate( V_lm )

    !WRITE(6,'(A,I3,A)') 'Effective Potential at the MT boundary for jatom=', jatom, '='
    !DO lm1=0,lmmax_ext-1
    !   l = abs(lm_ext(1,lm1))
    !   m = lm_ext(2,lm1)
    !   s = sign(1, lm_ext(1,lm1) )
    !   WRITE(6,'(I3,1x,I3,1x,3F20.10,1x)') l*s, m, Vs(1:3,s,m,l)
    !ENDDO
    !WRITE(6,*)

    fri(:)=0.d0
    DO lm1=0,lmmax_ext-1
       l = abs(lm_ext(1,lm1)) 
       m = lm_ext(2,lm1)
       s = sign(1, lm_ext(1,lm1) )
       fri(:) = fri(:) + dble(Vs(:,s,m,l)*rholm(lm1))
    enddo
    fri(:) = fri(:)*Rmt**2
    
    fra(:)=0.d0
    DO lm1=0,lmmax-1
       l = abs(lm(1,lm1)) 
       m = lm(2,lm1)
       s = sign(1, lm(1,lm1) )
       fra(:) = fra(:) + dble(Vs(:,s,m,l)*rhoalm(lm1))
    enddo
    fra(:) = fra(:)*Rmt**2
    deallocate( rholm )
    deallocate( Vs )

    fsur2(1:3) = fsur(1:3)
    fsur2(4:6) = fri(:)-fra(:)
    fsur2(7:9) = fkdc(:)        ! Kinetic energy discontinuity part two computed in the interstitials
    do i=1,9
       if (abs(fsur2(i)).LT.1e-10) fsur2(i)=0.d0
    enddo
    write(6,'(A,I3,A,3F15.7)') 'Fext_intersti[', jatom,']=', fri(:)*1000.d0
    write(6,'(A,I3,A,3F15.7)') 'Fext_mtsphere[', jatom,']=', fra(:)*1000.d0
    !WRITE(6, '(2i3,a7,4e15.7)') jatom,1,'+ SUR', fsur_norm, (fsur2(ik),ik=1,3)
    !WRITE(21,79) jatom,jatom, fsur_norm*1000, (fsur2(ik)*1000,ik=1,3)
    !79  FORMAT(':FSX',i3.3,':',1x,i3,'.ATOM',4F17.9) 
  END Function Force_surface_extra
  
End module Forces


