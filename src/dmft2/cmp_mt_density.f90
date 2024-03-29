! @Copyright 2007 Kristjan Haule
! 

SUBROUTINE cmp_MT_density(w_RHOLM, Aweight, wEpsw, weight, fsph, fnsp, fsph2, nemin, nemax, DM_nemin, DM_nemax, DM_nemaxx, nbands, nbandsx, nbands_dft, lm_max, n0, nnlo, coord, ikp, DM_EF, time_bl,time_bl_w, time_reduc, time_reduc_w, time_radprod, time_radprod_w,time_m, time_m_w, time_rad, time_rad_w, time_ilm, time_ilm_w, time_lo, time_lo_w, time_force, time_force_w)
  USE defs,  ONLY: pi, imag
  USE param, ONLY: lmax2, nume, nloat, iblock, nmat, nrad, CHUNK, lomax, ngau
  USE atspdt,ONLY: p, dp, pe, dpe, el
  USE com,   ONLY: nat, rel
  USE lo,    ONLY: nlo, nlov, nlon, lapw, ilo, loor, a1lo, b1lo, elo
  USE lohelp,ONLY: u21, ue21, u12, ue12, u22, sum12, sum21, sum22, sume21, sume12, zerotc_lohelp
  USE structure, ONLY: Rmt, vol, mult, rotij, br1, pos, tauij, jri, rotloc, natm
  USE xa,    ONLY: TC100, TCA100, TCB100, LM, E, SUMA, SUMB, SUMAB, SUMBA
  USE w_atpar, ONLY: w_xwt1, w_xwt1l, w_xwt1h, w_xwteh, w_xwtel, retrieve_w_atpar, Vts, lpv, lv, mpv, mv, Nhmx
  USE dmfts, ONLY: iso, shft
  USE charp, ONLY: zero_charp
  USE chard, ONLY: zero_chard
  USE charf, ONLY: zero_charf
  USE w_atpar,ONLY: ri_mat, ul_Rmt, dul_Rmt
  USE alm_blm_clm, ONLY: cmp_alm_blm_clm
  USE Forces, ONLY: Force1_corrected, Force1_DMFT, Force2, Force3, forcea, Qforce !, Force5_mine
  USE com_mpi, ONLY : myrank, master
  IMPLICIT NONE
  REAL*8, intent(out)    :: w_RHOLM(NRAD,LM_MAX,nat)!, w_vRHOLM(NRAD,LM_MAX,nat)
  INTEGER, intent(in)    :: nemin, nemax, DM_nemin, DM_nemax, DM_nemaxx, nbands, nbandsx, nbands_dft, lm_max, n0, nnlo, ikp
  COMPLEX*16, intent(in) :: Aweight(nbands,nbands), wEpsw(nbands_dft,nbands_dft)
  REAL*8, intent(in)     :: weight(nume), DM_EF
  REAL*8, intent(inout)  :: time_bl, time_bl_w, time_reduc, time_reduc_w, time_radprod, time_radprod_w, time_m, time_m_w, time_rad, time_rad_w, time_ilm, time_ilm_w, time_lo, time_lo_w, time_force, time_force_w
  CHARACTER*5, intent(in):: coord
  REAL*8,PARAMETER       :: Ry2eV= 13.60569253d0
  !
  REAL*8, intent(inout) :: fsph(3,natm)
  REAL*8, intent(inout) :: fnsp(3,natm)
  REAL*8, intent(inout) :: fsph2(3,natm)
  !
  Interface
     integer FUNCTION NOTRI(k,l,m)
       integer, intent(in) :: k,l,m
     end FUNCTION NOTRI
     real*8 FUNCTION GAUNT(l1,l2,l3,mm1,mm2,mm3)
       integer, intent(in) :: l1,l2,l3,mm1,mm2,mm3
     end FUNCTION GAUNT
     SUBROUTINE YLM(V,LMAX,Y)
       COMPLEX*16, intent(out)      :: Y(*)
       INTEGER, intent(in)          :: LMAX
       DOUBLE PRECISION, intent(in) :: V(3)
     end SUBROUTINE YLM
  end interface
  ! commons
  COMMON /RADFU/  RRAD1(NRAD,0:LMAX2),RADE1(NRAD,0:LMAX2),RRAD2(NRAD,0:LMAX2),RADE2(NRAD,0:LMAX2)
  REAL*8     :: RRAD1, RADE1, RRAD2, RADE2
  COMMON /UHELP/  UA(NRAD),UB(NRAD),UBA(NRAD),UAB(NRAD)
  REAL*8     :: UA, UB, UBA, UAB
  ! locals
  COMPLEX*16, ALLOCATABLE :: alm(:,:),blm(:,:),clm(:,:,:), dh_alm(:,:,:), dh_blm(:,:,:), dh_clm(:,:,:,:)
  COMPLEX*16 :: imag1, gint, wkp
  INTEGER :: jatom, lfirst, lmmax, mmax, latom, mu, is, i, l, m, lmx, l1, ilm, li, mi, lp, lp1, mpmax, imax
  INTEGER :: jlo, jlop, num, mp, ms, mps, mtest, ly, lpy, jatombad, lbad, lmax2lmax2, lomaxlomax
  REAL*8  :: CIN, test1, eqbad!, facv
  REAL*8  :: time1, time1_w, time2, time2_w, time3, time3_w, time4, time4_w
  COMPLEX*16 :: sa, sb, sab, sba, tsuma, tsumb, tsumab, tsumba!, vsa, vsb, vsab, vsba
  COMPLEX*16 :: s12(nloat), se12(nloat), s21(nloat), se21(nloat), s22(nloat,nloat)
  !COMPLEX*16 :: vs12(nloat), vse12(nloat), vs21(nloat), vse21(nloat), vs22(nloat,nloat)
  COMPLEX*16 :: tsum21(nloat), tsume21(nloat), tsum12(nloat), tsume12(nloat), tsum22(nloat,nloat)
  COMPLEX*16 :: tsa, tsb, tsab, tsba, ts12(nloat), tse12(nloat), ts21(nloat), tse21(nloat), ts22(nloat,nloat)
  REAL*8     :: tc_buf(nrad)!,vtc_buf(nrad)
  LOGICAL    :: Qforce_j
  COMPLEX*16, ALLOCATABLE :: aalm(:,:,:), bblm(:,:,:), cclm(:,:,:,:)
  REAL*8, allocatable :: GNT(:,:)
  CHARACTER*20 :: cstr
  INTEGER    :: j, ik
  INTEGER    :: lmax2lmax2_,nbands_dft_,dim_,lomaxlomax_,nloat_
  
  lmax2lmax2 = (LMAX2+1)*(LMAX2+1)
  lomaxlomax = (lomax+1)*(lomax+1)
  LMX=LMAX2                                                    
  CIN=1.d0/137.0359895d0**2
  IF (.NOT.REL) CIN=4.0*1.0D-22                                     

  !nbands_dft = nemax-nemin+1

  test1=0.0d0
  ALLOCATE(alm(lmax2lmax2,nume),blm(lmax2lmax2,nume),clm(lomaxlomax,nume,nloat))
  ALLOCATE( dh_alm(iso,lmax2lmax2,nume), dh_blm(iso,lmax2lmax2,nume), dh_clm(nloat,iso,lomaxlomax,nume) )

  if (Qforce) then
     lmax2lmax2_ = lmax2lmax2
     nbands_dft_ = nbands_dft
     dim_ = 3
     lomaxlomax_ = lomaxlomax
     nloat_ = nloat
  else
     ! the arrays need to be always allocated     
     lmax2lmax2_ = 1
     nbands_dft_ = 1
     dim_ = 1
     lomaxlomax_ = 1
     nloat_ = 1
  endif
  ALLOCATE( aalm(lmax2lmax2_,nbands_dft_,dim_), bblm(lmax2lmax2_,nbands_dft_,dim_), cclm(lomaxlomax_,nbands_dft_,nloat_,dim_) )

  do jatom=1,nat
     ! Quantities computed by "atpar" which do not depend on k-point
     CALL retrieve_w_atpar(jatom,lfirst,lmmax)
     ! zero koeff. for charge analysis
     CALL zero_charp(nemin,nemax)
     CALL zero_chard(nemin,nemax)
     CALL zero_charf(nemin,nemax)
     CALL zerotc_lohelp(nemin,nemax)
     TC100(0:lmax2,nemin:nemax)=0.0d0
     TCA100(0:lmax2,nemin:nemax)=0.0d0
     TCB100(0:lmax2,nemin:nemax)=0.0d0

     ! CALCULATE ALM, BLM                                                
     latom=lfirst-1                                                    

     Qforce_j = Qforce .AND. forcea(0,jatom)

     DO mu=1,mult(jatom)
        latom=latom+1
        DO is=1,iso !!! over spin: Mar 25, 2010

           CALL cputim(time1)
           CALL walltim(time1_w)
           !CALL cmp_alm_blm_clm(alm,blm,clm,aalm,bblm,cclm,Aweight,latom,lfirst,jatom,nemin,nemax,n0,nnlo,is,DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,Qforce_j,lmax2lmax2_,nbands_dft_,dim_,lomaxlomax_,nloat_,time_bl,time_bl_w,time_lo,time_lo_w)
           CALL cmp_alm_blm_clm(alm,blm,clm,aalm,bblm,cclm,Aweight,latom,lfirst,jatom,nemin,nemax,n0,nnlo,is,DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,Qforce_j,time_bl,time_bl_w,time_lo,time_lo_w)
           
           dh_alm(is,:,:) = alm(:,:)
           dh_blm(is,:,:) = blm(:,:)
           DO jlo=1,nloat
              dh_clm(jlo,is,:,:) = clm(:,:,jlo)
           ENDDO

           CALL cputim(time2)
           CALL walltim(time2_w)

           ! MAIN FORCE CALCULATIONS
           if (Qforce_j) then
              !fsph(:,latom)  = fsph(:,latom)  + Force1(nemin,nemax,DM_nemaxx,lmx,alm,blm,clm,aalm,bblm,cclm,weight,el,elo,E,ri_mat(:,:,:,jatom),ilo,lmax2,nume,nloat,lomax,nbands_dft)
              !fsph(:,latom)  = fsph(:,latom)  + Force1_corrected(nemin,nemax,DM_nemaxx,lmx,alm,blm,clm,aalm,bblm,cclm,weight,el,elo,E,ri_mat(:,:,:,jatom),ilo,lmax2,nume,nloat,lomax,nbands_dft)
              fsph(:,latom)  = fsph(:,latom)  + Force1_DMFT(nemin,nemax,DM_nemaxx,lmx,alm,blm,clm,aalm,bblm,cclm,weight,el,elo,wEpsw,ri_mat(:,:,:,jatom),ilo,lmax2,nume,nloat,lomax,nbands_dft,nbands)
              fnsp(:,latom)  = fnsp(:,latom)  + Force2(nemin,nemax,DM_nemaxx,Nhmx(jatom),alm,blm,clm,aalm,bblm,cclm,weight,Vts(:,:,:,jatom),lv(:,jatom),lpv(:,jatom),mv(:,jatom),mpv(:,jatom),ilo,lmax2,lomax,nume,nloat,ngau,nbands_dft)
              fsph2(:,latom) = fsph2(:,latom) + Force3(nemin,nemax,DM_nemaxx,lmx,alm,blm,clm,aalm,bblm,cclm,weight,ul_Rmt(:,:,jatom),dul_Rmt(:,:,jatom),Rmt(jatom),ilo,lmax2,lomax,nume,nloat,nbands_dft)
           endif
           
           CALL cputim(time3)
           CALL walltim(time3_w)
           time_reduc=time_reduc+time2-time1
           time_reduc_w=time_reduc_w+time2_w-time1_w
           time_force=time_force+time3-time2
           time_force_w=time_force_w+time3_w-time2_w
        enddo  !!! Mar 27, 2010: loop over spin
        
        !! Mar 27, 2010: we need alm for csplit below. Take the average between up and down
        if (iso.eq.2) then
           alm(:,:) = (dh_alm(1,:,:)+dh_alm(2,:,:))/2.
           blm(:,:) = (dh_blm(1,:,:)+dh_blm(2,:,:))/2.
           DO jlo=1,nloat
              clm(:,:,jlo) = (dh_clm(jlo,1,:,:)+dh_clm(jlo,2,:,:))/2.
           ENDDO
        endif

        !.....C(L,M) TO BE CALCULATED                                           
        CALL cputim(time1)
        CALL walltim(time1_w)

        DO ILM=1,LMMAX
           Li=IABS(LM(1,ILM))
           Mi=LM(2,ILM)
           !     SEE KURKI-SUONI  FACTOR FOR  LM+ (1,0), LM-(0,-1), M.NE.2N *(-1)
           IMAG1=(1.0d0,0.0d0)                                                   
           IF(LM(1,ILM).LT.0) IMAG1= IMAG   ! (1,i)*(-1)**m for s=(1,-1)
           IF(MOD(Mi,2).EQ.1) IMAG1=-IMAG1  ! (-1)^m*i

           l_sum: DO L=0,lmx
              L1=L+1
              MMAX=2*L+1
              lp_sum: DO LP=0,lmx
                 IF(NOTRI(Li,L,LP).LT.0) CYCLE

                 CALL cputim(time3)
                 CALL walltim(time3_w)
                 ! It is a good idea to recompute gaunt coefficients for this case as it takes some time to compute for each band.
                 allocate( GNT(-L:L, -LP:LP) )
                 GNT(:,:)=0.d0
                 do M=-L,L
                    do MP=-LP,LP
                       if (MP .EQ. M-Mi) GNT(M,MP) = GAUNT(L,Li,LP,M,Mi,MP)
                    enddo
                 enddo

                 LP1=LP+1
                 MPMAX=2*LP+1
                 IMAX=JRI(JATOM)
                 ! Radial functions are constructed
                 DO I=1,IMAX
                    UA(I) =RRAD1(I,L)*RRAD1(I,LP)+CIN*RRAD2(I,L)*RRAD2(I,LP)    
                    UB(I) =RADE1(I,L)*RADE1(I,LP)+CIN*RADE2(I,L)*RADE2(I,LP)    
                    UAB(I)=RRAD1(I,L)*RADE1(I,LP)+CIN*RRAD2(I,L)*RADE2(I,LP)   
                    UBA(I)=RRAD1(I,LP)*RADE1(I,L)+CIN*RRAD2(I,LP)*RADE2(I,L)   
                 ENDDO
                 ! Radial functions as product of local orbitals and u,dot{u}
                 DO jlo=1,ilo(l)
                    IF(loor(jlo,l)) THEN
                       DO I=1,IMAX
                          U21(I,jlo) =a1lo(I,jlo,L)*RRAD1(I,LP)+CIN*b1lo(I,jlo,L)*RRAD2(I,LP)    
                          Ue21(I,jlo)=a1lo(I,jlo,L)*RADE1(I,LP)+CIN*b1lo(I,jlo,L)*RADE2(I,LP)    
                       ENDDO
                    ENDIF
                 ENDDO
                 DO jlop=1,ilo(lp)
                    IF(loor(jlop,lp)) THEN
                       DO i=1,imax
                          u12(i,jlop) =RRAD1(i,l)*a1lo(i,jlop,lp)+CIN*RRAD2(i,l)*b1lo(i,jlop,lp)    
                          ue12(i,jlop)=RADE1(i,l)*a1lo(i,jlop,lp)+CIN*RADE2(i,l)*b1lo(i,jlop,lp) 
                       ENDDO
                    ENDIF
                 ENDDO
                 ! Radial functions for local orbitals
                 DO jlo=1,ilo(l)
                    DO jlop=1,ilo(lp) 
                       IF(loor(jlo,l).AND.loor(jlop,lp)) THEN
                          DO i=1,imax
                             u22(i,jlop,jlo)=a1lo(i,jlo,l)*a1lo(i,jlop,lp)+cin*b1lo(i,jlo,l)*b1lo(i,jlop,lp) 
                          ENDDO
                       ENDIF
                    ENDDO
                 ENDDO
                 CALL cputim(time4)
                 CALL walltim(time4_w)
                 time_radprod=time_radprod+time4-time3
                 time_radprod_w=time_radprod_w+time4_w-time3_w
                 !.....M SUM                

                 sa=0.0d0; sb=0.0d0; sab=0.0d0; sba=0.0d0; s12=0.0d0; 
                 se12=0.0d0; s21=0.0d0; se21=0.0d0; s22=0.0d0

                 CALL cputim(time3)
                 CALL walltim(time3_w)

!!! This loop takes really a lot of time, but the multithreading is dissabled, because the code can corredump.
!!! You should debug this very carefully.
                 !$OMP PARALLEL DO SHARED(suma,sumb,sumab,sumba,sum21,sume21,sum12,sume12,sum22,GNT)&
                 !$OMP& PRIVATE(num,M,MP,MS,MPS,MTEST,LY,LPY,is,gint,tsuma,tsumb,tsumab,tsumba,tsum21,tsume21,tsum12,tsume12,tsum22,tsa,tsb,tsab,tsba,ts21,tse21,ts12,tse12,ts22,wkp)&
                 !$OMP& SCHEDULE(STATIC,CHUNK) &
                 !$OMP& REDUCTION(+:sa,sb,sab,sba,s21,se21,s12,se12,s22)
                 DO num=nemin,DM_nemaxx
                    tsuma = 0.0d0; tsumb = 0.0d0; tsumab= 0.0d0; tsumba= 0.0d0
                    tsum21 (:ilo(l))  = 0.0d0;  tsume21(:ilo(l))  = 0.0d0
                    tsum12 (:ilo(lp)) = 0.0d0;  tsume12(:ilo(lp)) = 0.0d0
                    tsum22 (:ilo(lp),:ilo(l))=0.0d0
                    M=-(L+1) 
                    MP=-(LP+1)
                    DO MS=1,2*L+1
                       M=M+1
                       MP=-LP1
                       DO MPS=1,2*LP+1
                          MP=MP+1
                          MTEST=-M+Mi+MP
                          IF(MTEST.NE.0) CYCLE
                          LY  = L*(L+1)+M+1
                          LPY = LP*(LP+1)+MP+1
                          ! GNT=GAUNT(L,Li,LP,M,Mi,MP)
                          gint = IMAG1*GNT(M,MP)
                          do is=1,iso
                             tSA = dh_alm(is,LY,num)*dconjg(dh_alm(is,LPY,num))*gint
                             tSB = dh_blm(is,LY,num)*dconjg(dh_blm(is,LPY,num))*gint
                             tSAB= dh_alm(is,LY,num)*dconjg(dh_blm(is,LPY,num))*gint
                             tSBA= dh_blm(is,LY,num)*dconjg(dh_alm(is,LPY,num))*gint
                             DO jlo=1,ilo(l)
                                tS21(jlo)  = dh_clm(jlo,is,LY,num)*dconjg(dh_alm(is,LPY,num))*gint
                                tSe21(jlo) = dh_clm(jlo,is,LY,num)*dconjg(dh_blm(is,LPY,num))*gint
                             ENDDO
                             DO jlop=1,ilo(lp)
                                tS12(jlop)  = dh_alm(is,LY,num)*dconjg(dh_clm(jlop,is,LPY,num))*gint
                                tSe12(jlop) = dh_blm(is,LY,num)*dconjg(dh_clm(jlop,is,LPY,num))*gint
                             ENDDO
                             DO jlo=1,ilo(l)
                                DO jlop=1,ilo(lp)
                                   tS22(jlop,jlo) = dh_clm(jlo,is,LY,num)*dconjg(dh_clm(jlop,is,LPY,num))*gint
                                ENDDO
                             ENDDO

                             tSUMA  = tSUMA  + tSA/iso
                             tSUMB  = tSUMB  + tSB/iso
                             tSUMAB = tSUMAB + tSAB/iso
                             tSUMBA = tSUMBA + tSBA/iso
                             DO jlo=1,ilo(l)
                                tSUM21 (jlo) = tSUM21 (jlo) + tS21 (jlo)/iso
                                tSUMe21(jlo) = tSUMe21(jlo) + tSe21(jlo)/iso
                             ENDDO
                             DO jlop=1,ilo(lp)
                                tSUM12 (jlop) = tSUM12 (jlop) + tS12(jlop)/iso
                                tSUMe12(jlop) = tSUMe12(jlop) + tSe12(jlop)/iso
                             ENDDO
                             DO jlo=1,ilo(l)
                                DO jlop=1,ilo(lp)
                                   tSUM22(jlop,jlo) = tSUM22(jlop,jlo) + tS22(jlop,jlo)/iso
                                ENDDO
                             ENDDO
                          enddo

                       ENDDO
                    ENDDO

                    wkp = weight(num)
                    sa  = sa  + tsuma  * wkp
                    sb  = sb  + tsumb  * wkp
                    sab = sab + tsumab * wkp
                    sba = sba + tsumba * wkp
                    DO jlo=1,ilo(l)
                       s21(jlo)   = s21(jlo)  + tsum21 (jlo) * wkp
                       se21(jlo)  = se21(jlo) + tsume21(jlo) * wkp
                    ENDDO
                    DO jlop=1,ilo(lp)
                       s12  (jlop) = s12  (jlop) + tsum12 (jlop) * wkp
                       se12 (jlop) = se12 (jlop) + tsume12(jlop) * wkp
                    ENDDO
                    DO jlo=1,ilo(l)
                       DO jlop=1,ilo(lp)
                          s22 (jlop,jlo) = s22 (jlop,jlo) + tsum22(jlop,jlo) * wkp
                       ENDDO
                    ENDDO

                    suma  (num) = tsuma
                    sumb  (num) = tsumb
                    sumab (num) = tsumab
                    sumba (num) = tsumba
                    sum21 (num,:ilo(l))  = tsum21 (:ilo(l))
                    sume21(num,:ilo(l))  = tsume21(:ilo(l))
                    sum12 (num,:ilo(lp)) = tsum12 (:ilo(lp))
                    sume12(num,:ilo(lp)) = tsume12(:ilo(lp))
                    sum22 (num,:ilo(lp),:ilo(l)) = tsum22 (:ilo(lp),:ilo(l))
                 ENDDO
                 !$OMP END PARALLEL DO

                 deallocate( GNT )

                 CALL cputim(time4)
                 CALL walltim(time4_w)
                 time_m=time_m+time4-time3
                 time_m_w=time_m_w+time4_w-time3_w

                 IF(Li.eq.0) then
                    call csplit(nemin,DM_nemaxx,l,jatom,mu,alm,blm,clm,coord) 
                 ENDIF
                 CALL cputim(time3)
                 CALL walltim(time3_w)
                 IMAX=JRI(JATOM)

                 DO I=1,IMAX
                    TC_buf(i)=SA*UA(I)+SB*UB(I)+SAB*UAB(I)+SBA*UBA(I)
                 ENDDO
                 DO jlo=1,ilo(l)
                    DO i=1,imax
                       tc_buf(i)=tc_buf(i)+S21(jlo)*U21(I,jlo)+Se21(jlo)*Ue21(I,jlo)
                    ENDDO
                 ENDDO
                 DO jlop=1,ilo(lp)
                    DO i=1,imax
                       tc_buf(i)=tc_buf(i)+S12(jlop)*U12(I,jlop)+Se12(jlop)*Ue12(I,jlop)
                    ENDDO
                 ENDDO
                 DO jlo=1,ilo(l)
                    DO jlop=1,ilo(lp)
                       DO i=1,imax
                          tc_buf(i)=tc_buf(i)+s22(jlop,jlo)*u22(i,jlop,jlo)
                       ENDDO
                    ENDDO
                 ENDDO

                 DO i=1,imax
                    w_RHOLM(I,ILM,jatom)=w_RHOLM(I,ILM,jatom)+TC_buf(i)/MULT(JATOM)
                 ENDDO
                 CALL cputim(time4)
                 CALL walltim(time4_w)
                 time_rad=time_rad+time4-time3
                 time_rad_w=time_rad_w+time4_w-time3_w
              ENDDO lp_sum
           ENDDO l_sum
        ENDDO

        CALL cputim(time2)
        CALL walltim(time2_w)
        time_ilm=time_ilm+time2-time1
        time_ilm_w=time_ilm_w+time2_w-time1_w
     ENDDO
     call psplit(w_xwt1(:,jatom),w_xwteh(:,jatom),w_xwtel(:,jatom),w_xwt1h(:,jatom),w_xwt1l(:,jatom),jatom,nemin,DM_nemaxx,test1,EQBAD,jatombad,lbad)
  ENDDO


  IF(test1.GT.15.d0) THEN
     IF(myrank.eq.master) WRITE(6,1001)  test1,(EQBAD-DM_EF)*Ry2eV, jatombad,lbad,test1/100.
     IF(myrank.eq.master) WRITE(21,1001) test1,(EQBAD-DM_EF)*Ry2eV, jatombad,lbad,test1/100.
  else IF(test1.GT.2.d0) THEN
     IF(myrank.eq.master) WRITE(6,1002)  test1,(EQBAD-DM_EF)*Ry2eV, jatombad,lbad,test1/100.
     IF(myrank.eq.master) WRITE(21,1002) test1,(EQBAD-DM_EF)*Ry2eV, jatombad,lbad,test1/100.
  ENDIF



  DEALLOCATE( aalm, bblm, cclm )
  DEALLOCATE(alm,blm,clm)
  DEALLOCATE(dh_alm,dh_blm,dh_clm)


1001 FORMAT(//,'   QTL-B VALUE .EQ. ',F10.5,' in Band of energy ',F9.4,'eV from EF ATOM=',i5,2x,'L=',i3,/,'    Check for ghostbands or EIGENVALUES BELOW XX messages',/, &
              '    The projection of charge to \dot{u}(r) is',F10.5,', which is large. Adjust your Energy-parameters for this ATOM and L (or use -in1new switch), check RMTs  !!!',//)
1002 FORMAT(//,'   QTL-B VALUE .EQ. ',F10.5,' in Band of energy ',F9.4,'eV from EF ATOM=',i5,2x,'L=',i3,/, &
               '   The projection of charge to \dot{u}(r) is',F10.5,', which is large. Most likely no ghostbands, but adjust Energy-parameters for this ATOM and L',//)

END SUBROUTINE cmp_MT_density

