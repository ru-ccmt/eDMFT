SUBROUTINE  HNS(NV,NAT,NNSL,POS,Vol,Rmt,MULT,LQIND,LQNS,GFAC,ROTLOC,ROTIJ,VNS1,VNS2,VNS3)
  use matrices, only : HS, HSROWS, KZZ, Kn
  use albl, only : AL, BL, init_albl, end_albl
  use loabc, only : ALO, BLO, CLO, nloat_new
  use lolog, only : nlo, loor, ilo, nlov
  use loint, only : VNS1LO, VNS2LO, VNS3LO
  use structure, only : ndf
  use lapw_timer, only : READ_CPU_TIME, READ_WALL_TIME, time_pwpw, time_pwlo,time_lolo,time_hns30,  START_TIMER, STOP_TIMER, time_hns
  use param
  use mpi, only : Qprint
  IMPLICIT NONE
  !        Arguments
  INTEGER, intent(in)    :: NV, NAT, NNSL, mult(nat), LQIND(NAT), LQNS(6,NGAU,NAT)
  REAL*8,  intent(in)    :: POS(3,NDF), Vol, Rmt(NAT), ROTLOC(3,3,NAT), ROTIJ(3,3,NDF)
  COMPLEX*16, intent(in) :: GFAC(NGAU,NAT)
  REAL*8,  intent(in)    :: VNS1(0:NSLMAX-1,LMMX,0:NSLMAX-1,NAT), VNS2(0:NSLMAX-1,LMMX,0:NSLMAX-1,NAT), VNS3(0:NSLMAX-1,LMMX,0:NSLMAX-1,NAT)
  ! Locals
  INTEGER    :: jlo, jlo0, jlop, iml, i, i1, i2, j, jatom, mu, LINDEX, L0M0, LPMP, L0, LP, M0, MP, LMX, L0IDX, LMXIDX, LPIDX, M0IDX, MPIDX, latom, lfirst
  REAL*8     :: ARGX, RHOATM, PI, ROTV2(3), rotloc_x_BR1_x_rotij(3,3)
  COMPLEX*16 :: TMPC, xlmc, i_to_l, TUU, TDU, TUD, TDD, GNTFLD(NGAU)
  INTEGER    :: LMXPTR(NGAU+2), LPPTR(NGAU+2), L0PTR(NGAU+2), M0PTR(NGAU+2), MPPTR(NGAU+2), L0FLD(NGAU), LMXFLD(NGAU), LPFLD(NGAU), M0FLD(NGAU), MPFLD(NGAU)
  LOGICAL    :: loor_l(0:NNSL-1)
  !
  COMPLEX*16, allocatable :: TUU12(:), TUU21(:), TUU22(:,:), TUD21(:), TDU12(:), PHS(:), YL(:,:)
  COMPLEX*16, allocatable :: ALM(:,:), BLM(:,:), CLM(:,:)
  REAL*8,     allocatable :: A1R(:,:), A1I(:,:), A2R(:,:), A2I(:,:), B1R(:,:), B1I(:,:), B2R(:,:), B2I(:,:), C1R(:,:,:), C1I(:,:,:), C2R(:,:,:), C2I(:,:,:)
  INTEGER,    allocatable :: i2jlo(:), iind(:,:)
#ifndef _REAL  
  COMPLEX*16, ALLOCATABLE :: A1(:,:), A2(:,:)
#endif  
  !        Local Parameters
  COMPLEX*16 :: CIMAG, CONE
  REAL*8     :: ONE
  INTEGER    :: NSL2
  PARAMETER  (CIMAG  = (0.0D0,1.0D0), CONE=(1.0D0,0.0D0), ONE=1.0D0 )
  PARAMETER  (NSL2=NSLMAX*NSLMAX)
  !
  DOUBLE PRECISION    tot_time_alm, tot_time_vector, tot_time_vector2, tot_time_v1xv2,tot_cputime_copy2abclm,tot_cputime_PDSYR2M,tot_walltime_copy2abclm,tot_walltime_PDSYR2M
  !        External Subroutines
  EXTERNAL           LMSORT, YLM
  !
  CALL START_TIMER(time_hns)
  nloat=nloat_new
  PI = 4.d0*ATAN(1.d0)
  
  allocate( A1R(HSROWS,NSL2), A1I(HSROWS, NSL2), A2R(HSROWS, NSL2), A2I(HSROWS, NSL2) )
  allocate( B1R(HSROWS,NSL2), B1I(HSROWS, NSL2), B2R(HSROWS, NSL2), B2I(HSROWS, NSL2) )
  allocate( C1R(HSROWS,NSL2,nloat), C1I(HSROWS,NSL2,nloat), C2R(HSROWS,NSL2,nloat), C2I(HSROWS, NSL2,nloat) )
  allocate( ALM(HSROWS,0:NSL2-1), BLM(HSROWS, 0:NSL2-1), CLM(HSROWS, 0:NSL2-1) )
  !
  allocate( YL(0:NSL2-1,HSROWS),PHS(HSROWS) )
  allocate(i2jlo(HSROWS))
  !
  C1R(:,:,:)=0.0d0
  C1I(:,:,:)=0.0d0
  
#ifndef _REAL  
  allocate( A1(HSROWS, NSL2) )
  allocate( A2(HSROWS, NSL2) )
#endif
  
  call init_albl(NNSL-1,'HNS')
  
  allocate(TUU21(0:nloat),TUD21(0:nloat),TUU12(0:nloat),TDU12(0:nloat), TUU22(0:nloat,0:nloat))
  TUU21(0:nloat)=0.0d0
  TUD21(0:nloat)=0.0d0
  TUU12(0:nloat)=0.0d0
  TDU12(0:nloat)=0.0d0
  TUU22(0:nloat,0:nloat)=0.0d0

  tot_time_alm=0.0d0
  tot_time_vector=0.0d0
  tot_time_vector2=0.0d0
  tot_time_v1xv2=0.0d0
  tot_cputime_copy2abclm=0.0d0
  tot_walltime_copy2abclm=0.0d0
  tot_cputime_PDSYR2M=0.0d0
  tot_walltime_PDSYR2M=0.0d0
  
  i2jlo(:)=0
  
  allocate( iind(0:lomax,nloat) )
  
  jatom_loop: DO jatom=1,nat
     ! precompute index for local orbitals iind(l,jlo)
     iind(:,:)=0
     i=nlov(jatom) ! local orbitals before this atom
     DO l0=0,lomax
        do jlo=1,ilo(l0,jatom)
           iind(l0,jlo) = i+1    ! first local orbital for (jatom,l,jlo), i.e., and m=-l and mu=1
           i = i + mult(jatom)*(2*l0+1)
        end do
     END DO
     
     ! precompute al and bl
     lfirst = 1+sum(mult(1:jatom-1))
     call make_albl(jatom,Rmt(jatom),NNSL-1,NV,'HNS',1,0)
     CALL LMSORT(NAT,NGAU,JATOM,LQIND,LQNS,GFAC,L0PTR,L0FLD,M0PTR,M0FLD,LPPTR,LPFLD,MPPTR,MPFLD,LMXPTR,LMXFLD,GNTFLD)
     RHOATM = 4.d0*PI*Rmt(JATOM)**2/SQRT(Vol)

     loor_l(0:lomax) = loor(0:lomax,jatom)
     loor_l(lomax+1:)=.false.
     
     mu_loop: DO mu=1,mult(jatom)
        call start_timer(time_pwpw)
        latom = lfirst + mu -1
        rotloc_x_BR1_x_rotij(:,:) = matmul(rotloc(:,:,jatom), rotij(:,:,latom))
        do i=1,NV+NLO
           ARGX = 2.d0*Pi*( POS(1,latom)*KZZ(1,i) + POS(2,latom)*KZZ(2,i) + POS(3,latom)*KZZ(3,i))
           PHS(i) = dcmplx( dcos(ARGX), dsin(ARGX) )
           ROTV2 = matmul(rotloc_x_BR1_x_rotij, Kn(:,i))
           CALL YLM (ROTV2,NNSL-1,YL(0,i))
        enddo
        !  Precompute Alm and Blm and Clm
        ALM(:,:)=0.d0
        BLM(:,:)=0.d0
        CLM(:,:)=0.d0
        DO L0 = 0, NNSL - 1
           i_to_l = cimag**l0
           DO M0=-L0,L0
              iml = l0*(l0+1)+m0
              do i=1,NV
                 xlmc = RHOATM*i_to_l*PHS(i)*dconjg(YL(iml,i))
                 ALM(i,iml) = AL(i,l0,1)*xlmc
                 BLM(i,iml) = BL(i,l0,1)*xlmc
              enddo
              IF (L0.LE.LOMAX) THEN
                 do jlo=1,ilo(l0,jatom)
                    ! At fixed (jatom,l,jlo) we go over all [mu',m']
                    i1 = iind(l0,jlo)
                    i2 = i1 + (2*l0+1)*mult(jatom)-1
                    do i=NV+i1,NV+i2
                       xlmc = RHOATM*i_to_l*PHS(i)*dconjg(YL(iml,i))
                       ALM(i,iml) = ALO(l0,jlo,JATOM)*xlmc
                       BLM(i,iml) = BLO(l0,jlo,JATOM)*xlmc
                       CLM(i,iml) = CLO(l0,jlo,JATOM)*xlmc
                       i2jlo(i)=jlo
                    enddo
                 enddo
              ENDIF
           enddo
        enddo
        call stop_timer(time_pwpw)
        tot_time_alm=tot_time_alm+READ_CPU_TIME(time_pwpw)
        
        A1R(:,:)=0.d0; A1I(:,:)=0.d0; B1R(:,:)=0.d0; B1I(:,:)=0.d0; C1R(:,:,:)=0.d0; C1I(:,:,:)=0.d0
        call start_timer(time_pwlo)
        !  determine nonspherical contribution to Hamilton-matrix elements
        ! ************* initialize Loop-Index
        LINDEX=0
        DO L0IDX = L0PTR(1), L0PTR(2) - 1
           L0 = L0FLD(L0IDX)    
           DO M0IDX = M0PTR(L0IDX), M0PTR(L0IDX+1) - 1
              ! ************ increment Loop-Index 
              LINDEX=LINDEX+1
              M0 = M0FLD(M0IDX)
              L0M0 = L0*(L0+1) + M0
              !*** First set A1,B1,C1 matrices
              do j=1,NV+NLO
                 jlo=i2jlo(j)
                 A1R(j,LINDEX) = dble(ALM(J,L0M0))
                 A1I(j,LINDEX) = dimag(ALM(J,L0M0))
                 B1R(j,LINDEX) = dble(BLM(J,L0M0))
                 B1I(j,LINDEX) = dimag(BLM(J,L0M0))
                 if (jlo.gt.0) then
                    C1R(j,LINDEX,jlo) = dble(CLM(J,L0M0))
                    C1I(j,LINDEX,jlo) = dimag(CLM(J,L0M0))
                 endif
              enddo
              !*** Now working on A2,B2,C2
              A2R(1:NV+NLO,LINDEX) = 0.d0
              A2I(1:NV+NLO,LINDEX) = 0.d0
              B2R(1:NV+NLO,LINDEX) = 0.d0
              B2I(1:NV+NLO,LINDEX) = 0.d0
              C2R(1:NV+NLO,LINDEX,1:nloat) = 0.d0
              C2I(1:NV+NLO,LINDEX,1:nloat) = 0.d0
              
              DO LPIDX = LPPTR(M0IDX), LPPTR(M0IDX+1) - 1
                 LP = LPFLD(LPIDX)
                 DO MPIDX = MPPTR(LPIDX), MPPTR(LPIDX+1) - 1
                    MP = MPFLD(MPIDX)
                    LPMP = LP*(LP + 1) + MP
                    !
                    TUU = 0.d0; TDU = 0.d0; TUD = 0.d0; TDD = 0.d0
                    TUU21 = 0.d0; TUD21 = 0.d0; TUU12 = 0.d0; TDU12 = 0.d0; TUU22 = 0.d0
                    !
                    DO LMXIDX = LMXPTR(MPIDX), LMXPTR(MPIDX+1)-1
                       LMX = LMXFLD(LMXIDX)
                       TUU = TUU + GNTFLD(LMXIDX)*VNS1(LP,LMX,L0,JATOM)
                       TDU = TDU + GNTFLD(LMXIDX)*VNS3(L0,LMX,LP,JATOM)
                       TUD = TUD + GNTFLD(LMXIDX)*VNS3(LP,LMX,L0,JATOM)
                       TDD = TDD + GNTFLD(LMXIDX)*VNS2(LP,LMX,L0,JATOM)
                       !  local orbitals part
                       IF (loor_l(l0)) THEN
                          do jlo0=1,ilo(l0,jatom)
                             TUU12(jlo0) = TUU12(jlo0) + GNTFLD(LMXIDX)*VNS1LO(l0+1,LMX,lp+1,jlo0,jatom)
                             TDU12(jlo0) = TDU12(jlo0) + GNTFLD(LMXIDX)*VNS2LO(l0+1,LMX,LP+1,jlo0,JATOM)
                          enddo
                       ENDIF
                       IF (loor_l(lp)) THEN
                          do jlop=1,ilo(lp,jatom)
                             IF (loor_l(l0)) THEN
                                do jlo0=1,ilo(l0,jatom)
                                   TUU22(jlop,jlo0) = TUU22(jlop,jlo0) + GNTFLD(LMXIDX)*VNS3LO(LP+1,LMX,L0+1,jlop,jlo0,JATOM)
                                enddo
                             ENDIF
                             TUU21(jlop) = TUU21(jlop) + GNTFLD(LMXIDX)*VNS1LO(LP+1,LMX,L0+1,jlop,JATOM)
                             TUD21(jlop) = TUD21(jlop) + GNTFLD(LMXIDX)*VNS2LO(LP+1,LMX,L0+1,jlop,JATOM)
                          enddo
                       ENDIF
                    enddo
                    call start_timer(time_hns30)
                    
                    do i=1,NV+NLO
                       jlo = i2jlo(i)
                       TMPC = TUU*conjg(ALM(I,LPMP)) + TDU*conjg(BLM(I,LPMP)) + TUU21(jlo)*conjg(CLM(I,LPMP))
                       A2R(I,LINDEX)  = A2R(I,LINDEX) + dble(TMPC)
                       A2I(I,LINDEX)  = A2I(I,LINDEX) + dimag(TMPC)
                       !
                       TMPC = TUD*conjg(ALM(I,LPMP)) + TDD*conjg(BLM(I,LPMP)) + TUD21(jlo)*conjg(CLM(I,LPMP))
                       B2R(I,LINDEX) = B2R(I,LINDEX) + dble(TMPC)
                       B2I(I,LINDEX) = B2I(I,LINDEX) + dimag(TMPC)
                    enddo
                    !
                    if (loor_l(l0)) then
                       do jlo0=1,ilo(l0,jatom)
                          do i=1,NV+NLO
                             jlo = i2jlo(i)
                             TMPC = TUU12(jlo0)*conjg(ALM(I,LPMP)) + TDU12(jlo0)*conjg(BLM(I,LPMP)) + TUU22(jlo0,jlo)*conjg(CLM(I,LPMP))
                             C2R(I,LINDEX,jlo0) = C2R(I,LINDEX,jlo0) + dble(TMPC)
                             C2I(I,LINDEX,jlo0) = C2I(I,LINDEX,jlo0) + dimag(TMPC)
                          enddo
                       enddo
                    endif
                    call stop_timer(time_hns30)
                    tot_time_vector2=tot_time_vector2+READ_CPU_TIME(time_hns30)
                 enddo ! mpidx
              enddo    ! lpidx
           enddo       ! m0idx
        enddo          ! l0idx
        
        call stop_timer(time_pwlo)
        tot_time_vector=tot_time_vector+READ_CPU_TIME(time_pwlo)
        call start_timer(time_lolo)

#ifdef _REAL
        CALL DSYR2M('L','N','T',NV+NLO,LINDEX,ONE,A2R,HSROWS,A1R,HSROWS,ONE, HS,HSROWS)
        CALL DSYR2M('L','N','T',NV+NLO,LINDEX,-ONE,A2I,HSROWS,A1I,HSROWS,ONE, HS,HSROWS)
        CALL DSYR2M('L','N','T',NV+NLO,LINDEX,ONE,B2R,HSROWS,B1R,HSROWS,ONE, HS,HSROWS)
        CALL DSYR2M('L','N','T',NV+NLO,LINDEX,-ONE,B2I,HSROWS,B1I,HSROWS,ONE, HS,HSROWS)
        IF(NLO.GT.0) THEN
           do jlo=1,maxval(ilo(:,jatom))                     
              CALL DSYR2M('L','N','T',NLO,LINDEX,ONE,C2R(NV+1,1,jlo),HSROWS,C1R(NV+1,1,jlo),HSROWS,ONE, HS(NV+1,NV+1),HSROWS)
              CALL DSYR2M('L','N','T',NLO,LINDEX,-ONE,C2I(NV+1,1,jlo),HSROWS,C1I(NV+1,1,jlo),HSROWS,ONE, HS(NV+1,NV+1),HSROWS)
           enddo
        END IF
#else        
        A1(1:NV+NLO,1:LINDEX)=A1R(1:NV+NLO,1:LINDEX)+CIMAG*A1I(1:NV+NLO,1:LINDEX)
        A2(1:NV+NLO,1:LINDEX)=A2R(1:NV+NLO,1:LINDEX)+CIMAG*A2I(1:NV+NLO,1:LINDEX)
        !               Note: not 'C', but 'T' for second matrix
        CALL prod_zher('L','N','T',NV+NLO,LINDEX,CONE,A2,HSROWS,A1,HSROWS,CONE, HS,HSROWS)
        !                  Reuse A1 and A2 for B1 and B2
        A1(1:NV+NLO,1:LINDEX)=B1R(1:NV+NLO,1:LINDEX)+CIMAG*B1I(1:NV+NLO,1:LINDEX)
        A2(1:NV+NLO,1:LINDEX)=B2R(1:NV+NLO,1:LINDEX)+CIMAG*B2I(1:NV+NLO,1:LINDEX)
        CALL prod_zher('L','N','T',NV+NLO,LINDEX,CONE,A2,HSROWS,A1,HSROWS,CONE, HS,HSROWS)
        !                  Reuse A1 and A2 for C1 and C2
        do jlo=1,maxval(ilo(:,jatom))
           A1(NV+1:NV+NLO,1:LINDEX)=C1R(NV+1:NV+NLO,1:LINDEX,jlo)+CIMAG*C1I(NV+1:NV+NLO,1:LINDEX,jlo)
           A2(NV+1:NV+NLO,1:LINDEX)=C2R(NV+1:NV+NLO,1:LINDEX,jlo)+CIMAG*C2I(NV+1:NV+NLO,1:LINDEX,jlo)
           IF(NLO.GT.0) CALL prod_zher('L','N','T',NLO,LINDEX,CONE,A2(NV+1,1),HSROWS,A1(NV+1,1),HSROWS,CONE, HS(NV+1,NV+1),HSROWS)
        enddo
#endif
        call stop_timer(time_lolo)
        tot_time_v1xv2=tot_time_v1xv2+READ_CPU_TIME(time_lolo)
     enddo mu_loop
  enddo jatom_loop
  
  if (Qprint) then
     WRITE (6,'(a,f12.1)') 'Time for alm         (hns) : ', tot_time_alm
     WRITE (6,'(a,f12.1)') 'Time for vector      (hns) : ', tot_time_vector          
     WRITE (6,'(a,f12.1)') 'Time for vector2     (hns) : ', tot_time_vector2
     WRITE (6,'(a,f12.1)') 'Time for VxV         (hns) : ', tot_time_v1xv2
     WRITE (6,'(a,f12.1)') 'Wall Time for VxV    (hns) : ', READ_WALL_TIME(time_lolo)
  endif
  
  deallocate( A1R, A1I, A2R, A2I )
  deallocate( B1R, B1I, B2R, B2I )
  deallocate( C1R, C1I, C2R, C2I )
  deallocate( ALM, BLM, CLM )
  deallocate( YL )
  deallocate(i2jlo)
#ifndef _REAL  
  deallocate( A1 )
  deallocate( A2 )
#endif
  call end_albl('HNS')
  deallocate(TUU21,TUD21,TUU12,TDU12,TUU22)
  deallocate( iind )
  CALL STOP_TIMER(time_hns)
  RETURN
END SUBROUTINE HNS
