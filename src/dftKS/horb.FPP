#ifdef _REAL
#define tpcast(X) dble(X)  !! typecasting in fortran
#define Ttype REAL*8
#else
#define tpcast(X) X
#define Ttype COMPLEX*16
#endif

SUBROUTINE  Horb(NV,NAT,POS,Vol,R,MULT,PEI,ROTLOC,ROTIJ)
  use matrices, only: HSROWS, KZZ,HS, Kn
  use loabc, only : ALO, BLO, CLO, PI12LO, PE12LO,pilolo
  use lolog, only : NLO, ILO
  use orb, only   : VORB, NATORB, IATOM, NLORB, LORB
  use structure, only: ndf
  use albl, only : AL, BL,init_albl, end_albl
  use lapw_timer, only : READ_CPU_TIME, READ_WALL_TIME, time_pwpw, time_pwlo,time_lolo,  START_TIMER, STOP_TIMER, time_horb_setup, time_horb
  use param, only: blocksize, lmax, nslmax
  use mpi, only : Qprint
  IMPLICIT NONE
  !        Arguments
  INTEGER, intent(in) :: NV, NAT!, NNSL
  REAL*8,  intent(in) :: POS(3,NDF)
  REAL*8,  intent(in) :: Vol, R(NAT)
  INTEGER, intent(in) :: MULT(NAT)!, LQIND(NAT), LQNS(6,NGAU,NAT)
  REAL*8,  intent(in) :: PEI(LMAX,NAT)
  REAL*8,  intent(in) :: ROTLOC(3,3,NAT), ROTIJ(3,3,NDF)
  ! locals
  COMPLEX*16 :: TMP1,TMP2,TMP3, CIL
  Ttype      :: TMP4,TMP5,TMP6,TMP7,TMP8
  COMPLEX*16 :: pot(-3:3,-3:3)
  COMPLEX*16, allocatable :: alm(:,:,:),blm(:,:,:),clm(:,:,:), almt(:,:,:),blmt(:,:,:)
  COMPLEX*16, allocatable :: YL(:,:,:), PHS(:,:)
  !INTEGER,    allocatable :: i_end(:)
  !        Locals
  REAL*8     :: RHOATM, PI, TPI, ARG
  REAL*8     :: ROTV1(3), ROTV2(3)
  INTEGER    :: i1, i2, jlo, i12, i22, j, jlo2, mu, nlv, iat, indatm, jatom, lvor, l0m0, l0, m0, m1, m2, i
  INTEGER    :: iset_end,npe_iset,myid_iset,ldm,iset,INDATM0,ihelp,iouter,imin,imax
  REAL*8     :: tot_time_horb_setup, tot_time_pwpw, tot_time_pwlo, tot_time_lolo
  !        External Subroutines
  EXTERNAL           YLM
  !        Local Parameters
  COMPLEX*16  :: CIMAG, CZERO
  PARAMETER  (CIMAG  = (0.0D+0,1.0D+0), CZERO=(0.0D0,0.0D0) )
  !
  CALL START_TIMER(time_horb)
  !
  tot_time_horb_setup=0.0d0
  tot_time_pwpw=0.0d0
  tot_time_pwlo=0.0d0
  tot_time_lolo=0.0d0
  ldm=HSROWS
  iset_end=1
  allocate( ALM(-3:3,ldm,iset_end) )
  allocate( BLM(-3:3,ldm,iset_end) )
  allocate( CLM(-3:3,ldm,iset_end) )
  allocate( ALMt(-3:3,ldm,iset_end) )
  allocate( BLMt(-3:3,ldm,iset_end) )
  allocate( YL(0:NSLMAX*NSLMAX-1,ldm,iset_end),PHS(ldm,iset_end) )

  lvor=0
  do iat= 1,natorb
     do nlv=1,nlorb(iat)
        lvor=max(lvor,lorb(nlv,iat))
     enddo
  enddo
  call init_albl(lvor,'ORB')

  iset_end=1          

  PI = 4.0D+0*ATAN(1.0D+0)
  TPI= 2.0D0*PI

  j_neq_loop: DO iat = 1, NATORB
     jatom=iatom(iat)
     INDATM0 = 0 + sum(mult(1:jatom-1))
     RHOATM = 4.0D+0*PI*R(JATOM)*R(JATOM)/sqrt(Vol)
     lvor=maxval(lorb(1:nlorb(iat),iat))
     call make_albl(jatom,R(jatom),lvor,NV,'ORB',1,0)
     j_eq_loop: DO mu=1,MULT(JATOM)
        INDATM = INDATM0 + mu
        nvl_loop: do nlv=1,nlorb(iat)
           call start_timer(time_horb_setup)
           lvor=lorb(nlv,iat)
           do m1=-lvor,lvor
              do m2=-lvor,lvor
                 pot(m2,m1)=vorb(iat,lvor,m1,m2)
              enddo
           enddo
           alm(:,:,:)=czero
           blm(:,:,:)=czero
           clm(:,:,:)=czero
           PHS(:,:)=czero
           !do iset=1,1
           iset=1
           npe_iset=1
           myid_iset=0 !mycolhs
           IHELP = 0
           DO IOUTER = 0, (NV+NLO)/(NPE_iset*BLOCKSIZE)
              IMIN = (IOUTER * NPE_iset + MYID_iset) * BLOCKSIZE + 1
              IMAX = MIN(NV+NLO, (IOUTER * NPE_iset + MYID_iset + 1) * BLOCKSIZE)
              DO I = IMIN, IMAX
                 IHELP = IHELP + 1
                 !
                 !        calculate phase-factor for Alm, Blm
                 !
                 ARG = ( POS(1,INDATM)*KZZ(1,I) + POS(2,INDATM)*KZZ(2,I) + POS(3,INDATM)*KZZ(3,I))*TPI
                 PHS(Ihelp,iset) = DCMPLX( COS(ARG), SIN(ARG))
                 !
                 !        calculate spherical-harmonics for Alm, Blm
                 !
                 !CALL ROTATE (Kn(:,i),ROTIJ(1,1,INDATM),ROTV1)
                 !CALL ROTATE (ROTV1,ROTLOC(1,1,JATOM),ROTV2)
                 ROTV1 = matmul(ROTIJ(:,:,INDATM),Kn(:,i))
                 ROTV2 = matmul(ROTLOC(:,:,JATOM),ROTV1)
                 CALL YLM (ROTV2,lvor,YL(0,Ihelp,iset))
              enddo
           enddo
           !        Alm, Blm
           !
           l0m0=lvor**2
           cil=cimag**lvor
           l0=lvor
           do m0 = -l0, l0         
              ihelp = 0
              do iouter = 0, (nv)/(npe_iset*blocksize)
                 imin = (iouter * npe_iset + myid_iset) * blocksize + 1
                 imax = min(nv, (iouter * npe_iset + myid_iset + 1) * blocksize)
                 do i = imin, imax
                    ihelp = ihelp + 1
                    tmp1=rhoatm*cil*phs(ihelp,iset)*conjg(yl(l0m0,ihelp,iset))
                    alm(m0,ihelp,iset)=tmp1*al(ihelp,l0,iset)
                    blm(m0,ihelp,iset)=tmp1*bl(ihelp,l0,iset)
                 enddo
              enddo

              do jlo=1,ilo(lvor,jatom)
                 call lohns(jatom,mult,i1,i2,l0,jlo)                            
                 ihelp = 0
                 do iouter = 0, (nv+nlo)/(npe_iset*blocksize)
                    imin = (iouter * npe_iset + myid_iset) * blocksize + 1
                    imax = min(nv+nlo, (iouter * npe_iset + myid_iset + 1) * blocksize)
                    do i = imin, imax
                       ihelp = ihelp + 1
                       if ((i.ge.nv+i1).and.(i.le.nv+i2)) then
                          tmp1=rhoatm*cil*phs(ihelp,iset)*conjg(yl(l0m0,ihelp,iset))
                          blm(m0,ihelp,iset)=blo(l0,jlo,jatom)*tmp1
                          clm(m0,ihelp,iset)=clo(l0,jlo,jatom)*tmp1
                          alm(m0,ihelp,iset)=alo(l0,jlo,jatom)*tmp1
                       endif
                    enddo
                 enddo
              enddo
              l0m0=l0m0+1                                                  
           enddo
           !   end of alm, blm and clm calculation
           !   add orbital term to the hamiltonian
           ALMT(:,:,iset)=czero
           BLMT(:,:,iset)=czero
           ihelp=0

           do iouter = 0, (nv)/(npe_iset*blocksize)
              imin = (iouter * npe_iset + myid_iset) * blocksize + 1
              imax = min(nv, (iouter * npe_iset + myid_iset + 1) * blocksize)
              do i = imin, imax
                 ihelp = ihelp + 1
                 DO M0=-LVOR,LVOR
                    TMP1=czero
                    TMP2=czero
                    do M1=-LVOR,LVOR
                       TMP1=TMP1+pot(m1,m0)*alm(m1,ihelp,iset)
                       TMP2=TMP2+pot(m1,m0)*blm(m1,ihelp,iset)
                    ENDDO
                    ALMT(M0,ihelp,iset)=TMP1
                    BLMT(M0,ihelp,iset)=TMP2*pei(lvor+1,jatom)
                 ENDDO
              ENDDO
           enddo
           !enddo

           call stop_timer(time_horb_setup)
           tot_time_horb_setup=tot_time_horb_setup+READ_CPU_TIME(time_horb_setup)

           !   <APW||APW>        
           call start_timer(time_pwpw)
           DO J=1,NV
              do I=J,NV
                 TMP4=0.0D0
                 TMP5=0.0D0
                 do m0=-lvor,lvor
                    TMP4=TMP4+tpcast(ALMT(M0,J,1)*conjg(alm(m0,i,1)))
                    TMP5=TMP5+tpcast(BLMT(M0,J,1)*conjg(blm(m0,i,1)))
                 enddo
                 HS(i,j)=HS(i,j)+TMP4+TMP5
              enddo
           enddo
           call stop_timer(time_pwpw)
           tot_time_pwpw=tot_time_pwpw+READ_CPU_TIME(time_pwpw)

           do jlo=1,ilo(lvor,jatom)
              !  <APW||LO>
              call start_timer(time_pwlo)
              call lohns(jatom,mult,i1,i2,lvor,jlo)
              do i=nv+i1,nv+i2
                 do j=1,nv
                    TMP4=0.0D0
                    TMP5=0.0D0
                    TMP6=0.0D0
                    TMP7=0.0D0
                    do m0=-lvor,lvor
                       TMP1=czero
                       TMP2=czero
                       do m1=-lvor,lvor
                          TMP1=TMP1+alm(m1,j,1)*pot(m1,m0)
                          TMP2=TMP2+blm(m1,j,1)*pot(m1,m0)
                       enddo
                       TMP4=TMP4+tpcast(TMP1*conjg(alm(m0,i,1)))
                       TMP5=TMP5+tpcast(TMP2*conjg(blm(m0,i,1)))
                       TMP6=TMP6+tpcast(TMP1*conjg(clm(m0,i,1)))
                       TMP7=TMP7+tpcast(TMP2*conjg(clm(m0,i,1)))
                    enddo
                    HS(I,J)=HS(I,J)+TMP4+TMP5*pei(lvor+1,jatom)+TMP6*pi12lo(lvor,jlo,jatom)+TMP7*pe12lo(lvor,jlo,jatom)
                 enddo
              enddo
              call stop_timer(time_pwlo)
              tot_time_pwlo=tot_time_pwlo+READ_CPU_TIME(time_pwlo)
              !  <LO||LO>                   
              do jlo2=1,jlo
                 call start_timer(time_lolo)
                 call lohns(jatom,mult,i12,i22,lvor,jlo2)                      
                 do i=nv+i1,nv+i2
                    do j=nv+i12,min(i,nv+i22)
                       TMP4=0.0D0
                       TMP5=0.0D0
                       TMP6=0.0D0
                       TMP7=0.0D0
                       TMP8=0.0d0
                       do m0=-lvor,lvor
                          TMP1=czero
                          TMP2=czero
                          TMP3=czero
                          do m1=-lvor,lvor
                             TMP1=TMP1+pot(m1,m0)*alm(m1,j,1)
                             TMP2=TMP2+pot(m1,m0)*blm(m1,j,1)
                             TMP3=TMP3+pot(m1,m0)*clm(m1,j,1)
                          enddo
                          TMP4=TMP4+tpcast(TMP1*dconjg(alm(m0,i,1))+ dconjg(clm(m0,i,1))*TMP3)
                          TMP5=TMP5+tpcast(dconjg(blm(m0,i,1))*TMP2)
                          TMP6=TMP6+tpcast(dconjg(clm(m0,i,1))*TMP1+ dconjg(alm(m0,i,1))*TMP3)
                          TMP7=TMP7+tpcast(dconjg(blm(m0,i,1))*TMP3+ dconjg(clm(m0,i,1))*TMP2)
                          TMP8=TMP8+tpcast(dconjg(clm(m0,i,1))*TMP3)
                       enddo
                       HS(I,J)=HS(I,J)+TMP5*pei(lvor+1,jatom)+TMP6*pi12lo(lvor,jlo2,jatom)+TMP7*pe12lo(lvor,jlo2,jatom)+TMP4+TMP8*pilolo(lvor,jlo2,jlo,jatom)
                    enddo
                 enddo
                 call stop_timer(time_lolo)
                 tot_time_lolo=tot_time_lolo+READ_CPU_TIME(time_lolo)
              enddo
           enddo
        enddo nvl_loop
     enddo j_eq_loop
  enddo j_neq_loop

  if (Qprint) then
     WRITE (6,'(a,f12.1)') 'Time for setup      (horb) : ', tot_time_horb_setup
     WRITE (6,'(a,f12.1)') 'Time for pwpw       (horb) : ', tot_time_pwpw
     WRITE (6,'(a,f12.1)') 'Time for pwlo       (horb) : ', tot_time_pwlo
     WRITE (6,'(a,f12.1)') 'Time for lolo       (horb) : ', tot_time_lolo
  endif
  
  deallocate( ALM )
  deallocate( BLM )
  deallocate( CLM )
  deallocate( YL )
  deallocate( PHS )
  deallocate( ALMT )
  deallocate( BLMT )
  call end_albl('ORB')
  CALL STOP_TIMER(time_horb)
  
  RETURN

END SUBROUTINE Horb
