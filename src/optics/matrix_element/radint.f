SUBROUTINE radint(JATOM,is)
  use potnlc
  use intu
  ! The major and the minor components of the Dirac solution are u=(g/r, f/r).
  !   It turns out that f \approx (dg/dr-g/r)/2 = (du/dr)*r/2 and hence du/dr = 2*f/r
  !
  !   It then follows: <u_1|d/dr|u_2> = Int[ g_1/r  2/r f_2 r^2  dr] = 2 Int[ g_1 f_2 dr]   == 2<g_1|f_2>
  !                    <u_1|1/r |u_2> = Int[ g_1/r  1/r g_2/r r^2 dr] = Int[ g_1 g_2 /r dr] == <g_1/sqrt(r)|g_2/sqrt(r)>
  !
  !  The following components are used
  !     u   <- g   | up  <- dotg
  !     du  <- f   | dup <- dotf
  !     ul  <- go  for localized orbital  
  !     dul <- fo  for localized orbital
  !ad  calculation of radial integrals:
  !ad  U   = u
  !ad  DU  = du/dr
  !ad  UP  = du/dE
  !ad  DUP = dUP/dr
  !ad 
  !ad  for LO's: 
  !ad  UL  = u_lo
  !ad  DUL = du_lo/dr
  !ad
  !ad  the radial functions are calculated in ATPAR and LOMAIN
  !ad
  INCLUDE 'param.inc'
  !
  IMPLICIT REAL*8 (A-H,O-Z)
  LOGICAL  lapw(0:lmax2),loor(0:lomax),lloor(0:lmax2) 
  !ad
  COMMON /XA/     R(NRAD),BK(3)
  COMMON /RADFU/  U(NRAD,LMAX1),UP(NRAD,LMAX1),DU(NRAD,LMAX1),DUP(NRAD,LMAX1)    
  COMMON /UHELP/  UA(NRAD), UB(NRAD), UAB(NRAD), UBA(NRAD), UAA(NRAD), UBB(NRAD)
  common /loabc/  alo(0:lomax,nloat),blo(0:lomax,nloat),clo(0:lomax,nloat),elo(0:lomax,nloat),plo(0:lomax),dplo(0:lomax),pelo(0:lomax),dpelo(0:lomax),peilo(0:lomax),pi12lo(0:lomax),pe12lo(0:lomax),UL(nrad,lomax+1),DUL(nrad,lomax+1)
  common /lolog/  nlo,nlov,nlon,lapw,ilo(0:lomax),loor,lloor
!ad..................calculation of radial integrals....................
  TWO=2.d0
  IMAX=JRI(JATOM)
  DO I=1,IMAX ! 2
     R(I)=Rnot(JATOM)*EXP((I-1)*DX(JATOM))
  ENDDO       ! 2

!!! BRISI
  !if (is.eq.1) then
  !   open(991,FILE='radials_up.dat')
  !else
  !   open(991,FILE='radials_dn.dat')
  !endif
  !do i=1,imax
  !   WRITE(991,'(7f12.7)') R(i), U(i,1), U(i,2), U(i,3), dU(i,1), dU(i,2), dU(i,3)
  !enddo
  !close(991)
!!! BRISI

  

  
  !ad
  !ad..................RUU(L)..RUUP(L)..RUPU(L)..RUPUP(L).................
  !ad
  DO L=1,LMAX     !103
     DO I=1,IMAX  !105
        SR=SQRT(R(I))
        UA(I) =U(I,L)/SR       ! g_l/sqrt(r)
        UB(I) =U(I,L+1)/SR     ! g_{l+1}/sqrt(r)
        UAB(I)=UP(I,L)/SR      ! dotg_l/sqrt(r)
        UBA(I)=UP(I,L+1)/SR    ! dotg_{l+1}/sqrt(r)
     ENDDO
     ! Kristjan : This can not be correct. I think minor and major component are here taken to be the same
     !            I think one should get correct relativistic components from previous routines.
     !            and set .FALSE. to .TRUE. Right now minor components are neglected.
      CALL RINT13(.FALSE.,UA,UA,UB,UB,RUU(L,JATOM,is),JATOM)      ! ruu  =        <u_l|1/r|u_{l+1}>
      CALL RINT13(.FALSE.,UB,UB,UAB,UAB,RUUP(L,JATOM,is),JATOM)   ! ruup =    <u_{l+1}|1/r|dotu_l>
      CALL RINT13(.FALSE.,UBA,UBA,UA,UA,RUPU(L,JATOM,is),JATOM)   ! rupu = <dotu_{l+1}|1/r|u_l>
      CALL RINT13(.FALSE.,UBA,UBA,UAB,UAB,RUPUP(L,JATOM,is),JATOM)! rupup= <dotu_{l+1}|1/r|dotu_l>
   ENDDO         !103
   !
   !ad
   !ad...............DUU1(L)..DUUP1(L)..DUPU1(L)..DUPUP1(L)................
   !ad
   DO L=1,LMAX       ! 203
      DO I=1,IMAX    ! 205
         UA(I) =U(I,L+1)       ! g_{l+1}
         UB(I) =DU(I,L)*TWO    ! 2*f_{l}
         UAB(I)=DUP(I,L)*TWO   ! 2*dotf_l
         UBA(I)=UP(I,L+1)      ! dotg_{l+1}
      ENDDO          ! 205
      CALL RINT13(.FALSE.,UA,UA,UB,UB,DUU1(L,JATOM,is),JATOM)      ! duu1  = <g_{l+1}|2*f_l>       ==    <u_{l+1}|d/dr|u_l>
      CALL RINT13(.FALSE.,UA,UA,UAB,UAB,DUUP1(L,JATOM,is),JATOM)   ! duup1 = <g_{l+1}|2*dotf_l>    ==    <u_{l+1}|d/dr|dotu_l>
      CALL RINT13(.FALSE.,UBA,UBA,UB,UB,DUPU1(L,JATOM,is),JATOM)   ! dupu1 = <dotg_{l+1}|2*f_l>    == <dotu_{l+1}|d/dr|u_l>
      CALL RINT13(.FALSE.,UBA,UBA,UAB,UAB,DUPUP1(L,JATOM,is),JATOM)! dupup1= <dotg_{l+1}|2*dotf_l> == <dotu_{l+1}|d/dr|dotu_l>
   ENDDO             ! 203
   !ad...............DUU2(L)..DUUP2(L)..DUPU2(L)..DUPUP2(L)................
   DO L=1,LMAX             ! 303
      DO I=1,IMAX          ! 305
         UA(I) =U(I,L)               ! g_l
         UB(I) =DU(I,L+1)*TWO        ! 2*f_{l+1}
         UAB(I)=DUP(I,L+1)*TWO       ! 2*dotf_{l+1}
         UBA(I)=UP(I,L)              ! dotg_l
      ENDDO                ! 305     CONTINUE
      CALL RINT13(.FALSE.,UA,UA,UB,UB,DUU2(L,JATOM,is),JATOM)      ! duu2  = <g_l|2*f_{l+1}>       == <u_l|d/dr|u_{l+1}>
      CALL RINT13(.FALSE.,UA,UA,UAB,UAB,DUUP2(L,JATOM,is),JATOM)   ! duup2 = <g_l|2*dotf_{l+1}>    == <u_l|d/dr|dotu_{l+1}>
      CALL RINT13(.FALSE.,UBA,UBA,UB,UB,DUPU2(L,JATOM,is),JATOM)   ! dupu2 = <dotg_l|2*f_{l+1}>    == <dotu_l|d/dr|u_{l+1}>
      CALL RINT13(.FALSE.,UBA,UBA,UAB,UAB,DUPUP2(L,JATOM,is),JATOM)! dupup2= <dotg_l|2*dotf_{l+1}> == <dotu_l|d/dr|dotu_{l+1}>
   ENDDO                   ! 303
   !ad..............................for LO's...............................
   DO L=0,LOMAX  ! 503
      L1=L+1      ! l1 is actually the index to the l-th component of U(l), UP(l), UL(l)
      L11=L1+1    ! l11 is the index for U(l+1), UP(l+1), UL(l+1)
      IF (.not.loor(L)) THEN
         DUUL1(L1,JATOM,is) =  0
         DUPUL1(L1,JATOM,is) = 0
         DULU2(L1,JATOM,is)  = 0
         DULUP2(L1,JATOM,is) = 0
         RUUL(L1,JATOM,is)   = 0
         RUPUL(L1,JATOM,is)  = 0
      ELSE 
         DO I=1,IMAX     ! 505
            SR=SQRT(R(I))
            UA(I) =U (I,L11)/SR     ! g_{l+1}/sqrt(r)
            UB(I) =UP(I,L11)/SR     ! dotg_{l+1}/sqrt(r)
            UAB(I)=UL(I,L1) /SR     ! go_l/sqrt(r)
         ENDDO           ! 505
         !ad............DUUL1..DUPUL1..DULU2..DULUP2..RUUL..RUPUL................
         CALL RINT13(.FALSE.,U(1,L11),U(1,L11),DUL(1,L1),DUL(1,L1),DUUL1(L1,JATOM,is),JATOM)   ! <g_{1+1}|fo_l>
         DUUL1(L1,JATOM,is)=TWO*DUUL1(L1,JATOM,is)                                             ! duul1 = 2*<g_{1+1}|fo_l>     == <u_{l+1}|d/dr|ulo_l>
         CALL RINT13(.FALSE.,UP(1,L11),UP(1,L11),DUL(1,L1),DUL(1,L1),DUPUL1(L1,JATOM,is),JATOM)! <dotg_{l+1}|fo_l>
         DUPUL1(L1,JATOM,is)=TWO*DUPUL1(L1,JATOM,is)                                           ! dupul1 = 2*<dotg_{l+1}|fo_l> == <dotu_{l+1}|d/dr|ulo_l>
         CALL RINT13(.FALSE.,UL(1,L1),UL(1,L1),DU(1,L11),DU(1,L11),DULU2(L1,JATOM,is),JATOM)   ! <go_l|f_{l+1}>
         DULU2(L1,JATOM,is)=TWO*DULU2(L1,JATOM,is)                                             ! dulu2 = 2*<go_l|f_{l+1}>     == <ulo_l|d/dr|u_{l+1}>
         CALL RINT13(.FALSE.,UL(1,L1),UL(1,L1),DUP(1,L11),DUP(1,L11),DULUP2(L1,JATOM,is),JATOM)! <go_l|dotf_{l+1}> 
         DULUP2(L1,JATOM,is)=TWO*DULUP2(L1,JATOM,is)                                           ! dulup2 = 2*<go_l|dotf_{l+1}> == <ulo_l|d/dr|dotu_{l+1}>
         CALL RINT13(.FALSE.,UA,UA,UAB,UAB,RUUL(L1,JATOM,is),JATOM)                            ! ruul = <g_{l+1}|1/r|go_l>    == <u_{l+1}|1/r|ulo_l>
         CALL RINT13(.FALSE.,UB,UB,UAB,UAB,RUPUL(L1,JATOM,is),JATOM)                           ! rupul= <dotg_{l+1}|1/r|go_l> == <dotu_{l+1}|1/r|ulo_l>
      END IF
      IF (L+1.le.LOMAX) THEN
         IF (.not.loor(L+1)) THEN
            DULU1(L1,JATOM,is)  = 0
            DULUP1(L1,JATOM,is) = 0
            DUUL2(L1,JATOM,is)  = 0
            DUPUL2(L1,JATOM,is) = 0
            DULUL1(L1,JATOM,is) = 0
            DULUL2(L1,JATOM,is) = 0
            RULU(L1,JATOM,is)   = 0
            RULUP(L1,JATOM,is)  = 0
            RULUL(L1,JATOM,is)  = 0
         ELSE
            DO I=1,IMAX    ! 506
               SR=SQRT(R(I))
               UA(I)=U(I,L1)/SR    ! g_l/sqrt(r)  
               UB(I)=UP(I,L1)/SR   ! dotg_l/sqrt(r)
               UAB(I)=UL(I,L11)/SR ! go_{l+1}/sqrt(r)
               UBA(I)=UL(I,L1)/SR  ! go_{l}/sqrt(r)
            ENDDO          ! 506
            !ad
            !ad..........DULU1..DULUP1..DUUL2..DUPUL2..DULUP2..DULUL1..DULUL2.......
            !ad
            CALL RINT13(.FALSE.,UL(1,L11),UL(1,L11),DU(1,L1),DU(1,L1),DULU1(L1,JATOM,is),JATOM)    ! <go_{l+1}|f_l>
            DULU1(L1,JATOM,is)=TWO*DULU1(L1,JATOM,is)                                              ! dulu1 = 2*<go_{l+1}|f_l>      == <ulo_{l+1}|d/dr|u_l>
            CALL RINT13(.FALSE.,UL(1,L11),UL(1,L11),DUP(1,L1),DUP(1,L1),DULUP1(L1,JATOM,is),JATOM) ! <go_{l+1}|dotf_l>
            DULUP1(L1,JATOM,is)=TWO*DULUP1(L1,JATOM,is)                                            ! dulup1 = 2*<go_{l+1}|dotf_l>  == <ulo_{l+1}|d/dr|dotu_l>
            CALL RINT13(.FALSE.,U(1,L1),U(1,L1),DUL(1,L11),DUL(1,L11),DUUL2(L1,JATOM,is),JATOM)    ! <g_l|fo_{l+1}>
            DUUL2(L1,JATOM,is)=TWO*DUUL2(L1,JATOM,is)                                              ! duul2 = 2*<g_l|fo_{l+1}>      == <u_l|d/dr|ulo_{l+1}>
            CALL RINT13(.FALSE.,UP(1,L1),UP(1,L1),DUL(1,L11),DUL(1,L11),DUPUL2(L1,JATOM,is),JATOM) ! <dotg|fo_{l+1}>
            DUPUL2(L1,JATOM,is)=TWO*DUPUL2(L1,JATOM,is)                                            ! dupul2 = 2*<dotg_l|fo_{l+1}>  == <dotu_l|d/dr|ulo_{l+1}>
            CALL RINT13(.FALSE.,UL(1,L11),UL(1,L11),DUL(1,L1),DUL(1,L1),DULUL1(L1,JATOM,is),JATOM) ! <go_{l+1}|fo_l>
            DULUL1(L1,JATOM,is)=TWO*DULUL1(L1,JATOM,is)                                            ! dulul1 = 2*<go_{l+1}|fo_l>    == <ulo_{l+1}|d/dr|ulo_l>
            CALL RINT13(.FALSE.,UL(1,L1),UL(1,L1),DUL(1,L11),DUL(1,L11),DULUL2(L1,JATOM,is),JATOM) ! <go_l|fo_{l+1}>
            DULUL2(L1,JATOM,is)=TWO*DULUL2(L1,JATOM,is)                                            ! dulul2 = 2*<go_l|fo_{l+1}>== <ulo_l|d/dr|ulo_{l+1}>
            !ad...........................RULU..RULUP..RULUL........................
            CALL RINT13(.FALSE.,UAB,UAB,UA,UA,RULU(L1,JATOM,is),JATOM)      ! rulu = <go_{l+1}|1/r|g_l>    = <ulo_{l+1}|1/r|u_l>
            CALL RINT13(.FALSE.,UAB,UAB,UB,UB,RULUP(L1,JATOM,is),JATOM)     ! rulup = <go_{l+1}|1/r|dotg_l>= <ulo_{l+1}|1/r|dotu_l>
            CALL RINT13(.FALSE.,UAB,UAB,UBA,UBA,RULUL(L1,JATOM,is),JATOM)   ! rulul = <go_{l+1}|1/r|go_l>  = <ulo_{l+1}|1/r|ulo_l>
         END IF
      END IF
   ENDDO    ! 503
   
   !
   !....... END CALCULATION OF RADIAL INTEGRALS ...............
   RETURN 
 END SUBROUTINE radint
