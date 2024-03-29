SUBROUTINE MMatrix(JATOM,ivmax,NN_,N_,is)
  use mxyz
  use ablm, a=>alm,b=>blm,c=>clm
  use intu
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
  !......COMPUTING : < PSI'| GRAD | PSI >...............................
  !................. | PSI > := (ALM(NB) UL + BLM(NB) UPL) YLM .........
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  !kh Note that radint.f has computed the following matrix elements
  !  ruu   =        <u_l|1/r|u_{l+1}>
  !  ruup  =    <u_{l+1}|1/r|dotu_l>
  !  rupu  = <dotu_{l+1}|1/r|u_l>
  !  rupup = <dotu_{l+1}|1/r|dotu_l>
  !  duu1  =    <u_{l+1}|d/dr|u_l>
  !  duup1 =    <u_{l+1}|d/dr|dotu_l>
  !  dupu1 = <dotu_{l+1}|d/dr|u_l>
  !  dupup1= <dotu_{l+1}|d/dr|dotu_l>
  !  duu2  =        <u_l|d/dr|u_{l+1}>
  !  duup2 =        <u_l|d/dr|dotu_{l+1}>
  !  dupu2 =     <dotu_l|d/dr|u_{l+1}>
  !  dupup2=     <dotu_l|d/dr|dotu_{l+1}>
  !
  ! duul1  = <u_{l+1}|d/dr|ulo_l>
  ! dupul1 = <dotu_{l+1}|d/dr|ulo_l>
  ! dulu2  = <ulo_l|d/dr|u_{l+1}>
  ! dulup2 = <ulo_l|d/dr|dotu_{l+1}>
  ! ruul   = <u_{l+1}|1/r|ulo_l>
  ! rupul  = <dotu_{l+1}|1/r|ulo_l>
  ! dulu1  = <ulo_{l+1}|d/dr|u_l>
  ! dulup1 = <ulo_{l+1}|d/dr|dotu_l>
  ! duul2  = <u_l|d/dr|ulo_{l+1}>
  ! dupul2 = <dotu_l|d/dr|ulo_{l+1}>
  ! dulul1 = <ulo_{l+1}|d/dr|ulo_l>
  ! dulul2 = <ulo_l|d/dr|ulo_{l+1}>
  ! rulu   = <ulo_{l+1}|1/r|u_l>
  ! rulup  = <ulo_{l+1}|1/r|dotu_l>
  ! rulul  = <ulo_{l+1}|1/r|ulo_l>
  !
  IMPLICIT REAL*8 (A-H,O-Z)
  INCLUDE 'param.inc'
  INTEGER JATOM,J,N_(ivmax),NN_(ivmax)
  COMPLEX*16  CZERO,IMAG
  LOGICAL     lapw(0:lmax2),loor(0:lomax),lloor(0:lmax2)
  common /lolog/   nlo,nlov,nlon,lapw,ilo(0:lomax),loor,lloor
  !ad 
  DATA    CZERO/(0.0D0,0.0D0)/IMAG/(0.0D+0,1.0D+0)/ 
  DATA    ONE/1.0D+0/TWO/2.0D+0/THREE/3.0D+0/
  !........INLINE FUNCTIONS FOR YLM OTHOGONALITY ....................
  !      
  F0(L)   = (TWO*L+ONE)*(TWO*L+THREE)
  F1(L,M) = - SQRT(( L + M + ONE ) * ( L + M + TWO ) / F0(L))
  F3(L,M) =   SQRT(( L - M + ONE ) * ( L - M + TWO ) / F0(L))
  F5(L,M) =   SQRT(( L - M + ONE ) * ( L + M + ONE ) / F0(L))
  !..................................................................
  SX_(1:ivmax)=0.d0
  SY_(1:ivmax)=0.d0
  SZ_(1:ivmax)=0.d0
  J=JATOM
  IDEX=0 
  DO L=0,LMAX-1 
     L1=L+1  ! This is becuase Duu and Ruu are stored starting by index 1, rather then zero.
     L2=L+2                           ! Here we combine together the following matrix elements:
     H1  = (Duu1  (L1,J,is) - L * Ruu  (L1,J,is))  ! h1=   <u_{l+1}|d/dr - l/r  |u_l>
     H2  = (Duup1 (L1,J,is) - L * Ruup (L1,J,is))  ! h2=   <u_{l+1}|d/dr - l/r  |dotu_l>
     H3  = (Dupu1 (L1,J,is) - L * Rupu (L1,J,is))  ! h3=<dotu_{l+1}|d/dr - l/r  |u_l>
     H4  = (Dupup1(L1,J,is) - L * Rupup(L1,J,is))  ! h4=<dotu_{l+1}|d/dr - l/r  |dotu_l>
     H5  = (Duu2  (L1,J,is) + L2* Ruu  (L1,J,is))  ! h5=       <u_l|d/dr+(l+2)/r|u_{l+1}>
     H6  = (Duup2 (L1,J,is) + L2* Rupu (L1,J,is))  ! h6=       <u_l|d/dr+(l+2)/r|dotu_{l+1}>
     H7  = (Dupu2 (L1,J,is) + L2* Ruup (L1,J,is))  ! h7=    <dotu_l|d/dr+(l+2)/r|u_{l+1}>
     H8  = (Dupup2(L1,J,is) + L2* Rupup(L1,J,is))  ! h8=    <dotu_l|d/dr+(l+2)/r|dotu_{l+1}>
     !ad   for LO's
     if (L1.LE.LOMAX+1) THEN
        HL1 = (Duul1 (L1,J,is) - L * Ruul(L1,J,is))  ! <u_{l+1}|d/dr-l/r|ulo_l>
        HL2 = (Dupul1(L1,J,is) - L * Rupul(L1,J,is)) ! <dotu_{l+1}|d/dr-l/r|ulo_l>
        HL3 = (Dulu2 (L1,J,is) + L2* Ruul(L1,J,is))  ! <ulo_l|d/dr+(l+2)/r|u_{l+1}>
        HL4 = (Dulup2(L1,J,is) + L2* Rupul(L1,J,is)) ! <ulo_l|d/dr+(l+2)/r|dotu_{l+1}>
        HL5 = (Dulu1 (L1,J,is) - L * Rulu(L1,J,is))  ! <ulo_{l+1}|d/dr-l/r|u_l>
        HL6 = (Dulup1(L1,J,is) - L * Rulup(L1,J,is)) ! <ulo_{l+1}|d/dr-l/r|dotu_l>
        HL7 = (Duul2 (L1,J,is) + L2* Rulu(L1,J,is))  ! <u_l|d/dr+(l+2)/r|ulo_{l+1}>
        HL8 = (Dupul2(L1,J,is) + L2* Rulup(L1,J,is)) ! <dotu_l|d/dr+(l+2)/r|ulo_{l+1}>
        HL9 = (Dulul1(L1,J,is) - L * Rulul(L1,J,is)) ! <ulo_{l+1}|d/dr-l/r|ulo_l>
        HL0 = (Dulul2(L1,J,is) + L2* Rulul(L1,J,is)) ! <ulo_l|d/dr+(l+2)/r|ulo_{l+1}>
     end if
     !ad..........................................................M=-L,L
     DO M=-L,L 
        !ad
        !ad......IDEX  = (l,m)    -> l**2+l+m+1=l*(l+1)+m+1
        !ad......IDEX1 = (l+1,m)
        !ad......IDEX11=(l+1,m+1)
        !ad......IDEX10=(l+1,m-1)
        !ad
        !...IDEX->l,m..IDEX11->l+1,m+1..IDEX10->l+1,m-1..IDEX1->l+1,m
        IDEX  = IDEX+1
        IDEX11= IDEX+2*L+3 
        IDEX10= IDEX+2*L+1 
        IDEX1 = IDEX+2*L+2
        !ad   for LAPW's
        do iv=1,ivmax
           n=n_(iv)
           nn=nn_(iv)
           !ad   x+iy component
           ! S(x+iy) <= f1(l,m) * <a(j,l+1,m+1)u_{l+1}+b(j,l+1,m+1)dotu_{l+1}|d/dr  -  l/r|a(i,l,  m  )u_l    +b(i,l,  m  )dotu_l> +
           !            f3(l,m) * <a(j,l,  m  )u_l    +b(j,l,  m  )dotu_l    |d/dr+(l+2)/r|a(i,l+1,m-1)u_{l+1}+b(i,l+1,m-1)dotu_{l+1}>
           SX_(iv)=SX_(iv)+(H1*CONJG(A(nn,IDEX11))*A(n,IDEX) + H2*CONJG(A(nn,IDEX11))*B(n,IDEX)  + H3*CONJG(B(nn,IDEX11))*A(n,IDEX) + H4*CONJG(B(nn,IDEX11))*B(n,IDEX))*F1(L,M)  &
                + (H5*CONJG(A(nn,IDEX))  *A(n,IDEX10)+ H6*CONJG(A(nn,IDEX))  *B(n,IDEX10) + H7*CONJG(B(nn,IDEX))  *A(n,IDEX10) + H8*CONJG(B(nn,IDEX))  *B(n,IDEX10) )*F3(L,M)
           !ad   x-iy component
           ! S(x-iy) <= f3(l,m) * <a(j,l+1,m-1)u_{l+1}+b(j,l+1,m-1)dotu_{l+1}|d/dr  -  l/r|a(i,l,m)u_l+b(i,l,m)dotu_l>+
           !            f1(l,m) * <a(j,l,  m  )u_l    +b(j,l  ,m  )dotu_l    |d/dr+(l+2)/r|a(i,l+1,m+1)u_{l+1}+b(i,l+1,m+1)dotu_{l+1}>
           SY_(iv)=SY_(iv) +(H1*CONJG(A(nn,IDEX10))*A(n,IDEX) + H2*CONJG(A(nn,IDEX10))*B(n,IDEX)+ H3*CONJG(B(nn,IDEX10))*A(n,IDEX)+ H4*CONJG(B(nn,IDEX10))*B(n,IDEX))*F3(L,M) &
                +(H5*CONJG(A(nn,IDEX))  *A(n,IDEX11) + H6*CONJG(A(nn,IDEX))  *B(n,IDEX11) + H7*CONJG(B(nn,IDEX))  *A(n,IDEX11) + H8*CONJG(B(nn,IDEX))  *B(n,IDEX11))*F1(L,M)
           !ad   z component
           ! S(z)   <= f5(l,m) * <a(j,l+1,m)u_{l+1}+b(j,l+1,m)dotu_{l+1}|d/dr - l/r  |a(i,l,  m)u_l    +b(i,l,  m)dotu_l> +
           !           f5(l,m) * <a(j,l,  m)u_l    +b(j,l,  m)dotu_l    |d/dr+(l+2)/r|a(i,l+1,m)u_{l+1}+b(i,l+1,m)dotu_{l+1}>
           SZ_(iv)=SZ_(iv) +(H1*CONJG(A(nn,IDEX1)) *A(n,IDEX) + H2*CONJG(A(nn,IDEX1)) *B(n,IDEX)+ H3*CONJG(B(nn,IDEX1)) *A(n,IDEX) + H4*CONJG(B(nn,IDEX1)) *B(n,IDEX) &
                + H5*CONJG(A(nn,IDEX))  *A(n,IDEX1) + H6*CONJG(A(nn,IDEX))  *B(n,IDEX1) + H7*CONJG(B(nn,IDEX))  *A(n,IDEX1) + H8*CONJG(B(nn,IDEX))  *B(n,IDEX1))*F5(L,M)
        enddo
        !ad   for LO's
        IF (lloor(l)) THEN
           do iv=1,ivmax
              n=n_(iv)
              nn=nn_(iv)
              !ad   x+iy component
              ! S(x+iy)  = f1(l,m)*<a(j,l+1,m+1)u_{l+1}+b(j,l+1,m+1)dotu_{l+1}|d/dr  - l /r|c(i,l,m)ulo_l>
              ! S(x+iy) += f3(l,m)*<c(j,l,m)ulo_l                             |d/dr+(l+2)/r|a(i,l+1,m-1)u_{l+1}+b(i,l+1,m-1)dotu_{l+1}>
              SX_(iv)=SX_(iv)+(HL1*CONJG(A(nn,IDEX11))*C(n,IDEX)+HL2*CONJG(B(nn,IDEX11))*C(n,IDEX))*F1(L,M) + (HL3*CONJG(C(nn,IDEX))*A(n,IDEX10)+ HL4*CONJG(C(nn,IDEX))*B(n,IDEX10))*F3(L,M)
              !ad   x-iy component
              ! S(x-iy) = f3(l,m)*<a(j,l+1,m-1)u_{l+1}+b(j,l+1,m-1)dotu_{l+1}|d/dr  - l /r|c(i,l,m)ulo_l>
              ! S(x-iy)+= f1(l,m)*<c(j,l,m)ulo_l                             |d/dr+(l+2)/r|a(i,l+1,m+1)u_{l+1}+b(i,l+1,m+1)dotu_{l+1}>
              SY_(iv)=SY_(iv)+(HL1*CONJG(A(nn,IDEX10))*C(n,IDEX)+HL2*CONJG(B(nn,IDEX10))*C(n,IDEX))*F3(L,M) + (HL3*CONJG(C(nn,IDEX))*  A(n,IDEX11)+ HL4*CONJG(C(nn,IDEX))*  B(n,IDEX11))*F1(L,M)
              !ad   z component
              ! S(z) =  f5(l,m)*<a(j,l+1,m)u_{l+1}+b(j,l+1,m)dotu_{l+1}|d/dr  -l  /r|c(i,l,m)ulo_l>
              ! S(z) += f5(l,m)*<c(j,l,m)ulo_l                         |d/dr+(l+2)/r|a(i,l+1,m)u_{l+1}+b(i,l+1,m)dotu_{l+1}>
              SZ_(iv)=SZ_(iv)+(HL1*CONJG(A(nn,IDEX1)) *C(n,IDEX)+HL2*CONJG(B(nn,IDEX1))*C(n,IDEX)+HL3*CONJG(C(nn,IDEX))*A(n,IDEX1)+HL4*CONJG(C(nn,IDEX))*B(n,IDEX1))*F5(L,M)
           enddo
        END IF
        IF (lloor(l+1)) THEN   
           do iv=1,ivmax
              n=n_(iv)
              nn=nn_(iv)
              !ad   x+iy component
              SX_(iv)=SX_(iv) +(HL5*CONJG(C(nn,IDEX11))*A(n,IDEX) + HL6*CONJG(C(nn,IDEX11))*B(n,IDEX))*F1(L,M)+(HL7*CONJG(A(nn,IDEX))  *C(n,IDEX10)+ HL8*CONJG(B(nn,IDEX))*C(n,IDEX10))*F3(L,M)
              !ad   x-iy component
              SY_(iv)=SY_(iv)+(HL5*CONJG(C(nn,IDEX10))*A(n,IDEX)+HL6*CONJG(C(nn,IDEX10))*B(n,IDEX))*F3(L,M)+(HL7*CONJG(A(nn,IDEX))  *C(n,IDEX11)+ HL8*CONJG(B(nn,IDEX))  *C(n,IDEX11))*F1(L,M)
              !ad   z component
              SZ_(iv)=SZ_(iv)+(HL5*CONJG(C(nn,IDEX1)) *A(n,IDEX)+HL6*CONJG(C(nn,IDEX1)) *B(n,IDEX)+HL7*CONJG(A(nn,IDEX))  *C(n,IDEX1)+HL8*CONJG(B(nn,IDEX))*C(n,IDEX1))*F5(L,M)
           enddo
        END IF
        IF ( (lloor(l)).AND.(lloor(l+1)) ) THEN
           do iv=1,ivmax
              n=n_(iv)
              nn=nn_(iv)
              !ad   x+iy component
              ! S(x+iy)  = f1(l,m) * <c(j,l+1,m+1)ulo_{l+1}|d/dr-l/r|c(i,l,m)ulo_l>
              ! S(x+iy) += f3(l,m) * <c(j,l,m)ulo_l|d/dr+(l+2)/r|c(i,l+1,m-1)ulo_{l+1}>
              SX_(iv)=SX_(iv)+(HL9*CONJG(C(nn,IDEX11))*C(n,IDEX))*F1(L,M)+(HL0*CONJG(C(nn,IDEX))  *C(n,IDEX10))*F3(L,M)
              !ad   x-iy component
              SY_(iv)=SY_(iv)+(HL9*CONJG(C(nn,IDEX10))*C(n,IDEX))*F3(L,M)+(HL0*CONJG(C(nn,IDEX))  *C(n,IDEX11))*F1(L,M)
              !ad   z component
              SZ_(iv)=SZ_(iv)+(HL9*CONJG(C(nn,IDEX1)) *C(n,IDEX)+HL0*CONJG(C(nn,IDEX))  *C(n,IDEX1))*F5(L,M)
           enddo
        END IF
        !ad    LO's done
     END DO
     !ad................................................end.......M=-L,L
  END DO
  !ad................................................end...L=0,LMAX-1
  !ad    Transformation to Cartesian coordinates
  do iv=1,ivmax
     MX_(iv)=(SX_(iv)+SY_(iv)) / TWO
     MY_(iv)=(SX_(iv)-SY_(iv)) /IMAG  / TWO
     MZ_(iv)=SZ_(iv)
  enddo
  
  RETURN
END SUBROUTINE MMatrix
