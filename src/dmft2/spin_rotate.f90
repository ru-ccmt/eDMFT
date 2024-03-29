! @Copyright 2007 Kristjan Haule
! 

SUBROUTINE GetSpinRotation(Rspin,rotij,crotloc,norbitals,natom,natm,iso,lmaxp,iatom,nl,cix,ll,iorbital)
  USE com_mpi, ONLY: Qprint
  USE structure, ONLY: BR1, rot_spin_quantization
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: Rspin(2,2,norbitals)
  REAL*8, intent(in)  :: rotij(3,3,natm), crotloc(3,3,natom)
  INTEGER, intent(in) :: norbitals, natom, natm, iso, lmaxp
  INTEGER, intent(in) :: iatom(natom), nl(natom), cix(natom,4), ll(natom,4), iorbital(natom,lmaxp+1)
  ! locals
  INTEGER :: icase, latom, lcase, icix, l, nind, iorb, i, j
  REAL*8  :: Det, phi1, the1, psi1
  COMPLEX*16 :: Rtmp(2,2)
  REAL*8 :: tmp3(3,3), Trans3(3,3), BR1inv(3,3), rotij_cartesian(3,3), Id(3,3), crotloc_x_rotij_cartesian(3,3)
  ! external
  REAL*8 :: detx
  INTEGER            STDERR
  PARAMETER          (STDERR = 99)
  !
  call inv_3x3(BR1,BR1inv)
  Id(:,:)=0.d0
  Id(1,1)=1.d0
  Id(2,2)=1.d0
  Id(3,3)=1.d0
  DO icase=1,natom
     latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
     do lcase=1,nl(icase)
        icix = cix(icase,lcase)
        if ( icix.EQ.0 ) CYCLE
        l = ll(icase,lcase)
        nind = (2*l+1)*iso
        iorb = iorbital(icase,lcase)
        !!  local_axis_defined_by_locrot  <- local_axis_of_equivalent_atom <- from_spin_quantization_to_global_cartesian
        !!* Trans3 = crotloc(:,:,icase) * rotij_cartesian * rot_spin_quantization
        tmp3 = matmul(BR1,rotij(:,:,latom))
        rotij_cartesian = matmul(tmp3,BR1inv)
        
        crotloc_x_rotij_cartesian = matmul(crotloc(:,:,icase),rotij_cartesian)
        Trans3 = matmul(crotloc_x_rotij_cartesian, rot_spin_quantization)
        !
        Det = detx(Trans3)
        Trans3 = transpose(Trans3*Det)
        if (sum(abs(Trans3-Id)).LT.1e-10) Trans3(:,:)=Id(:,:)
        CALL Angles_zxz(phi1,the1,psi1, Trans3 )
        CALL Spin_Rotation(Rtmp,phi1,the1,psi1)
        
        Rspin(:,:,iorb) = Rtmp
        
        if (Qprint) then
           WRITE(6,'(A,I3)') 'Spin-rotations for iorb=', iorb
           WRITE(6,*) 'angles=', phi1, the1, psi1
           WRITE(6,*) 'Rspin='
           do i=1,2
              write(6,'(2f12.6,2x,2f12.6)') (Rspin(i,j,iorb),j=1,2)
           enddo
        endif
     enddo
  ENDDO
END SUBROUTINE GetSpinRotation

REAL*8 FUNCTION acos_(x)
  IMPLICIT NONE
  REAL*8, intent(in) :: x
  if (abs(x).lt.1.d0) then
     acos_ = acos(x)
     return
  endif
  acos_ = 0.d0
  if (x.ge.1.d0) return
  if (x.le.-1.d0) acos_ = acos(-1.d0)
  return
END FUNCTION acos_


SUBROUTINE Angles_zxz(phi,the,psi,U)
  ! Routine gets rotation matrix U -- (only proper rotations are allowed)
  ! It returns three Euler angles (phi,the,psi), which generate rotation U by
  ! U = Rz(psi)*Rx(the)*Rz(phi)
  IMPLICIT NONE
  REAL*8, intent(out):: phi,the,psi
  REAL*8, intent(in) :: U(3,3)
  !
  REAL*8 :: acos_
  REAL*8 :: PI
  PARAMETER (PI = 3.141592653589793D0)

  the = acos_(U(3,3))
  
  if (abs(abs(U(3,3))-1)>1e-4) then
     phi = atan2(U(3,1),U(3,2))
     psi = PI-atan2(U(1,3),U(2,3))
  else
     psi=0.
     if (U(2,1)==0 .and. U(1,2)==0) then
        phi = acos_(U(1,1))-cos(the)*psi
     else
        phi = atan2( U(2,1)*cos(the), U(2,2)*cos(the) )-psi*cos(the)
     endif
  endif

  !WRITE(6,*) 'Routine Angles called for rotation='
  !do i=1,3
  !   write(6,'(2f15.9,2x,2f15.9,2x,2f15.9)') (U(i,j),j=1,3)
  !enddo
  !WRITE(6,*) 'angles are', phi,the,psi
  
END SUBROUTINE Angles_zxz

!!!!!! Rotations of real space vectors !!!!!!!!!!!!!!!
SUBROUTINE Rz_3(R, psi)
  IMPLICIT NONE
  REAL*8, intent(out) :: R(3,3)
  REAL*8, intent(in)  :: psi
  !
  R(1,1)=cos(psi); R(1,2)=-sin(psi); R(1,3)=0
  R(2,1)=sin(psi); R(2,2)=cos(psi); R(2,3)=0
  R(3,1)=0; R(3,2)=0; R(3,3)=1
END SUBROUTINE Rz_3

SUBROUTINE Ry_3(R,the)
  IMPLICIT NONE
  REAL*8, intent(out) :: R(3,3)
  REAL*8, intent(in)  :: the
  !
  R(1,1)=cos(the); R(1,2)=0.0; R(1,3)=sin(the)
  R(2,1)=0.0; R(2,2)=1.0; R(2,3)=0.0
  R(3,1)=-sin(the); R(3,2)=0.0; R(3,3)=cos(the)
END SUBROUTINE Ry_3

SUBROUTINE Rx_3(R, phi)
  IMPLICIT NONE
  REAL*8, intent(out) :: R(3,3)
  REAL*8, intent(in)  :: phi
  !
  R(1,1)=1.0; R(1,2)=0.0; R(1,3)=0.0
  R(2,1)=0.0; R(2,2)=cos(phi); R(2,3)=-sin(phi)
  R(3,1)=0.0; R(3,2)=sin(phi); R(3,3)=cos(phi)
END SUBROUTINE Rx_3

!!!!!! Rotations of spinors !!!!!!!!!!!!!!!
SUBROUTINE Rz_2(R, psi)
  IMPLICIT NONE
  COMPLEX*16, intent(out):: R(2,2)
  REAL*8, intent(in)     :: psi
  ! locals
  COMPLEX*16 :: i
  i = (0,1)
  R(1,1)=exp(psi/2.*i); R(1,2)=0.0;
  R(2,1)=0.0;  R(2,2)=exp(-psi/2.*i)
END SUBROUTINE Rz_2

SUBROUTINE Ry_2(R, the)
  IMPLICIT NONE
  COMPLEX*16, intent(out):: R(2,2)
  REAL*8, intent(in)     :: the
  ! 
  R(1,1)=cos(the/2.); R(1,2)=sin(the/2.)
  R(2,1)=-sin(the/2.); R(2,2)=cos(the/2.)
END SUBROUTINE Ry_2
  
SUBROUTINE Rx_2(R, phi)
  IMPLICIT NONE
  COMPLEX*16, intent(out):: R(2,2)
  REAL*8, intent(in)     :: phi
  ! locals
  COMPLEX*16 :: i
  i = (0,1)
  R(1,1)=cos(phi/2.);  R(1,2)=sin(phi/2.)*i
  R(2,1)=sin(phi/2.)*i; R(2,2)=cos(phi/2.)
END SUBROUTINE Rx_2

SUBROUTINE Space_Rotation( R, phi,the,psi )
  IMPLICIT NONE
  REAL*8, intent(in)  :: phi, the, psi
  REAL*8, intent(out) :: R(3,3)
  ! locals
  REAL*8 :: R1(3,3), R2(3,3), R3(3,3)
  CALL Rz_3(R1, psi)
  CALL Rx_3(R2, the)
  CALL Rz_3(R3, phi)
  !R = matmul(matmul(R1,R2),R3)
  R = matmul(R1,R2)
  R = matmul(R,R3)
END SUBROUTINE Space_Rotation

SUBROUTINE Spin_Rotation( R, phi,the,psi )
  IMPLICIT NONE
  REAL*8, intent(in)     :: phi, the, psi
  COMPLEX*16, intent(out):: R(2,2)
  ! locals
  COMPLEX*16 :: R1(2,2), R2(2,2), R3(2,2)
  CALL Rz_2(R1, psi)
  CALL Rx_2(R2, the)
  CALL Rz_2(R3, phi)
  R = matmul(matmul(R1,R2),R3)
END SUBROUTINE Spin_Rotation

REAL*8 FUNCTION detx(a)
  IMPLICIT NONE
  REAL*8, intent(in) :: a(3,3)
  !
  REAL*8 :: cc
  cc = a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(3,1)*a(2,2)*a(1,3)-a(1,1)*a(3,2)*a(2,3)-a(2,1)*a(1,2)*a(3,3)
  detx = cc
  return
END FUNCTION detx



SUBROUTINE Angles_zyz2(a1,a2,a3,U)
  ! Routine gets rotation matrix U -- (only proper rotations are allowed)
  ! It returns three Euler angles (a1,a2,a3), which generate rotation U by
  ! U = Rz(a3)*Ry(a2)*Rz(a1)
  IMPLICIT NONE
  REAL*8, intent(out):: a1,a2,a3
  REAL*8, intent(in) :: U(3,3)
  ! locals
  REAL*8 :: sy
  sy = sqrt(U(3,2)**2 + U(3,1)**2)
  if (sy .gt. 1e-6) then
    a1 = -atan2( U(3,2),  U(3,1));
    a2 = -atan2( sy,      U(3,3));
    a3 = -atan2( U(2,3), -U(1,3));
 else
    a1 = -atan2(-U(2,1),  U(2,2));
    a2 = -atan2( sy,      U(3,3));
    a3 = 0.0;
  endif
end SUBROUTINE Angles_zyz2

SUBROUTINE Angles_zyz(a1,a2,a3,U)
  ! Routine gets rotation matrix U -- (only proper rotations are allowed)
  ! It returns three Euler angles (a1,a2,a3), which generate rotation U by
  ! U = Rz(a3)*Ry(a2)*Rz(a1)
  IMPLICIT NONE
  REAL*8, intent(out):: a1,a2,a3
  REAL*8, intent(in) :: U(3,3)
  ! locals
  REAL*8 :: PI, sy
  PI = acos(-1.0d0)
  sy = sqrt(U(3,2)**2 +  U(3,1)**2)
  if (sy .gt. 1e-6) then
     a1 = PI-atan2( U(3,2), U(3,1))
     a2 = atan2(sy,      U(3,3))
     a3 = PI-atan2(U(2,3), -U(1,3))
  else 
     a1 = PI-atan2(-U(2,1), U(2,2))
     a2 = atan2(sy,      U(3,3))
     a3 = PI
  endif
End SUBROUTINE Angles_zyz

REAL*8 Function fact(a)
  IMPLICIT NONE
  ! factorial
  INTEGER, intent(in) :: a
  ! local
  INTEGER :: i
  REAL*8  :: r
  r = 1.0
  do i=a,2,-1
     r = r * i
  enddo
  fact = r
End Function fact

REAL*8 Function Wigner_reduced(J, m, n, beta)
  IMPLICIT NONE
  REAL*8, intent(in) :: J, m, n, beta
  !  ///////////////////////////////////////////////////////////////////////////////// 
  !  //  J                                     1/2
  !  // d   (beta) = [(J+n)!(J-n)!(J+m)!(J-m)!]      \times
  !  //  m,n
  !  //                           k                        2J+n-m-2k              m-n+2k
  !  //        ---            (-1)            [           ]         [            ]
  !  //        \   __________________________ |cos(beta/2)|         |-sin(beta/2)|
  !  //        /   (J-m-k)!(J+n-k)!(k+m-n)!k! [           ]         [            ]
  !  //        ---
  !  //         k
  !  //
  !  // Input		J     : Rank
  !  // 		      m,n     : Momentum index
  !  // 			n     : Momentum index
  !  // 			beta  : Euler angle (degrees)
  !  //
  !  // see : https://en.wikipedia.org/wiki/Wigner_D-matrix
  REAL*8 :: fact
  REAL*8  :: cosb, msinb, dj_, signf, cos_, sin_
  INTEGER :: k, d1, d2, d3
  if( (J<0) .or. (abs(m)>J) .or. (abs(n)>J) ) then  ! 1st insure resaonable {J,m,n}  
     print *, "Unable to Determine Reduced Wigner Element for j=", J, " m'=", m, " m=", n
     call exit(1)
  endif
  cosb  =  cos(beta/2.0d0) ! Cosine factors needed
  msinb = -sin(beta/2.0d0) ! Sine factors needed
  dj_ = 0.0d0              ! Value of Wigner element
  signf = 1.0
  do k=0,int(2*J)
     d1 = int(J-m)-k;
     d2 = int(J+n)-k;
     d3 = int(m-n) + k;
     if( d1>=0  .and.  d2>=0  .and.  d3>=0 ) then
        cos_ = cosb**(int( 2*J - 2*k + n - m))  ! Cosine to power
        sin_ = msinb**(int( 2*k + m - n))       ! Sine to power
        dj_ = dj_ + signf*cos_*sin_/( fact(d1)*fact(d2)*fact(d3)*fact(k) )
     endif
     signf = -signf                             ! Switch sign 1 <-> -1
  enddo
  Wigner_reduced = dj_ * sqrt(fact(int(J+m))*fact(int(J-m))*fact(int(J+n))*fact(int(J-n)))
End Function Wigner_reduced

COMPLEX*16 Function Wigner_zyz(J, m, n, alpha, beta, gamma)
  IMPLICIT NONE
  REAL*8, intent(in) :: J, m, n, alpha, beta, gamma
  COMPLEX*16, PARAMETER:: IMAG=   (0.0D0,1.0D0)
  REAL*8 :: Wigner_reduced
  Wigner_zyz = exp(-(m*alpha+n*gamma)*IMAG) * Wigner_reduced(J,m,n,beta)
End Function Wigner_zyz


SUBROUTINE Test_Angles(Nw)
  IMPLICIT NONE
  INTEGER, intent(in) :: Nw
  REAL*8  :: PI, alpha,beta,gamma,a1,a2,a3, diff
  REAL*8  :: R1(3,3), R2(3,3), R3(3,3), U(3,3), V(3,3)
  INTEGER :: i, k, l
  PI = acos(-1.0d0)
!!! Testing Angles
  do i=1,Nw
     do k=1,Nw
        do l=1,Nw
           alpha = (2*PI*(i-1.))/(Nw-1.)
           beta =  (2*PI*(k-1.))/(Nw-1.)
           gamma = (2*PI*(l-1.))/(Nw-1.)
           CALL Rz_3(R1, alpha)
           CALL Ry_3(R2, beta)
           CALL Rz_3(R3, gamma)
           U = matmul( matmul(R3,R2),R1 )

           CALL Angles_zyz(a1,a2,a3,U)
           !CALL Angles_zyz(a1,a2,a3,U)
           !print *, 'angles=', a1/PI, a2/PI, a3/PI
  
           CALL Rz_3(R1, a1)
           CALL Ry_3(R2, a2)
           CALL Rz_3(R3, a3)
  
           V = matmul( matmul(R3,R2),R1 )
           R1 = U-V
           diff = sum(abs(R1))
           if (diff > 1e-6) print *, 'diff=', sum(abs(R1))
           !do i=1,3
           !   do k=1,3
           !      WRITE(6,'(F12.7,1x)',advance='no') R1(i,k)
           !   enddo
           !   WRITE(6,*)
           !enddo
        enddo
     enddo
  enddo
  return 
END SUBROUTINE Test_Angles


SUBROUTINE GetWignerSpinMatrix_compatible(WM, Rotation)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: WM(2,2)
  REAL*8, intent(in)      :: Rotation(3,3)
  ! locals
  COMPLEX*16 :: Wigner_zyz
  REAL*8 :: a1, a2, a3, j
  INTEGER :: im, in, n, m
  COMPLEX*16 :: dz
  
  ! Wigner rotation uses R_z(a3) * R_y(a2) * R_z(a1) order, hence we need corresponding angles
  !CALL Angles_zyz(a1,a2,a3, Rotation)
  CALL Angles_zyz2(a1,a2,a3, Rotation)

  !print *, 'angles=', a1/PI, a2/PI, a3/PI
  !CALL Rz_3(R1, a1)
  !CALL Ry_3(R2, a2)
  !CALL Rz_3(R3, a3)
  !V = matmul( matmul(R3,R2),R1 )
  !print *, 'diff=', sum(abs(V-Rotation))

  j = 0.5d0
  WM=0
  im=0
  ! Our convention is [[M(up,up),M(up,dn)],[M(dn,up),M(dn,dn)]]
  do m=int(2*j),-int(2*j),-2
     im = im+1; in=0
     do n=int(2*j),-int(2*j),-2
        in = in+1
        ! Wigner rotation is passive e^{-i*a3*Jz} e^{-i*a2*Jy} e^{-i*a1*Jz}
        ! we prefer active rotation, i.e., D = e^{i*a3*Jz} e^{i*a2*Jy} e^{i*a1*Jz}
        dz = Wigner_zyz(j, dble(m/2.), dble(n/2.), -a1, -a2, -a3)
        ! When applying Wigner rotation, the convention is T^new_i  = \sum_j T_j * D_{j,i}
        ! but we would like to use it as matrix-vector multiplication, T^new = WM * T
        ! hence we transpose
        WM(in,im) = dz   
     enddo
  enddo
END SUBROUTINE GetWignerSpinMatrix_Compatible

SUBROUTINE GetWignerSpinMatrix(WM, j, Rotation)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: WM(2,2)
  REAL*8, intent(in)      :: Rotation(3,3)
  REAL*8, intent(in)      :: j ! we expect half-integer j, although it should work for all
  ! locals
  COMPLEX*16 :: Wigner_zyz
  REAL*8 :: a1, a2, a3
  INTEGER :: im, in, n, m
  COMPLEX*16 :: dz
  ! Wigner rotation uses R_z(a1) * R_y(a2) * R_z(a3) order, hence we need corresponding angles
  !CALL Angles_zyz(a3,a2,a1, Rotation)
  CALL Angles_zyz2(a3,a2,a1, Rotation)
  WM=0
  im=0
  ! Our convention is [[M(up,up),M(up,dn)],[M(dn,up),M(dn,dn)]]
  do m=int(2*j),-int(2*j),-2
     im = im+1; in=0
     do n=int(2*j),-int(2*j),-2
        in = in+1
        ! Wigner active rotation is <m|e^{-i*a1*Jz} e^{-i*a2*Jy} e^{-i*a3*Jz}|n>
        dz = Wigner_zyz(j, dble(m/2.), dble(n/2.), a1, a2, a3)
        ! We have : Y_{jm} = \sum_n Y_{jn} D_{n,m}
        ! hence you should use  Y(Rr) = Y(r) * WM   or
        ! Y_{jm} = \sum_n Y_{j,n} W_{n,m}
        WM(im,in) = dz   
     enddo
  enddo
END SUBROUTINE GetWignerSpinMatrix
SUBROUTINE GetWignerOrbitalMatrix(WM, l, Rotation)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: WM(2*l+1,2*l+1)
  REAL*8, intent(in)      :: Rotation(3,3)
  INTEGER, intent(in)     :: l
  ! locals
  COMPLEX*16 :: Wigner_zyz
  REAL*8 :: a1, a2, a3, j
  INTEGER :: im, in, n, m
  COMPLEX*16 :: dz
  ! Wigner rotation uses R_z(a1) * R_y(a2) * R_z(a3) order, hence we need corresponding angles
  !CALL Angles_zyz(a3,a2,a1, Rotation)
  CALL Angles_zyz2(a3,a2,a1, Rotation)
  j = l
  WM=0
  im=0
  ! Our convention is [[M(-2,-2),M(-2,-1).....],....,[M(2,1),M(2,2)]]
  do m=-int(2*j),int(2*j),2
     im = im+1; in=0
     do n=-int(2*j),int(2*j),2
        in = in+1
        ! Wigner active rotation is <m|e^{-i*a1*Jz} e^{-i*a2*Jy} e^{-i*a3*Jz}|n>
        dz = Wigner_zyz(j, dble(m/2.), dble(n/2.), a1, a2, a3)
        ! We have : Y_{jm} = \sum_n Y_{jn} D_{n,m}
        ! hence you should use  Y(Rr) = Y(r) * WM  or
        ! Y_{jm} = \sum_n Y_{j,n} WM_{n,m}
        WM(im,in) = dz
     enddo
  enddo
END SUBROUTINE GetWignerOrbitalMatrix

!program sprot
!  IMPLICIT NONE
!  REAL*8 PI
!  PARAMETER (PI = 3.141592653589793D0)
!  REAL*8 :: j, beta, alpha, gamma, dj
!  COMPLEX*16 :: dz
!  INTEGER :: n, m, i, k, l, Nw, im, in
!  REAL*8 :: fact, Wigner_reduced
!  COMPLEX*16 :: Wigner_zyz
!  REAL*8 :: R1(3,3), R2(3,3), R3(3,3), U(3,3), V(3,3)
!  REAL*8 :: a1,a2,a3, diff
!  COMPLEX*16 :: RS2(2,2)
!  COMPLEX*16, allocatable :: WM(:,:)
!  
!  alpha = 0.5d0*PI
!  beta  = 0.5d0*PI
!  gamma = 0.0d0*PI
!
!  CALL Rz_3(R1, alpha)
!  CALL Ry_3(R2, beta)
!  CALL Rz_3(R3, gamma)
!  
!  U = matmul( matmul(R1,R2),R3 )
!  WRITE(6,*) 'The real space transformation'
!  do i=1,3
!     do k=1,3
!        WRITE(6,'(F12.7,1x)',advance='no') U(i,k)
!     enddo
!     WRITE(6,*)
!  enddo
!  
!  allocate( WM(2,2) )
!  Call GetWignerSpinMatrix(WM, 1.d0/2.d0, U)
!  
!  j=1.d0/2.d0
!  print *, "Wigner spin matrix="
!  do im=1,int(2*j)+1
!     do in=1,int(2*j)+1
!        WRITE(6,'(2F12.7,1x)',advance='no') WM(im,in)
!     enddo
!     WRITE(6,*)
!  enddo
!  Call GetWignerSpinMatrix_compatible(WM, U)
!  print *, "Old compatible Wigner spin matrix="
!  do im=1,int(2*j)+1
!     do in=1,int(2*j)+1
!        WRITE(6,'(2F12.7,1x)',advance='no') WM(im,in)
!     enddo
!     WRITE(6,*)
!  enddo
!  
!  deallocate(WM)
!
!  j=1.d0
!  allocate( WM(3,3) )
!  Call GetWignerOrbitalMatrix(WM, 1, U)
!  
!  print *, "Wigner l=1 matrix="
!  do im=1,int(2*j)+1
!     do in=1,int(2*j)+1
!        WRITE(6,'(2F12.7,1x)',advance='no') WM(im,in)
!     enddo
!     WRITE(6,*)
!  enddo
!  deallocate( WM )
!  
!  !CALL Angles_zxz(a1,a2,a3,U)
!  !
!  !print *, 'angles=', a1/PI, a2/PI, a3/PI
!  !CALL Rz_3(R1, a1)
!  !CALL Rx_3(R2, a2)
!  !CALL Rz_3(R3, a3)
!  !V = matmul( matmul(R3,R2),R1 )
!  !print *, 'diff=', sum(abs(U-V))
!  !
!  !CALL Spin_Rotation(RS2,a1,a2,a3)
!  !print *, 'old spin rotation'
!  !do i=1,2
!  !   do k=1,2 
!  !      WRITE(6,'(2F12.7,1x)',advance='no') RS2(i,k)
!  !   end do
!  !   print *
!  !end do
!  
!end program sprot
