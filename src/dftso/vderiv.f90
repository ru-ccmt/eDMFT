subroutine vderiv(lpot,jatom,lspin,rx,vru,vrr,vder,jri,jspin)
  USE param, ONLY: nrad, nato
  implicit none
  integer, intent(in) :: jspin
  real*8, intent(in)  :: vru(nrad,nato,jspin)       ! potential as written in the file
  real*8, intent(out) :: vrr(nrad,2), vder(nrad,2)  ! real potential and its derivative
  real*8, intent(in)  :: rx(nrad)                   ! real space mesh
  integer, intent(in) :: lpot, lspin, jri, jatom
  ! locals
  real*8 :: der
  real*8 :: x(5),y(5) 
  integer isi,i
  !
  do isi=1,lspin   ! 30
     !.....use 5-order-langrange-interpretaion to calculate dV/dr.
     ! Notice that Vru from file is stored in a  special form Vru==V_file = V_{Rydberg-units}(r)*sqrt(4*pi)/r
     do i=1,jri
        if (lpot.eq.3) then
           vrr(i,1) = ( (vru(i,jatom,1) + vru(i,jatom,jspin))/2. ) / rx(i)
           vrr(i,2) = vrr(i,1)
        else 
           vrr(i,isi) = vru(i,jatom,isi) / rx(i)
        endif
     enddo
     !.....VRR is the real potential in units of Ryd.
     DO I=1,jri    ! 30
        IF(I.EQ.1) THEN
           X(1)=RX(1)
           X(2)=RX(2)
           X(3)=RX(3)
           X(4)=RX(4)
           X(5)=RX(5)
           Y(1)=VRR(1,isi)
           Y(2)=VRR(2,isi)
           Y(3)=VRR(3,isi)
           Y(4)=VRR(4,isi)
           Y(5)=VRR(5,isi)
        ELSE IF(I.EQ.2) THEN
           X(1)=RX(2)
           X(2)=RX(1)
           X(3)=RX(3)
           X(4)=RX(4)
           X(5)=RX(5)
           Y(1)=VRR(2,isi)
           Y(2)=VRR(1,isi)
           Y(3)=VRR(3,isi)
           Y(4)=VRR(4,isi)
           Y(5)=VRR(5,isi)
        ELSE IF(I.EQ.jri-1) THEN
           X(1)=RX(jri-1)
           X(2)=RX(jri)  
           X(3)=RX(jri-2)
           X(4)=RX(jri-3)
           X(5)=RX(jri-4)
           Y(1)=VRR(jri-1,isi)
           Y(2)=VRR(jri,isi)  
           Y(3)=VRR(jri-2,isi)
           Y(4)=VRR(jri-3,isi)
           Y(5)=VRR(jri-4,isi)
        ELSE IF(I.EQ.jri) THEN
           X(1)=RX(jri)
           X(2)=RX(jri-1)
           X(3)=RX(jri-2)
           X(4)=RX(jri-3)
           X(5)=RX(jri-4)
           Y(1)=VRR(jri,isi)
           Y(2)=VRR(jri-1,isi)
           Y(3)=VRR(jri-2,isi)
           Y(4)=VRR(jri-3,isi)
           Y(5)=VRR(jri-4,isi)
        ELSE
           X(1)=RX(I)
           X(2)=RX(I-1)
           X(3)=RX(I-2)
           X(4)=RX(I+1)
           X(5)=RX(I+2)
           Y(1)=VRR(I,isi)
           Y(2)=VRR(I-1,isi)
           Y(3)=VRR(I-2,isi)
           Y(4)=VRR(I+1,isi)
           Y(5)=VRR(I+2,isi)
        END IF
        CALL LAGDER(5,X,Y,DER)
        !  vder is dV/dr in units of Ryd/a.u.
        if(lpot.eq.3)then
           vder(i,1)=der  
           vder(i,2)=der  
        else
           VDER(i,isi)=DER
        endif
     END DO
  end do
  return 
end subroutine vderiv
