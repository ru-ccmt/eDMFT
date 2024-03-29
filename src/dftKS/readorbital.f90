! @Copyright 2007 Kristjan Haule
!
SUBROUTINE ReadOrbitalPotential(nat,runorb)
  use orb, only     : BEXT, VORB, NMOD, NSP, NATORB, IATOM, NLORB, LORB, init_orb
  IMPLICIT NONE
  INTEGER, intent(in) :: nat
  LOGICAL, intent(in) :: runorb
  !locals
  character*25 :: charop(3)
  character*5  :: charsp (-1:1)
  INTEGER          :: jat, nl, l, m, m1, i, j
  DOUBLE PRECISION :: rval,cval
  
  call init_orb(NAT)
  ! read the orbital potential
  if(runorb) then
     charop(1)=' LDA+U potential'
     charop(2)=' Orbital polarization'
     charop(3)=' Interaction with Bext'
     charsp(-1)=' down'
     charsp(0)= ' para'
     charsp(1)= ' up  '
     nmod=0
     read(7,*,end=555,err=555)nmod,nsp,natorb,Bext
     DO JAT = 1, NATORB
        read(7,*)iatom(jat),nlorb(jat)
        do nl=1,nlorb(jat)
           read(7,*)lorb(nl,jat)
           l=lorb(nl,jat)
           write(6,556)charop(nmod),iatom(jat),lorb(nl,jat),charsp(nsp)
           write(21,556)charop(nmod),iatom(jat),lorb(nl,jat),charsp(nsp)
556        format(a22,' added for atom type',i3,' L=',i3,' spin',a5)
           do i=-l,l
              do j=-l,l
                 read(7,*) rval,cval
                 vorb(jat,l,i,j)=dcmplx(rval,cval)
              END DO
           enddo
102        format(' M=',i3,7f11.5)
           write(6,*)
           write(6,*)' Orbital potential real part'
           do m=-l,l
              write(6,102)m,(dble(vorb(jat,l,m,m1)),m1=-l,l)
           enddo
           write(6,*)
           write(6,*)' Orbital potential imaginary part'
           do m=-l,l
              write(6,102)m,(dimag(vorb(jat,l,m,m1)),m1=-l,l)
           enddo
           write(6,*)
        END DO
     END DO
  end if
555 continue
END SUBROUTINE ReadOrbitalPotential
