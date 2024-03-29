subroutine read_orbital_potential()
  USE orb, ONLY: iorbpot, nmod, nsp, natorb, nonso, Bext
  USE mpi, ONLY: Qprint
  IMPLICIT NONE
  ! locals
  INTEGER      :: isi, i, j, index, natorb0, jat, iatom, nlorb, nl, l, m, m1
  REAL*8       :: rval, cval
  character*25 :: charop(3)
  character*5  :: charsp (-1:2)
  complex*16   :: vorb(-3:3,-3:3)

  !  orbital potential part: start
  natorb=0
  natorb0=0
  nonso=0
  nmod=0
  if (iorbpot.eq.1) then
     ! read the orbital potential
     charop(1)=' LDA+U potential'
     charop(2)=' Orbital polarization'
     charop(3)=' Interaction with Bext'
     charsp(-1)=' down'
     charsp(0)= ' para'
     charsp(1)= ' up  '
     charsp(2)= ' dnup'
     do isi=1,3
        index=0
        read(10+isi,*,end=555)nmod,nsp,natorb,Bext  ! 11 and 12 are orbital potential files
        DO jat=1,natorb
           read(10+isi,*)iatom,nlorb
           do nl=1,nlorb
              read(10+isi,*) l
              if (Qprint) then
                 write(6,556)charop(nmod),iatom,l,charsp(nsp)
                 write(8,556)charop(nmod),iatom,l,charsp(nsp)
              endif
              do i=-l,l
                 do j=-l,l
                    read(10+isi,*) rval,cval
!!!!! BUG BUG BUG
                    !vorb(i,j) = cmplx(rval,cval)!!!! HERE IS W2K bug
                    vorb(i,j) = dcmplx(rval,cval)!!!! correct
                 enddo
              enddo
              call vorblo(index,iatom,l,isi,vorb)
              if(isi.eq.1)NATORB0=index
              if (Qprint) then
                 write(6,*)
                 write(6,*)' Orbital potential real part'
                 do m=-l,l
                    write(6,102)m,(dble(vorb(m,m1)),m1=-l,l)
                 enddo
                 write(6,*)
                 write(6,*)' Orbital potential imaginary part'
                 do m=-l,l
                    write(6,102)m,(dimag(vorb(m,m1)),m1=-l,l)
                 enddo
                 write(6,*)
              end if
           enddo
           ! end of orbital potential input
        endDO
        if(nonso.gt.0 .and. Qprint)then
           write(6,*)
           write(6,*) ' **** spin-orbit is excluded from the calculation ****'
           write(6,*)
           write(8,*)
           write(8,*) ' **** spin-orbit is excluded from the calculation ****'
           write(8,*)
        endif
555     continue
     enddo
     !  orbital potential part: end
  endif
  natorb=natorb0
  
102 format(' M=',i3,7f11.5)
556 format(a22,' added for atom type',i3,' L=',i3,' spin block',a5)
end subroutine read_orbital_potential
