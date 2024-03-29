SUBROUTINE RECFIL(FNAME1,GMAX,RCFILE,KXMAX,KYMAX,KZMAX,NWAVE)
  !
  USE param, ONLY: nsym
  use defs,  ONLY: PI
  use structure, ONLY: ALPHA, ortho, lattic, AA, BB, CC, PIA, BR2
  USE sym2,      ONLY: iord
  use reclat
  USE com_mpi,ONLY: myrank, master, Qprint
  !IMPLICIT REAL*8 (A-H,O-Z)
  IMPLICIT NONE
  CHARACTER*180, intent(in) :: FNAME1
  CHARACTER*4,   intent(in) :: RCFILE
  REAL*8,        intent(in) :: GMAX
  INTEGER,       intent(out):: KXMAX, KYMAX, KZMAX, NWAVE
  ! locals
  LOGICAL    :: NEWRC
  INTEGER    :: INDMAX, ITAPE, itape2, iord1, i
  REAL*8     :: ALPHA1(3), AA1, BB1, CC1, GMAX1, t1, t2

  newrc = .FALSE.

  AA1=AA
  BB1=BB
  CC1=CC
  ALPHA1(1:3)=ALPHA(1:3)

  OPEN(13,FILE=FNAME1,STATUS='UNKNOWN',FORM='UNFORMATTED')
  if (RCFILE.eq.'FILE') then
     read(13,iostat=itape) NWAVE,INDMAX,GMAX1,iord1
     if(itape.eq.0) THEN
        if (Qprint) write(6,*) 'reading recprlist from file'
        read(13,iostat=itape2) AA1,BB1,CC1
        if(itape2.ne.0) then
           NEWRC=.TRUE.
           goto 101
        endif
        read(13,iostat=itape2) ALPHA1
        if(itape2.ne.0) then
           NEWRC=.TRUE.
           goto 101
        endif
        read(13,iostat=itape2) KXMAX,KYMAX,KZMAX
        if(itape2.ne.0) then
           NEWRC=.TRUE.
           goto 101
        endif
        allocate (kzz(3,nwave),inst(nwave),tauk(nwave*NSYM))
        read(13,iostat=itape2) (INST(i),I=1,NWAVE)
        if(itape2.ne.0) then
           deallocate (kzz,inst,tauk)
           NEWRC=.TRUE.
           goto 101
        endif
        read(13,iostat=itape2) (KZZ(1,I),KZZ(2,I),KZZ(3,I),I=1,NWAVE)
        if(itape2.ne.0) then
           deallocate (kzz,inst,tauk)
           NEWRC=.TRUE.
           goto 101
        endif
        read(13,iostat=itape2) (TAUK(I),I=1,INDMAX)
        if(itape2.ne.0) then
           deallocate (kzz,inst,tauk)
           NEWRC=.TRUE.
           goto 101
        endif
        
        !           check if we are allowed to use stored recprlist
        if(abs(GMAX1-GMAX).gt.1.d-8)         NEWRC=.TRUE.
        if(abs(AA1-AA).gt.1.d-8)             NEWRC=.TRUE.
        if(abs(BB1-BB).gt.1.d-8)             NEWRC=.TRUE.
        if(abs(CC1-CC).gt.1.d-8)             NEWRC=.TRUE.
        if(abs(ALPHA1(1)-ALPHA(1)).gt.1.d-8) NEWRC=.TRUE.
        if(abs(ALPHA1(2)-ALPHA(2)).gt.1.d-8) NEWRC=.TRUE.
        if(abs(ALPHA1(3)-ALPHA(3)).gt.1.d-8) NEWRC=.TRUE.
        if(iord.NE.iord1)                    NEWRC=.TRUE.
        
        if (NEWRC .and. Qprint) write (6,*) GMAX1,GMAX,AA1,AA,BB1,BB,CC1,CC, ALPHA1(1),ALPHA(1),ALPHA1(2),ALPHA(2),ALPHA1(3), ALPHA(3),iord1,iord
        if(newrc) deallocate (kzz,inst,tauk)                           
     else
        NEWRC=.TRUE.
     endif
  else
     !        if we don't use the data from file we always calculate recpr
     NEWRC=.TRUE.
  endif
101 continue
  close(13)

  if(NEWRC) then
     !        generate a new recprlist 
     if (Qprint) write(6,*) 'generate new recprlist'
     if (ORTHO.or.lattic(1:1).eq.'R') then
        call deter(GMAX,PIA,BR2,KXMAX,KYMAX,KZMAX,LATTIC)
     else
        KXMAX=GMAX*AA/2.d0/PI+1
        KYMAX=GMAX*BB/2.d0/PI+1
        KZMAX=GMAX*CC/2.d0/PI+1
     endif
     if (lattic(1:1).eq.'R') then
        KXMAX=max(kxmax,int(GMAX*AA/2.d0/PI+1))
        KYMAX=max(kymax,int(GMAX*BB/2.d0/PI+1))
        KZMAX=max(kzmax,int(GMAX*CC/6.d0/PI+1))
     endif
     call cputim(t1)
     CALL RECPR (NWAVE,INDMAX,KXMAX,KYMAX,KZMAX,GMAX)
     call cputim(t2)
     if (Qprint) write(6,*) ' time in recpr: ',t2-t1
     if(RCFILE.EQ.'FILE') then
        !           write new recprlist into file
        OPEN(13,FILE=FNAME1,STATUS='UNKNOWN',FORM='UNFORMATTED')
        write(13) NWAVE,INDMAX,GMAX,iord
        write(13) AA,BB,CC
        write(13) ALPHA
        write(13) KXMAX,KYMAX,KZMAX
        write(13) (INST(i),I=1,NWAVE)
        write(13) (KZZ(1,I),KZZ(2,I),KZZ(3,I),I=1,NWAVE)
        write(13) (TAUK(I),I=1,INDMAX)
        close(13)
     endif
  endif
END SUBROUTINE RECFIL
