! @Copyright 2007 Kristjan Haule
! 
SUBROUTINE ReadThirdPart_of_in1(ITAPE,Emin,Emax,nband,fixed_nbands, &
     sproj_limit,eshift_iter,fh_in1)
  IMPLICIT NONE
  REAL*8,  intent(out) :: Emin, Emax, sproj_limit, eshift_iter
  INTEGER, intent(out) :: ITAPE, nband
  INTEGER, intent(in)  :: fh_in1
  LOGICAL, intent(out) :: fixed_nbands
  ! locals
  CHARACTER*80 :: line, keyword_line
  INTEGER      :: ios, i, code
  REAL*8       :: test1, test2
  
  read(fh_in1,'(a)')line
  ! Preserve the legacy numerical input exactly.  Fixed-band selection is
  ! enabled only by an explicit, case-insensitive trailing keyword, e.g.
  !   K-VECTORS FROM UNIT:4 -10.0 2.0 66 FIXED_NBANDS
  ! List-directed input ignores the trailing word after reading the four
  ! numerical fields, so the line remains readable by an unpatched WIEN2k.
  keyword_line=line
  do i=1,len(keyword_line)
     code=iachar(keyword_line(i:i))
     if (code.ge.iachar('a') .and. code.le.iachar('z')) &
          keyword_line(i:i)=achar(code-iachar('a')+iachar('A'))
  enddo
  fixed_nbands=index(keyword_line,'FIXED_NBANDS').gt.0
  nband=99999999
  read(line(21:),*,IOSTAT=ios) ITAPE,Emin,Emax,nband  ! ITAPE==klist unit ; Emin==min_E ; Emax==max_E ; nband==maximum number of bands
  if (ios.ne.0) then
     read(line(21:),*,IOSTAT=ios) ITAPE,Emin,Emax   ! maybe nbands is missing ; can be in two formats, either [int,real,real] or [int,real,real,it]
  endif
  read(fh_in1,'(a)',IOSTAT=ios)line              ! maybe case.in1 contains also sproj_limit ?
  if (ios.eq.0) then
     read(line,*,IOSTAT=ios) test1,test2   ! parse possible sproj_limit and eshift_iter
     if (ios.ne.0) then                    ! did not find both numbers. maybe just one?
        read(line,*,IOSTAT=ios) test1
        if (ios.ne.0) then
           sproj_limit=1.0d-15
           eshift_iter=0.d0
        else
           if(test1.gt.0.d0 .and. test1.lt.1.d-6) then
              sproj_limit=test1
              eshift_iter=0.d0
           else if(abs(test1).lt.4.d0 .and. abs(test1).gt.1.d-6) then
              sproj_limit=1.0d-15
              eshift_iter=test1
           else
              sproj_limit=1.0d-15
              eshift_iter=0.d0 
           endif
        endif
     else
        if(test1.gt.0.d0 .and. test1.lt.1.d-6) then
           sproj_limit=test1
           eshift_iter=0.d0
           if(abs(test2).lt.4.d0 .and. abs(test2).gt.1.d-6) eshift_iter=test2
        else if(abs(test1).lt.4.d0 .and. abs(test1).gt.1.d-6) then
           sproj_limit=1.0d-15
           eshift_iter=test1
           if(test2.gt.0.d0 .and. test2.lt.1.d-6) sproj_limit=test2
        else
           sproj_limit=1.0d-15
           eshift_iter=0.d0 
        endif
     endif
  else
     sproj_limit=1.0d-15   ! default value 
     eshift_iter=0.d0      ! default value
  endif
  !
END SUBROUTINE ReadThirdPart_of_in1
