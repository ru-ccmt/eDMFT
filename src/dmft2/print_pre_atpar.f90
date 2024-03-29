SUBROUTINE Print_pre_atpar(jatom,lfirst, el, lmax2)
  use structure,  ONLY: iatnr, mult, pos, rotij, tauij, rotloc, aname
  IMPLICIT NONE
  INTEGER, intent(in) :: jatom, lfirst, lmax2
  REAL*8, intent(in)  :: el(0:lmax2)
  ! locals
  INTEGER :: i, jrow, jcol, index, m, jc, jr
  !
  WRITE(6,13)  JATOM,IATNR(JATOM),(POS(I,lfirst),I=1,3),MULT(JATOM)        
  WRITE(6,1060) ANAME(JATOM)                                    
  DO JROW=1,3                                               
     WRITE(6,1070) ( ROTLOC(JCOL,JROW,JATOM),JCOL=1,3 )         
  ENDDO
  INDEX=lfirst-1                                                 
  DO M=1,MULT(JATOM)                                        
     INDEX=INDEX+1                                               
     WRITE(6,1080) M,( POS(JC,INDEX),JC=1,3 )                    
     DO JR=1,3                                               
        WRITE(6,1090)(ROTIJ(JC,JR,INDEX),JC=1,3),TAUIJ(JR,INDEX)   
     ENDDO
  ENDDO
  WRITE(6,7) ANAME(JATOM)                                           
  WRITE(6,5) el                                                      
  WRITE(6,14)
13 FORMAT(////,':POS',i3.3,':',1x,'AT.NR.',I4,1X,'POSITION =',3F8.5,2X,'MULTIPLICITY =',I3)
1060 FORMAT(//,3X,'NOT EQUIV ATOM ',A10,'  LOCAL ROTATION MATRIX')     
1070 FORMAT(30X,3F10.5)                                                
1080 FORMAT(13X,'EQUIV ATOM ',I3,3X,'POSITION: ',3F8.3)                
1090 FORMAT(30X,3F10.5,5X,F10.5)                                       
7 FORMAT(/10X,'ATOMIC PARAMETERS FOR ',A10/)                        
5 FORMAT(10X,' ENERGY PARAMETERS ARE',7F7.2)
14 FORMAT(/11X,'L',5X,'U(R)',10X,'U''(R)',9X,'DU/DE',8X,'DU''/DE',6X,'NORM-U''')
END SUBROUTINE Print_pre_atpar
