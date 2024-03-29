SUBROUTINE PrintForces2(forcea,fsph,fsph2,fnsp,fvdrho,forb,fextra)
  !     prints partial forces
  use structure, only : mult, natm
  USE com,   only : nat
  IMPLICIT NONE
  !IMPLICIT REAL*8 (A-H,O-Z)
  logical, intent(in) :: forcea(0:3,nat)
  real*8, intent(in)  :: fsph(3,natm),fsph2(3,natm),fnsp(3,natm),fvdrho(3,nat), forb(3,natm), fextra(3,nat)
  ! local variables
  real*8 :: fpsi(3,nat), fequ(3), ftot(3,nat)
  real*8 :: fequ_nrm, fnsp_nrm, fsph_nrm, fsph2_nrm, forb_nrm, fvdrho_nrm, fpsi_nrm, ftot_nrm, fextra_nrm
  INTEGER :: jatom, ia, ia1, ia2, ik
  
  fpsi(:,:)=0.d0
  DO jatom=1,nat
     ia2 = sum(mult(1:jatom))
     ia1 = ia2-mult(jatom)+1
     DO ia=ia1,ia2
        fequ(:)=0.d0
        DO ik=1,3
           if (forcea(ik,jatom)) fequ(ik)=(fsph(ik,ia)+fsph2(ik,ia)+forb(ik,ia)+fnsp(ik,ia))/mult(jatom)
        ENDDO
        fpsi(:,jatom) = fpsi(:,jatom) + fequ(:) ! averaged force in local coordiante system
        fequ_nrm = sqrt(dot_product(fequ,fequ))
        fnsp_nrm = sqrt(dot_product(fnsp(:,ia),fnsp(:,ia)))     ! Force2
        fsph_nrm = sqrt(dot_product(fsph(:,ia),fsph(:,ia)))     ! Force1
        fsph2_nrm= sqrt(dot_product(fsph2(:,ia),fsph2(:,ia)))   ! Force3
        forb_nrm = sqrt(dot_product(forb(:,ia),forb(:,ia)))     ! ForceU
        write(6,77) jatom,ia-ia1+1,'+ SPH', fsph_nrm*1000,  (fsph(ik,ia)*1000,ik=1,3)
        write(6,77) jatom,ia-ia1+1,'+ SP2', fsph2_nrm*1000, (fsph2(ik,ia)*1000,ik=1,3)
        write(6,77) jatom,ia-ia1+1,'+ ORB', forb_nrm*1000,  (forb(ik,ia)*1000,ik=1,3)
        write(6,77) jatom,ia-ia1+1,'+ NSP', fnsp_nrm*1000,  (fnsp(ik,ia)*1000,ik=1,3)
        write(6,77) jatom,ia-ia1+1,'> EQU', fequ_nrm*1000,  (fequ(ik)*1000,ik=1,3)
     ENDDO
     !     calculate partial forces 
     ftot(:,jatom) = fvdrho(:,jatom) + fpsi(:,jatom) + fextra(:,jatom)
     fvdrho_nrm    = sqrt(dot_product(fvdrho(:,jatom),fvdrho(:,jatom)))  ! Force4
     fpsi_nrm      = sqrt(dot_product(fpsi(:,jatom),fpsi(:,jatom)))
     fextra_nrm    = sqrt(dot_product(fextra(:,jatom),fextra(:,jatom)))  ! Extra force due to discountinuity of the basis functions
     ftot_nrm      = sqrt(dot_product(ftot(:,jatom),ftot(:,jatom)))
     write(6,77) jatom,1,'= PSI', fpsi_nrm*1000,   (fpsi(ik,jatom)*1000,ik=1,3)
     write(6,77) jatom,1,'+ VDR', fvdrho_nrm*1000, (fvdrho(ik,jatom)*1000,ik=1,3)
     write(6,77) jatom,1,'+ EXT', fextra_nrm*1000, (fextra(ik,jatom)*1000,ik=1,3)
     write(6,77) jatom,1,'> TOT', ftot_nrm*1000,   (ftot(ik,jatom)*1000,ik=1,3)
     write(21,78) 
     write(21,79) jatom,jatom, ftot_nrm*1000, (ftot(ik,jatom)*1000,ik=1,3)
     write(6,*)
     
77   format(2i3,a7,4f17.9)
78   FORMAT (7x,'VALENCE-FORCE IN mRy/a.u. = |F|',3x,'Fx',13x,'Fy',13x,'Fz')
79   FORMAT (':FVA',i3.3,':',1x,i3,'.ATOM',4f17.9)
  ENDDO
  RETURN
END SUBROUTINE PrintForces2

SUBROUTINE PrintForces(forcea,fsph,fsph2,fnsp,fvdrho,forb)
  !     prints partial forces
  use structure, only : mult, natm
  USE com,   only : nat
  IMPLICIT NONE
  !IMPLICIT REAL*8 (A-H,O-Z)
  logical, intent(in) :: forcea(0:3,nat)
  real*8, intent(in)  :: fsph(3,natm),fsph2(3,natm),fnsp(3,natm),fvdrho(3,nat), forb(3,natm)
  ! local variables
  real*8 :: fpsi(3,nat), fequ(3), ftot(3,nat)
  real*8 :: fequ_nrm, fnsp_nrm, fsph_nrm, fsph2_nrm, forb_nrm, fvdrho_nrm, fpsi_nrm, ftot_nrm
  INTEGER :: jatom, ia, ia1, ia2, ik
  
  fpsi(:,:)=0.d0
  DO jatom=1,nat
     ia2 = sum(mult(1:jatom))
     ia1 = ia2-mult(jatom)+1
     DO ia=ia1,ia2
        fequ(:)=0.d0
        DO ik=1,3
           if (forcea(ik,jatom)) fequ(ik)=(fsph(ik,ia)+fsph2(ik,ia)+forb(ik,ia)+fnsp(ik,ia))/mult(jatom)
        ENDDO
        fpsi(:,jatom) = fpsi(:,jatom) + fequ(:) ! averaged force in local coordiante system
        fequ_nrm = sqrt(dot_product(fequ,fequ))
        fnsp_nrm = sqrt(dot_product(fnsp(:,ia),fnsp(:,ia)))     ! Force2
        fsph_nrm = sqrt(dot_product(fsph(:,ia),fsph(:,ia)))     ! Force1
        fsph2_nrm= sqrt(dot_product(fsph2(:,ia),fsph2(:,ia)))   ! Force3
        forb_nrm = sqrt(dot_product(forb(:,ia),forb(:,ia)))     ! ForceU
        write(6,77) jatom,ia-ia1+1,'+ SPH', fsph_nrm*1000,  (fsph(ik,ia)*1000,ik=1,3)
        write(6,77) jatom,ia-ia1+1,'+ SP2', fsph2_nrm*1000, (fsph2(ik,ia)*1000,ik=1,3)
        write(6,77) jatom,ia-ia1+1,'+ ORB', forb_nrm*1000,  (forb(ik,ia)*1000,ik=1,3)
        write(6,77) jatom,ia-ia1+1,'+ NSP', fnsp_nrm*1000,  (fnsp(ik,ia)*1000,ik=1,3)
        write(6,77) jatom,ia-ia1+1,'> EQU', fequ_nrm*1000,  (fequ(ik)*1000,ik=1,3)
     ENDDO
     !     calculate partial forces 
     ftot(:,jatom) = fvdrho(:,jatom) + fpsi(:,jatom)
     fvdrho_nrm    = sqrt(dot_product(fvdrho(:,jatom),fvdrho(:,jatom)))  ! Force4
     fpsi_nrm      = sqrt(dot_product(fpsi(:,jatom),fpsi(:,jatom)))
     ftot_nrm      = sqrt(dot_product(ftot(:,jatom),ftot(:,jatom)))
     write(6,77) jatom,1,'= PSI', fpsi_nrm*1000,   (fpsi(ik,jatom)*1000,ik=1,3)
     write(6,77) jatom,1,'+ VDR', fvdrho_nrm*1000, (fvdrho(ik,jatom)*1000,ik=1,3)
     write(6,77) jatom,1,'> TOT', ftot_nrm*1000,   (ftot(ik,jatom)*1000,ik=1,3)
     write(21,78) 
     write(21,79) jatom,jatom, ftot_nrm*1000, (ftot(ik,jatom)*1000,ik=1,3)
     write(6,*)
     
77   format(2i3,a7,4f15.7)
78   FORMAT (7x,'VALENCE-FORCE IN mRy/a.u. = |F|',3x,'Fx',13x,'Fy',13x,'Fz')
79   FORMAT (':FVA',i3.3,':',1x,i3,'.ATOM',4f17.9)
  ENDDO
  RETURN
END SUBROUTINE PrintForces
