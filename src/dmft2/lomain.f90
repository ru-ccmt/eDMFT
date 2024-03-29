
module mod_lo ! We need moule to be able to pass allocatable arrays
CONTAINS
  SUBROUTINE lomain(XROTLOC,is,nemin,nemax,lfirst,latom,n,jatom,alm,blm,clm,Qforce_j,aalm,bblm,cclm,lmax2)
    USE param,ONLY: lomax, nume, nloat
    USE defs, ONLY: TWO, PI, IMAG
    USE structure,ONLY: mult, rotij, tauij, BR1, POS, ROTLOC
    USE lo, ONLY: nlov, ilo, Alo, Blo, Clo !, nlo, nlov, nlon, 
    USE xa, ONLY: BKROT, BKRLOC!, BK
    USE xa3,ONLY: BK3lo, k3lo, As_lo
    USE com,ONLY: 
    USE dmfts, ONLY: shft
    IMPLICIT NONE
    REAL*8,  INTENT(in) :: XROTLOC(3,3)
    INTEGER, INTENT(in) :: is, nemin,nemax,lfirst,latom,n,jatom, lmax2
    COMPLEX*16, INTENT(inout) :: alm((LMAX2+1)*(LMAX2+1),nume),blm((LMAX2+1)*(LMAX2+1),nume),clm((lomax+1)*(lomax+1),nume,nloat)
    LOGICAL, intent(in) :: Qforce_j
    COMPLEX*16, allocatable, intent(inout) :: aalm(:,:,:), bblm(:,:,:), cclm(:,:,:,:)
    !
    COMPLEX*16, allocatable :: YL(:),  phs(:)
    COMPLEX*16 :: PHSHEL                     
    INTEGER    :: l,jlo,jneq,i,m,m1,num,index, nmult, nbands_dft, lmax2lmax2, l2l2
    real*8     :: Krot(3), twopi,arg1,arg2,argt2, BKROT2(3), rotloc_x_BR1(3,3), rotloc_x_BR1_x_rotij(3,3), xrotloc_x_BR1(3,3)
    complex*16 :: yp,ayp,byp,cyp
    INTEGER, allocatable :: iind(:,:,:)

    allocate(YL((lomax+1)*(lomax+1)), phs(nume) )

    if (Qforce_j) then
       lmax2lmax2 = (lmax2+1)*(lmax2+1)
       if (size(aalm,3).NE.3 .or. size(cclm,4).NE.3) print*, 'ERROR lomain: Third dimension of aalm should be 3 but is', size(aalm,3)
       if (size(aalm,1).NE.lmax2lmax2 .or. size(cclm,1).NE.(lomax+1)*(lomax+1)) print*, 'ERROR lomain: first dimension of aalm should be', lmax2lmax2, 'but is', size(aalm,1)
       if (size(cclm,3).NE.nloat) print*, 'ERROR lomain: Second dimension of cclm should be nloat'
       if (size(aalm,2).NE.size(cclm,2)) print*, 'ERROR lomain: nbands_dft should be dimension in aalm and cclm', size(aalm,3)
       nbands_dft = size(aalm,2)
       if (nbands_dft.ne.(nemax-nemin+1)) print*, 'ERROR lomain: nbands_dft!=nemax-nemin+1 !'
    endif

    TWOPI=TWO*PI
    nmult = maxval(mult)
    ALLOCATE( iind((lomax+1)*(lomax+1),nloat,nmult) )
    iind(:,:,:)=0
    i=nlov
    DO L=0,LoMAX
       do jlo=1,ilo(l)
          do jneq=1,mult(jatom)
             DO M1=-l,l
                i=i+1
                iind(l*(l+1)+m1+1,jlo,jneq)=i
             END DO
          end do
       end do
    END DO

    rotloc_x_BR1(:,:) = matmul(rotloc(:,:,jatom), BR1)
    rotloc_x_BR1_x_rotij = matmul( rotloc_x_BR1, rotij(:,:,latom) )
    xrotloc_x_BR1(:,:) = matmul(xrotloc, BR1)
    
    !$OMP PARALLEL DO PRIVATE(jlo,jneq,M,M1,i,bkrot,bkrot2,bkrloc,Krot,yl,arg1,argt2,arg2,phshel,num,phs,index,yp,ayp,byp,cyp)&
    !$OMP& SHARED(alm,blm,clm,aalm,bblm,cclm,alo,blo,clo)&
    !$OMP& SCHEDULE(STATIC)
    DO l=0,min(lomax,lmax2)
       do jlo=1,ilo(l)
          do jneq=1,mult(jatom)
             DO m1=-l,+l
                i = iind(l*(l+1)+m1+1,jlo,jneq)
                bkrot(:)=BK3lo(:,i)
                ! BKROT2 = R_n.R_a.(k+K), transformation from the first atom to an equivalent atom
                bkrot2 = matmul(rotij(:,:,latom), bkrot)
                !---- BR1 transforms integer reciprocal lattice vectors, as given in the VECTORLIST of LAPW1, into cartesian system ----!
                ! BKROT3 = R_n.R_a.(k+K), but in cartesian coordinate system
                !bkrot3 = matmul(BR1, bkrot2)
                !---- BKRLOC = crotloc.R_n.R_a.(k+K),  rotates according to the user specified local coordinate system.
                !bkrloc = matmul(XROTLOC, bkrot3)
                bkrloc = matmul( xrotloc_x_BR1, bkrot2)
                !---- YLM = Y_{L}(Rotloc.R_g.(k+K))
                CALL YLM(BKRLOC,lomax,YL,(lomax+1)*(lomax+1))
                ! (R_n.R_a.(k+K)) *  R(first) * 2pi
                ARG1 = dot_product( BKROT2, POS(:,LFIRST) )*TWOPI
                ! ARGT2 = (R_a.(k+K)).tau_n * 2pi
                ARGT2= dot_product(BKROT,TAUIJ(:,LATOM))*TWOPI
                !
                ARG2 = dot_product(BKROT,shft(latom,:))*TWOPI
                ! PHSEHL = e^{I*2pi*( (R_a.(k+K))*tau_n + (K+k)*tau(isym) + (R_n.R_a.(k+K)*R(first)))}
                PHSHEL=EXP( IMAG*(ARG1+ARG2+ARGT2) )
                !print *, 'p=', ARG1, ARG2, ARG3, ARGT
                DO NUM=NEMIN,NEMAX
                   PHS(NUM) = PHSHEL*As_lo(I,NUM,is)   
                ENDDO
                if(Qforce_j) then
                   !bkrot2 = matmul(rotij(:,:,latom), k3lo(:,i))
                   !bkrot3 = matmul(BR1, bkrot2)
                   !Krot = matmul(rotloc(:,:,jatom), bkrot3) ! Note that the force should be computed in local coordinate system specified in the case.struct input file. All other forces are computed in this coordinate system!
                   Krot = matmul( rotloc_x_BR1_x_rotij, k3lo(:,i) ) ! Note that the force should be computed in local coordinate system specified in the case.struct input file. All other forces are computed in this coordinate system!
                endif
                do num=nemin,nemax
                   DO m=-l,+l
                      index=l*(l+1)+m+1
                      yp = conjg(YL(index))*PHS(num)
                      ayp = Alo(l,jlo)*yp
                      byp = Blo(l,jlo)*yp
                      cyp = Clo(l,jlo)*yp
                      
                      alm(index,num)     = alm(index,num)     + ayp
                      blm(index,num)     = blm(index,num)     + byp
                      clm(index,num,jlo) = clm(index,num,jlo) + cyp
                      
                      if (Qforce_j) then
                         aalm(index,num-nemin+1,:)     = aalm(index,num-nemin+1,:)     + ayp*Krot(:)
                         bblm(index,num-nemin+1,:)     = bblm(index,num-nemin+1,:)     + byp*Krot(:)
                         cclm(index,num-nemin+1,jlo,:) = cclm(index,num-nemin+1,jlo,:) + cyp*Krot(:)
                      endif
                   ENDDO
                end do
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    DEALLOCATE( iind )
    deallocate( YL, phs )
    
    return
  END SUBROUTINE lomain
END module mod_lo
