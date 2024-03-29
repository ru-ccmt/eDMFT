subroutine diracout_old(rel,v,rnot,dstep,nmax,eh,nql,nqk,val,slo,nodes,z)
  !         Integration of Dirac equation.
  ! 
  !  Input:
  ! 
  !    rel    switch for relativ. - nonrelativ. calculation
  !    v      rad.sym. potential in Hartree
  !    rnot   first radial meshpoint
  !    dstep  log. step
  !    nmax    number of radial meshpoints
  !    eh     energy in hartree
  !    nql    angular momentum 
  !    nqk    relativistic quantum number kappa 
  !    z      charge of nucleus
  !
  !  Output: 
  !
  !    val,slo:  Wellenfunktion und Steigung am Kugelrand
  !    nodes:    nomber of nodes
  !
  ! ----------------------------------------------------------------
  USE param, ONLY: nrad, clight
  USE ams, ONLY: atom_mass
  IMPLICIT NONE
  logical, intent(in) :: rel
  REAL*8,  intent(in) :: V(NRAD), Rnot, dstep, Eh, Z   ! V = potential*r
  INTEGER, intent(in) :: nmax, nql, nqk
  REAL*8,  intent(out):: val, slo
  INTEGER, intent(out):: nodes
  ! common blocks
  REAL*8 :: dp(nrad),dq(nrad) !  dp/dq    =  large/small component of the solution of the dirac equation
  REAL*8 :: AP(NRAD),BP(NRAD),AE(NRAD),BE(NRAD)
  COMMON  /work/     dp,dq,AP,BP,AE,BE
  SAVE    /work/
  REAL*8 :: DEP(5), DEQ(5), DB, DVC, DSAL, DK, DM
  COMMON /PS1/ DEP,DEQ,DB,DVC,DSAL,DK,DM
  save   /ps1/
  ! locals
  REAL*8, PARAMETER :: DKOEF = 0.1388888888888888D-2  ! 1/720
  REAL*8 :: dr(NRAD), dv(nrad)  ! radial mesh, V/r
  REAL*8 :: d1, dfl, dq1, dval, rnuc, test
  INTEGER:: i, nstop, nuc
  !
  !
  ! DEP,DEQ DERIVEES DE DP ET DQ   DB=ENERGIE/DVC    DVC VITESSE DE LA
  ! LUMIERE EN U.A.   DSAL=2.*DVC   DK NOMBRE QUANTIQUE KAPPA
  ! DM=PAS EXPONENTIEL/720., DKOEF=1./720.
  
  !  The name of dm should be changed to avoid a name collision
  !  with dm in inouh
  
  !   Set up radial mesh.
  do i=1,nmax
     DR(i)=RNOT*(exp(DSTEP*(i-1.d0)))
  enddo
  
  if (rel) then
     dvc = clight
  else 
     dvc = 1.d+10
  endif
  dsal = 2.d0*dvc
  db = eh/dvc
  dk = nqk
  dm=dstep*dkoef
  
  do i=1,nmax
     dv(i) = v(i)/dr(i)
  enddo
  !  Behavior of the solution at the origin
  
  !  finite size of the nucleus
  rnuc=2.2677D-05*(atom_mass(int(z))**(1/3.))
  do i=1,nmax
     d1=rnot*exp(DSTEP*(i-1.d0))
     if (d1.ge.rnuc) exit
  enddo
  nuc=I
  if (nuc.le.0) then
     dfl = sqrt(nqk*nqk-z*z/(dvc*dvc))
  else
     dfl=nqk*nqk
     do i=1,nuc
        dv(i)=dv(i)+z/dr(i)+z*((dr(i)/dr(nuc))**2-3.)/(2*dr(nuc))
     enddo
  end if
  dq1 = nqk/iabs(nqk)
  
  !  Determine expansion of the potential at the origin.
  test =1.e-8
  CALL INOUH_old (dp,dq,dr,dq1,dfl,dv(1),Z,TEST,nuc,NSTOP)

  !  Set up boundary conditions using a small r expansion.
  nodes = 0
  do i=1,5
     dval=dr(i)**dfl
     if (i.ne.1) then
        if (dp(i-1).ne.0.) then
           if ((dp(i)/dp(i-1)).le.0.) then
              nodes=nodes+1
           endif
        endif
     endif
     dp(i) = dp(i)*dval
     dq(i) = dq(i)*dval
     dep(i)=dep(i)*dval
     deq(i)=deq(i)*dval
  enddo
  
  !    Perform outward integration of dirac equation
  do i = 6, nmax
     dp(i) = dp(i-1)
     dq(i) = dq(i-1)
     call inth_old (dp(i),dq(i),dv(i),dr(i))
     if (dp(i-1).ne.0.) then
        if ((dp(i)/dp(i-1)).gt.0.) then
           nodes=nodes+1
        endif
     endif
  enddo
  val = dp(nmax)/dr(nmax)
  slo = dep(5)/(dstep*dr(nmax))/dvc*2.d0
  slo = (slo-val)/dr(nmax) 
  RETURN
END subroutine diracout_old
