subroutine diracout(dp,dq,rel,V,rnot,dstep,nmax,EH,nql,nqk,val,slo,nodes,z,nrad)
  !
  !         Integration of Dirac equation.
  ! 
  !  Input:
  ! 
  !    rel    switch for relativistic - nonrelativ. calculation
  !    V      rad.sym. potential in Hartree's
  !    rnot   first radial meshpoint
  !    dstep  log step for radial mesh
  !    nmax   number of radial meshpoints
  !    EH     energy in hartree
  !    nql    angular momentum 
  !    nqk    relativistic quantum number kappa 
  !    z      charge of nucleus
  !
  !  Output: 
  !
  !    val,slo:  Wave function and slope at the MT-boundary
  !    nodes:    number of nodes
  !
  !
  ! ----------------------------------------------------------------
  !USE defs,  ONLY: clight, test
  !USE param, ONLY: nrad
  !USE ams, ONLY: atom_mass
  !implicit real*8 (a-h,o-z)
  IMPLICIT NONE
  ! dp    =  large component of the solution of the dirac equation
  ! dq    =  small component of the solution
  REAL*8, intent(out)   :: val, slo
  INTEGER, intent(out)  :: nodes
  REAL*8, intent(out)   :: dp(nrad),dq(nrad)
  LOGICAL, intent(in)   :: rel
  REAL*8, intent(in)    :: V(NRAD)             ! V = potential*r
  REAL*8, intent(in)    :: rnot, dstep, EH, Z
  INTEGER, intent(in)   :: nmax, nql, nqk, nrad
  !DIMENSION  AP(NRAD),BP(NRAD),AE(NRAD),BE(NRAD)
  !COMMON  /uhelp/     dp,dq,AP,BP,AE,BE
  !SAVE    /uhelp/
  !
  ! DEP,DEQ DERIVEES DE DP ET DQ   DB=ENERGIE/DVC    DVC VITESSE DE LA
  ! LUMIERE EN U.A.   DSAL=2.*DVC   DK NOMBRE QUANTIQUE KAPPA
  ! DM=PAS EXPONENTIEL/720., DKOEF=1./720.

  !  The name of dm should be changed to avoid a name collision
  !  with dm in inouh and inth
  !COMMON /PS1/ DEP(5),DEQ(5),DB,DVC,DSAL,DK,DM  ! exchanges data with inouh and 
  !save   /ps1/
  !REAL*8 :: dep, deq, db, dvc, dsal, dk, dm
  !
  REAL*8,PARAMETER       :: CLIGHT= 137.0359895d0
  REAL*8,PARAMETER       :: TEST=   1.D-12
  !
  REAL*8 :: dep(5), deq(5), db, dvc, dsal, dk, dm
  DATA DKOEF/.1388888888888888D-2/
  ! locals
  REAL*8  :: dv(nrad)
  REAL*8  :: dr(NRAD) ! radial mesh
  REAL*8  :: d1, dfl, dkoef, dq1, dval, rnuc
  INTEGER :: i, nstop, nuc
  !
  REAL*8          :: atom_mass(103)
  atom_mass=[1.,4.,6.9,9.,10.8,12.,14.,16.,19.,20.2, &
       23.,24.3,27.,28.1,31.,32.,35.4,40.,39.1,40.,45., &
       47.9,50.9,52.,54.9,55.8,58.9,58.7,63.5,65.4,69.7, &
       72.6,74.9,79.,79.9,83.8,85.5,87.6,88.9,91.2,92.9, &
       95.9,98.,101.1,102.9,106.4,107.9,112.4,114.8, &
       118.7,121.8,127.6,126.9,131.3,132.9,137.3,138.9,140.1, &
       140.9,144.2,145.,150.4,152.,157.3,158.9,162.5, &
       164.9,167.3,168.9,173.,175.,178.5,180.9,183.8,186.2, &
       190.2,192.2,195.1,197.,200.6,204.4,207.2,209.,209., &
       210.,222.,223.,226.,227.,232.,231.,238.,237.,244.,243., &
       247.,247.,251.,252.,257.,258.,259.,262.]
  
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

  !jk   finite size of the nucleus
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
  !CALL INOUH (dp,dq,dr,dq1,dfl,dv(1),Z,TEST,nuc,NSTOP)
  CALL inouh (dp,dq, DEP,DEQ, dr,dq1,dfl,dv(1),Z,test,nuc,nstop,  DB,DVC,DSAL,DK, nrad)
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
     !call inth (dp(i),dq(i),dv(i),dr(i))
     CALL INTH(dp(i),dq(i),dep,deq,db,dvc,dsal,dk,dm,dv(i),dr(i))
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
END subroutine diracout

