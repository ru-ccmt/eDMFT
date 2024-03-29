REAL*8 FUNCTION GAUNT1(LP,L,LS,MP,M,MS)
  IMPLICIT NONE
  !
  !        Arguments
  !
  INTEGER, intent(in) :: L, LP, LS, M, MP, MS
  !..................................................................
  !
  !   GAUNT computes the integral of
  !      CONJG(Y(LP,MP))*Y(L,M)*Y(LS,MS) for LP+L+LS .LE. 23
  !   using gaussian quadrature with N=12 as given by
  !      M. Abramowitz and I.A. Stegun,
  !      'Handbook of Mathematical Functions',
  !      NBS Applied Mathematics Series 55 (1968), pages 887 and 916
  !   written by Bruce Harmon based on suggestion by W. Rudge
  !   Iowa State Sept.1973
  !
  !   extended by M. Weinert and E. Wimmer
  !   Northwestern University March 1980
  !
  !..................................................................
  INTEGER            MAXDIM, N
  PARAMETER          (MAXDIM = 81, N = 6)
  DOUBLE PRECISION   ZERO
  PARAMETER          (ZERO = 0.0D+0)
  !        Common blocks
  DOUBLE PRECISION   YR(N,MAXDIM)
  COMMON  /ASSLEG/   YR
  SAVE    /ASSLEG/
  !
  !        Local Scalars
  !
  INTEGER            I, IL, ILP, ILS
  DOUBLE PRECISION   S
  DOUBLE PRECISION   W(N)
  !
  !        Data statements
  !
  DATA W /0.24914704581340D+0, 0.23349253653836D+0,0.20316742672307D+0, 0.16007832854335D+0,0.10693932599532D+0, 0.04717533638651D+0/
  !
  IL = L*(L+1) + M + 1
  ILP = LP*(LP+1) + MP + 1
  ILS = LS*(LS+1) + MS + 1
  S = ZERO
  DO I = 1, N
     S = S + W(I)*YR(I,ILP)*YR(I,IL)*YR(I,ILS)
  ENDDO
  GAUNT1 = S
  !
  RETURN
END FUNCTION GAUNT1
