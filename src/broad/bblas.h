// @Copyright 2007 Kristjan Haule
// 
#ifndef _SBLAS
#define _SBLAS
#include <string>
#include <complex>
//#include "complex.h"
//#include "sutil.h"

#ifdef NO_APPEND_FORTRAN
# define FNAME(x) x
#else
# define FNAME(x) x##_
#endif


typedef void (*compfunc)(int* b, std::complex<double>* w, int* N);
typedef void (*compfund)(int* b, double* wr, double* wi, int* N);

extern "C" {
  void FNAME(dptsv)(const int* N, const int* NRHS, double* D, double* E, double* B, const int* LDB, int* INFO);
  void FNAME(zptsv)(const int* N, const int* NRHS, double* D, std::complex<double>* E, std::complex<double>* B, const int* LDB, int* INFO);
  void FNAME(dgetrf)(int* n1, int* n2, double* a, int* lda, int* ipiv,int* info);
  void FNAME(zgetrf)(int* n1, int* n2, std::complex<double>* a, int* lda, int* ipiv,int* info);
  void FNAME(dgetrs)(const char* trans, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
  void FNAME(zgetrs)(const char* trans, int* n, int* nrhs, std::complex<double>* a, int* lda, int* ipiv, std::complex<double>* b, int* ldb, int* info);
  void FNAME(zgemm)(const char* transa, const char* transb, const int* m, const int* n, const int* k, const std::complex<double>* alpha, const std::complex<double>* A, const int* lda, const std::complex<double>* B, const int* ldb, const std::complex<double>* beta, std::complex<double>* C, const int* ldc);
  void FNAME(dgemm)(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha, const double* A, const int* lda, const double* B, const int* ldb, const double* beta, double* C, const int* ldc);
  void FNAME(dsymm)(const char* side, const char* uplo, const int* m, const int* n, const double* alpha, const double* A, const int* lda, const double* B, const int* ldb, const double* beta, double* C, const int* ldc);
  void FNAME(dgeev)(const char* jobvl, const char* jobvr,  const int* n,  double* A, const int* lda, double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr, double* work, const int* lwork, int* info);
  void FNAME(dsyev)(const char* jobz,  const char* uplo,   const int* n,  double* A, const int* lda, double* w, double* work, const int* lwork, int* info);
  void FNAME(zgesdd)(const char* jobz, const int* m, const int* n, std::complex<double>* A, const int* lda, double* S, std::complex<double>* U, const int* ldu, std::complex<double>* Vt, const int* ldvt, std::complex<double>* work, const int* lwork, double* rwork, int* iwork, int* info);
  void FNAME(dgesdd)(const char* jobz, const int* m, const int* n, double* A, const int* lda, double* S, double* U, const int* ldu, double* Vt, const int* ldvt, double* work, const int* lwork, int* iwork, int* info);
  void FNAME(zgeev)(const char* jobvl, const char* jobvr,  const int* n,  std::complex<double>* A,const int* lda, std::complex<double>* w, std::complex<double>* vl, const int* ldvl, std::complex<double>* vr, const int* ldvr, std::complex<double>* work, const int* lwork, double* rwork, int* info);
  void FNAME(zheevd)(const char* job, const char* uplo,  const int* n,  std::complex<double>* A,const int* lda, double* w, std::complex<double>* work, const int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info);
  void FNAME(mzgees)(const char* jobvs, const char* sort, compfunc select, int* N, std::complex<double>* A, int* lda, int* sdim, std::complex<double>* W, std::complex<double>* VS, int* ldvs,  std::complex<double>* work,  int* lwork,  double* rwork,  int* bwork, int* info);
  void FNAME(mdgees)(const char* jobvs, const char* sort, compfund select, int* N, double* A, int* lda, int* sdim, double* wr, double* wi, double* VS, int* ldvs,  double* work,  int* lwork, int* bwork, int* info);
  double FNAME(dznrm2)(const int* N, const std::complex<double>* x, const int* incx);
  double FNAME(dnrm2)(const int* N, const double* x, const int* incx);
  void FNAME(zdscal)(const int* N, const double* alpha, std::complex<double>* zx, const int* incx);
  void FNAME(dscal)(const int* N, const double* alpha, double* zx, const int* incx);
  void FNAME(zgemv)(const char* trans, const int* m, const int* n, const std::complex<double>* alpha, const std::complex<double>* A, int* lda, const std::complex<double>* x, const int* incx, const std::complex<double>* beta, std::complex<double>* y, const int* incy);
  void FNAME(dgemv)(const char* trans, const int* m, const int* n, const double* alpha, const double* A, int* lda, const double* x, const int* incx, const double* beta, double* y, const int* incy);
  double FNAME(ddot)(const int* N, double* x, const int* incx, double* y, const int* incy);
  void FNAME(zsytri)(const char* uplo, const int* N, std::complex<double>* A, const int* lda, int* ipiv, std::complex<double>* work, int* info);
  void FNAME(zsytrf)(const char* uplo, const int* N, std::complex<double>* A, const int* lda, int* ipiv, std::complex<double>* work, int* lwork, int* info);
  void FNAME(dsytri)(const char* uplo, const int* N, double* A, const int* lda, int* ipiv, double* work, int* info);
  void FNAME(dsytrf)(const char* uplo, const int* N, double* A, const int* lda, int* ipiv, double* work, int* lwork, int* info);
  void FNAME(zaxpy)(const int* N, const std::complex<double>* alpha, const std::complex<double>* x, const int* incx, std::complex<double>* y, const int* incy);
  _Complex double FNAME(zdotc)(const int* N, const std::complex<double>* zx, const int* incx, const std::complex<double>* zy, const int* incy);
  //std::complex<double> zdotc_(const int* N, const std::complex<double>* zx, const int* incx, const std::complex<double>* zy, const int* incy);
  void FNAME(dger)(const int* n1, const int* n2, const double* alpha, const double* x, const int* incx, const double* y, const int* incy, double* A, const int* lda);  
}

inline void xptsv_(const int* N, const int* NRHS, double* D, double* E, double* B, const int* LDB, int* INFO)
{  FNAME(dptsv)(N, NRHS, D, E, B, LDB, INFO);}

inline void xptsv_(const int* N, const int* NRHS, double* D, std::complex<double>* E, std::complex<double>* B, const int* LDB, int* INFO)
{  FNAME(zptsv)(N, NRHS, D, E, B, LDB, INFO);}

inline void xgemm(const std::string& transa, const std::string& transb, const int m, const int n,
		  const int k, const double alpha, const double* A,
		  const int lda, const double* B, const int ldb, const double beta,
		  double* C, const int ldc)
{
  /*
  cout<<"m="<<m<<" n="<<n<<" k="<<k<<" lda="<<lda<<" ldb="<<ldb<<" ldc="<<ldc;
  cout<<" transa="<<transa<<" transb="<<transb<<endl;
  cout<<"A[0,0]="<<A[0]<<" "<<A[1]<<" "<<A[2]<<" "<<A[3]<<endl;
  cout<<"B[0,0]="<<B[0]<<" "<<B[1]<<" "<<B[2]<<" "<<B[3]<<endl;
  */
  FNAME(dgemm)(transa.c_str(), transb.c_str(), &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

inline void xsymm(const char* side, const char* uplo, int m, int n, double alpha, const double* A, int lda, const double* B, int ldb, double beta, double* C, int ldc)
{
  FNAME(dsymm)(side, uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

inline void xgemm(const std::string& transa, const std::string& transb, const int m, const int n,
		  const int k, const std::complex<double>& alpha, const std::complex<double>* A,
		  const int lda, const std::complex<double>* B, const int ldb, const std::complex<double>& beta,
		  std::complex<double>* C, const int ldc)
{
  FNAME(zgemm)(transa.c_str(), transb.c_str(), &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

inline void xaxpy(int n, const std::complex<double>& alpha, const std::complex<double>* x, std::complex<double>* y, int incx=1, int incy=1)
{ FNAME(zaxpy)(&n, &alpha, x, &incx, y, &incy);}

//inline std::complex<double> xdotc(int n, const std::complex<double>* zx, const std::complex<double>* zy, int incx=1, int incy=1)
//{ return reinterpret_cast<std::complex<double>>(zdotc_(&n, zx, &incx, zy, &incy));}

inline int xgetrf(int n, double* a, int lda, int* ipiv)
{
  int info = 0;
  FNAME(dgetrf)(&n, &n, a, &lda, ipiv, &info);
  if (info){
    std::cerr << "Something wrong in LU (real) decomposition! " << info << std::endl;
  }  
  return info;
}

inline int xgetrf(int n, std::complex<double>* a, int lda, int* ipiv)
{
  int info = 0;
  FNAME(zgetrf)(&n, &n, a, &lda, ipiv, &info);
  if (info){
    std::cerr << "Something wrong in LU (complex) decomposition! " << info << std::endl;
  }  
  return info;
}

inline int xgetrs(int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb)
{
  int info = 0;
  FNAME(dgetrs)("T", &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  if (info){
    std::cerr << "Something wrong with the system of (real) equations! " << info << std::endl;
  }  
  return info;
}

inline int xgetrs(int n, int nrhs, std::complex<double>* a, int lda, int* ipiv, std::complex<double>* b, int ldb)
{
  int info = 0;
  FNAME(zgetrs)("T", &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  if (info){
    std::cerr << "Something wrong with the system of (complex) equations! " << info << std::endl;
  }
  return info;
}

template <class T>
inline int xgesv(int n, int nrhs, T* a, int lda, int* ipiv, T* b, int ldb)
{
  int info = xgetrf(n, a, lda, ipiv);
  if (info) return info;
  return xgetrs(n, nrhs, a, lda, ipiv, b, ldb);
}

// Eigenvalues of real symmetric matrix
inline int xsyev(const int N, double* A, const int lda, double* w, double* work, const int lwork, std::string job="N")
{
  int info = 0;
  FNAME(dsyev)(job.c_str(), "U",  &N,  A, &lda, w, work, &lwork, &info);
  if (info)
    std::cerr << "Can't compute eigenvalues! " << info << std::endl;
  return info;
}

inline int xgeev(const int N, double* A, const int lda, std::complex<double>* w, double* wr, double* wi, double* work,
		 const int lwork, double*)
{
  int ena = 1, info = 0;
  FNAME(dgeev)("N", "N",  &N,  A, &lda, wr, wi, NULL, &ena, NULL, &ena, work, &lwork, &info);
  if (info)
    std::cerr << "Can't compute eigenvalues! " << info << std::endl;
  for (int i=0; i<N; i++) {
    //w[i].real()=wr[i]; w[i].imag()=wi[i];
    w[i] = std::complex<double>(wr[i],wi[i]);
  }
  return info;
}

inline int xgeev_(const int N, double* A, const int lda, double* wr, double* wi, double* work,
		 const int lwork)
{
  int ena = 1, info = 0;
  FNAME(dgeev)("N", "N",  &N,  A, &lda, wr, wi, NULL, &ena, NULL, &ena, work, &lwork, &info);
  if (info)
    std::cerr << "Can't compute eigenvalues! " << info << std::endl;
  return info;
}

inline int xgeev(int N, std::complex<double>* A, int lda, std::complex<double>* w, std::complex<double>* AL, int ld_AL,
		 std::complex<double>* AR, int ld_AR,  std::complex<double>* work, int lwork, double* rwork)
{
  int info = 0;
  FNAME(zgeev)("V", "V",  &N,  A, &lda, w, AL, &ld_AL, AR, &ld_AR, work, &lwork, rwork, &info);
  if (info) std::cerr << "Can't compute eigenvalues and eigenvectors! " << info << std::endl;
  return info;
}

inline int xheevd(int N, std::complex<double>* A, int lda, double* w, std::complex<double>* work, int lwork, double* rwork, int lrwork,
		  int* iwork, int liwork)
{
  int info=0;
  FNAME(zheevd)("V", "U",  &N,  A, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  if (info) std::cerr << "Can't compute eigenvalues and eigenvectors! " << info << std::endl;
  return info;
}

inline int xgesdd(bool vect, int M, int N, std::complex<double>* A, int lda, double* S, std::complex<double>* U, int ldu, std::complex<double>* Vt, int ldvt,
		  std::complex<double>* work, int lwork, double* rwork, int *iwork)
{
  int info = 0;
  std::string job = vect ? "S" : "N";
#ifndef _HP_  
  FNAME(zgesdd)(job.c_str(), &M, &N, A, &lda, S, U, &ldu, Vt, &ldvt, work, &lwork, rwork, iwork, &info);
#endif  
  if (info) {
    std::cerr << "Can't compute SVD of the kernel! " << info << std::endl;
  }
  return info;
}

inline int xgesdd(bool vect, int M, int N, double* A, int lda, double* S, double* U, int ldu, double* Vt, int ldvt,
		  double* work, int lwork, double*, int *iwork)
{
  int info = 0;
  std::string job = vect ? "S" : "N";
#ifndef _HP_  
  FNAME(dgesdd)(job.c_str(), &M, &N, A, &lda, S, U, &ldu, Vt, &ldvt, work, &lwork, iwork, &info);
#endif  
  if (info) {
    std::cerr << "Can't compute SVD of the kernel! " << info << std::endl;
  }
  return info;
}

inline void comparec(int* b, std::complex<double>* w, int* N)
{
  double dmin = 1e10;
  int imin = 0;
  for (int i=0; i<(*N); i++){
    if (fabs(w[i].real()) < dmin){
      dmin = fabs(w[i].real());
      imin = i;
    }
  }
  for (int i=0; i<(*N); i++) b[i] = 0;
  b[imin] = 1;
}

inline void compared(int* b, double* wr, double* wi, int* N)
{
  double dmin = 1e10;
  int imin = 0;
  for (int i=0; i<(*N); i++){
    if (fabs(wr[i]) < dmin){
      dmin = fabs(wr[i]);
      imin = i;
    }
  }
  for (int i=0; i<(*N); i++) b[i] = 0;
  b[imin] = 1;
}

inline int xgees(bool vect, int N, std::complex<double>* A, int lda, std::complex<double>* w, double*, double*, std::complex<double>* Z, int ldz, std::complex<double>* work,
		 int lwork, double* rwork, int* bwork)
{
  int sdim, info = 0;
  std::string job = vect ? "V" : "N";
  FNAME(mzgees)(job.c_str(), "S", &comparec, &N, A, &lda, &sdim, w, Z, &ldz, work, &lwork, rwork, bwork, &info);
  if (info) std::cerr << "Can't perform Schur decomposition! " << info << std::endl;
  return info;
}

inline int xgees(bool vect, int N, double* A, int lda, std::complex<double>* w, double* wr, double* wi, double* Z, int ldz, double* work,
		 int lwork, double*, int* bwork)
{
  int sdim, info = 0;
  std::string job = vect ? "V" : "N";
  FNAME(mdgees)(job.c_str(), "S", &compared, &N, A, &lda, &sdim, wr, wi, Z, &ldz, work, &lwork, bwork, &info);
  if (info) std::cerr << "Can't perform Schur decomposition! " << info << std::endl;
  for (int i=0; i<N; i++) {
    //w[i].real()=wr[i]; w[i].imag()=wi[i];
    w[i] = std::complex<double>(wr[i],wi[i]);
  }
  return info;
}

// template <enum TypeOfMatrix TN>
// inline void xgemv(int m, int n, const std::complex<double>& alpha, const std::complex<double>* A, int lda,
// 		  const std::complex<double>* x, const std::complex<double>& beta, std::complex<double>* y)

// {
//   std::cerr << "Not implemented yet! " << std::endl;
// }

// template <>
// inline void xgemv<_Normal>(int m, int n, const std::complex<double>& alpha, const std::complex<double>* A, int lda,
// 			  const std::complex<double>* x, const std::complex<double>& beta, std::complex<double>* y)
// {
//   int inc = 1;
//   zgemv_("T", &m, &n, &alpha, A, &lda, x, &inc, &beta, y, &inc);
// }

// template <>
// inline void xgemv<_Transpose>(int m, int n, const std::complex<double>& alpha, const std::complex<double>* A, int lda,
// 			     const std::complex<double>* x, const std::complex<double>& beta, std::complex<double>* y)
// {
//   int inc = 1;
//   zgemv_("N", &m, &n, &alpha, A, &lda, x, &inc, &beta, y, &inc);
// }

// template <enum TypeOfMatrix TN>
// inline void xgemv(int m, int n, double alpha, const double* A, int lda,
// 		  const double* x, const double beta, double* y)
// {
//   std::cerr << "Not implemented yet! " << std::endl;
// }

// template <>
// inline void xgemv<_Normal>(int m, int n, double alpha, const double* A, int lda,
// 			  const double* x, const double beta, double* y)
// {
//   int inc = 1;
//   dgemv_("T", &m, &n, &alpha, A, &lda, x, &inc, &beta, y, &inc);
// }

// template <>
// inline void xgemv<_Transpose>(int m, int n, double alpha, const double* A, int lda,
// 			     const double* x, const double beta, double* y)
// {
//   int inc = 1;
//   dgemv_("N", &m, &n, &alpha, A, &lda, x, &inc, &beta, y, &inc);
// }

inline double xnrm(int N, const std::complex<double>* x, int incx = 1)
{
  return FNAME(dznrm2)(&N, x, &incx);
}

inline double xnrm(int N, const double* x, int incx = 1)
{
  return FNAME(dnrm2)(&N, x, &incx);
}

inline void xdscal(int N, double alpha, std::complex<double>* zx, int incx = 1)
{
  FNAME(zdscal)(&N, &alpha, zx, &incx);
}

inline void xdscal(int N, double alpha, double* zx, int incx = 1)
{
  FNAME(dscal)(&N, &alpha, zx, &incx);
}

inline void xsytrf(const char* uplo, int N, std::complex<double>* A, int lda, int* ipiv, std::complex<double>* work, int lwork, int& info)
{
  FNAME(zsytrf)(uplo, &N, A, &lda, ipiv, work, &lwork, &info);
}
inline void xsytrf(const char* uplo, int N, double* A, int lda, int* ipiv, double* work, int lwork, int& info)
{
  FNAME(dsytrf)(uplo, &N, A, &lda, ipiv, work, &lwork, &info);
}
inline void xsytri(const char* uplo, int N, std::complex<double>* A, int lda, int* ipiv, std::complex<double>* work, int& info)
{
  FNAME(zsytri)(uplo, &N, A, &lda, ipiv, work, &info);
}
inline void xsytri(const char* uplo, int N, double* A, int lda, int* ipiv, double* work, int& info)
{
  FNAME(dsytri)(uplo, &N, A, &lda, ipiv, work, &info);
}

inline int xgeev(const int N, std::complex<double>* A, const int lda, std::complex<double>* w, std::complex<double>* work,
		 const int lwork, double* rwork)
{
  int ena = 1, info = 0;
  FNAME(zgeev)("N", "N",  &N,  A, &lda, w, NULL, &ena, NULL, &ena, work, &lwork, rwork, &info);
  if (info)
    std::cerr << "Can't compute eigenvalues! " << info << std::endl;
  return info;
}

#endif //_SBLAS
