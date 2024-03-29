// @Copyright 2007 Kristjan Haule
// 
#ifndef _BLAS_
#define _BLAS_
#include <string>
#include "complex.h"
#include "util.h"

#ifdef _HP_
#define dgetrf_  dgetrf
#define zgetrf_	 zgetrf
#define dgetrs_	 dgetrs
#define zgetrs_	 zgetrs
#define zgemm_	 zgemm
#define dgemm_	 dgemm
#define dgeev_	 dgeev
#define zgesdd_	 zgesdd
#define dgesdd_	 dgesdd
#define zgeev_	 zgeev
#define mzgees_	 mzgees
#define mdgees_	 mdgees
#define dznrm2_	 dznrm2
#define dnrm2_	 dnrm2
#define zdscal_	 zdscal
#define dscal_	 dscal
#define zgemv_	 zgemv
#define dgemv_	 dgemv
#endif

typedef void (*compfunc)(int* b, dcomplex* w, int* N);
typedef void (*compfund)(int* b, double* wr, double* wi, int* N);

extern "C" {
  void dgetrf_(int* n1, int* n2, double* a, int* lda, int* ipiv,int* info);
  void zgetrf_(int* n1, int* n2, dcomplex* a, int* lda, int* ipiv,int* info);
  void dgetrs_(const char* trans, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
  void zgetrs_(const char* trans, int* n, int* nrhs, dcomplex* a, int* lda, int* ipiv, dcomplex* b, int* ldb, int* info);
  void zgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, const dcomplex* alpha, const dcomplex* A, const int* lda, const dcomplex* B, const int* ldb, const dcomplex* beta, dcomplex* C, const int* ldc);
  void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha, const double* A, const int* lda, const double* B, const int* ldb, const double* beta, double* C, const int* ldc);
  void dsymm_(const char* side, const char* uplo, const int* m, const int* n, const double* alpha, const double* A, const int* lda, const double* B, const int* ldb, const double* beta, double* C, const int* ldc);
  void dgeev_(const char* jobvl, const char* jobvr,  const int* n,  double* A, const int* lda, double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr, double* work, const int* lwork, int* info);
  void dsyev_(const char* jobz,  const char* uplo,   const int* n,  double* A, const int* lda, double* w, double* work, const int* lwork, int* info);
  void zgesdd_(const char* jobz, const int* m, const int* n, dcomplex* A, const int* lda, double* S, dcomplex* U, const int* ldu, dcomplex* Vt, const int* ldvt, dcomplex* work, const int* lwork, double* rwork, int* iwork, int* info);
  void dgesdd_(const char* jobz, const int* m, const int* n, double* A, const int* lda, double* S, double* U, const int* ldu, double* Vt, const int* ldvt, double* work, const int* lwork, int* iwork, int* info);
  void  zgeev_(const char* jobvl, const char* jobvr,  const int* n,  dcomplex* A,const int* lda, dcomplex* w, dcomplex* vl, const int* ldvl, dcomplex* vr, const int* ldvr, dcomplex* work, const int* lwork, double* rwork, int* info);
  void mzgees_(const char* jobvs, const char* sort, compfunc select, int* N, dcomplex* A, int* lda, int* sdim, dcomplex* W, dcomplex* VS, int* ldvs,  dcomplex* work,  int* lwork,  double* rwork,  int* bwork, int* info);
  void mdgees_(const char* jobvs, const char* sort, compfund select, int* N, double* A, int* lda, int* sdim, double* wr, double* wi, double* VS, int* ldvs,  double* work,  int* lwork, int* bwork, int* info);
  double dznrm2_(const int* N, const dcomplex* x, const int* incx);
  double dnrm2_(const int* N, const double* x, const int* incx);
  void zdscal_(const int* N, const double* alpha, dcomplex* zx, const int* incx);
  void dscal_(const int* N, const double* alpha, double* zx, const int* incx);
  void zgemv_(const char* trans, const int* m, const int* n, const dcomplex* alpha, const dcomplex* A, int* lda, const dcomplex* x, const int* incx, const dcomplex* beta, dcomplex* y, const int* incy);
  void dgemv_(const char* trans, const int* m, const int* n, const double* alpha, const double* A, int* lda, const double* x, const int* incx, const double* beta, double* y, const int* incy);
  double ddot_(const int* N, double* x, const int* incx, double* y, const int* incy);
  void zsytri_(const char* uplo, const int* N, dcomplex* A, const int* lda, int* ipiv, dcomplex* work, int* info);
  void zsytrf_(const char* uplo, const int* N, dcomplex* A, const int* lda, int* ipiv, dcomplex* work, int* lwork, int* info);
  void dsytri_(const char* uplo, const int* N, double* A, const int* lda, int* ipiv, double* work, int* info);
  void dsytrf_(const char* uplo, const int* N, double* A, const int* lda, int* ipiv, double* work, int* lwork, int* info);
}

inline void xgemm(const char* transa, const char* transb, const int m, const int n,
		  const int k, const double alpha, const double* A,
		  const int lda, const double* B, const int ldb, const double beta,
		  double* C, const int ldc)
{
  dgemm_(transa, transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

inline void xsymm(const char* side, const char* uplo, int m, int n, double alpha, const double* A, int lda, const double* B, int ldb, double beta, double* C, int ldc)
{
  dsymm_(side, uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

inline void xgemm(const char* transa, const char* transb, const int m, const int n,
		  const int k, const dcomplex& alpha, const dcomplex* A,
		  const int lda, const dcomplex* B, const int ldb, const dcomplex& beta,
		  dcomplex* C, const int ldc)
{
  zgemm_(transa, transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

inline int xgetrf(int n, double* a, int lda, int* ipiv, int ldb)
{
  int info = 0;
  dgetrf_(&n, &n, a, &lda, ipiv, &info);
  if (info){
    std::cerr << "Something wrong in LU (real) decomposition! " << info << std::endl;
  }  
  return info;
}

inline int xgetrf(int n, dcomplex* a, int lda, int* ipiv, int ldb)
{
  int info = 0;
  zgetrf_(&n, &n, a, &lda, ipiv, &info);
  if (info){
    std::cerr << "Something wrong in LU (complex) decomposition! " << info << std::endl;
  }  
  return info;
}

inline int xgetrs(int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb)
{
  int info = 0;
  dgetrs_("T", &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  if (info){
    std::cerr << "Something wrong with the system of (real) equations! " << info << std::endl;
  }  
  return info;
}

inline int xgetrs(int n, int nrhs, dcomplex* a, int lda, int* ipiv, dcomplex* b, int ldb)
{
  int info = 0;
  zgetrs_("T", &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  if (info){
    std::cerr << "Something wrong with the system of (complex) equations! " << info << std::endl;
  }
  return info;
}

template <class T>
inline int xgesv(int n, int nrhs, T* a, int lda, int* ipiv, T* b, int ldb)
{
  int info = xgetrf(n, a, lda, ipiv, ldb);
  if (info) return info;
  return xgetrs(n, nrhs, a, lda, ipiv, b, ldb);
}

inline int xgeev(const int N, double* A, const int lda, dcomplex* w, double* wr, double* wi, double* work,
		 const int lwork, double*)
{
  int ena = 1, info = 0;
  dgeev_("N", "N",  &N,  A, &lda, wr, wi, NULL, &ena, NULL, &ena, work, &lwork, &info);
  if (info)
    std::cerr << "Can't compute eigenvalues! " << info << std::endl;
  for (int i=0; i<N; i++) { w[i].real()=wr[i]; w[i].imag()=wi[i];}
  return info;
}

// Eigenvalues of real symmetric matrix
inline int xsyev(const int N, double* A, const int lda, double* w, double* work, const int lwork)
{
  int info = 0;
  dsyev_("N", "U",  &N,  A, &lda, w, work, &lwork, &info);
  if (info)
    std::cerr << "Can't compute eigenvalues! " << info << std::endl;
  return info;
}

inline int xgeev(const int N, dcomplex* A, const int lda, dcomplex* w, double* wr, double* wi, dcomplex* work,
		 const int lwork, double* rwork)
{
  int ena = 1, info = 0;
  zgeev_("N", "N",  &N,  A, &lda, w, NULL, &ena, NULL, &ena, work, &lwork, rwork, &info);
  if (info)
    std::cerr << "Can't compute eigenvalues! " << info << std::endl;
  return info;
}

inline int xgesdd(bool vect, int N, dcomplex* A, int lda, double* S, dcomplex* U, int ldu, dcomplex* Vt, int ldvt,
		  dcomplex* work, int lwork, double* rwork, int *iwork)
{
  int info = 0;
  std::string job = vect ? "A" : "N";
#ifndef _HP_  
  zgesdd_(job.c_str(), &N, &N, A, &lda, S, U, &ldu, Vt, &ldvt, work, &lwork, rwork, iwork, &info);
#endif  
  if (info) {
    std::cerr << "Can't compute SVD of the kernel! " << info << std::endl;
  }
  return info;
}

inline int xgesdd(bool vect, int N, double* A, int lda, double* S, double* U, int ldu, double* Vt, int ldvt,
		  double* work, int lwork, double*, int *iwork)
{
  int info = 0;
  std::string job = vect ? "A" : "N";
#ifndef _HP_  
  dgesdd_(job.c_str(), &N, &N, A, &lda, S, U, &ldu, Vt, &ldvt, work, &lwork, iwork, &info);
#endif  
  if (info) {
    std::cerr << "Can't compute SVD of the kernel! " << info << std::endl;
  }
  return info;
}

inline void comparec(int* b, dcomplex* w, int* N)
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

inline int xgees(bool vect, int N, dcomplex* A, int lda, dcomplex* w, double*, double*, dcomplex* Z, int ldz, dcomplex* work,
		 int lwork, double* rwork, int* bwork)
{
  int sdim, info = 0;
  std::string job = vect ? "V" : "N";
  mzgees_(job.c_str(), "S", &comparec, &N, A, &lda, &sdim, w, Z, &ldz, work, &lwork, rwork, bwork, &info);
  if (info) std::cerr << "Can't perform Schur decomposition! " << info << std::endl;
  return info;
}

inline int xgees(bool vect, int N, double* A, int lda, dcomplex* w, double* wr, double* wi, double* Z, int ldz, double* work,
		 int lwork, double*, int* bwork)
{
  int sdim, info = 0;
  std::string job = vect ? "V" : "N";
  mdgees_(job.c_str(), "S", &compared, &N, A, &lda, &sdim, wr, wi, Z, &ldz, work, &lwork, bwork, &info);
  if (info) std::cerr << "Can't perform Schur decomposition! " << info << std::endl;
  for (int i=0; i<N; i++) {w[i].real()=wr[i]; w[i].imag()=wi[i];}
  return info;
}

template <enum TypeOfMatrix TN>
inline void xgemv(int m, int n, const dcomplex& alpha, const dcomplex* A, int lda,
		  const dcomplex* x, const dcomplex& beta, dcomplex* y)

{
  std::cerr << "Not implemented yet! " << std::endl;
}

template <>
inline void xgemv<_Normal>(int m, int n, const dcomplex& alpha, const dcomplex* A, int lda,
			  const dcomplex* x, const dcomplex& beta, dcomplex* y)
{
  int inc = 1;
  zgemv_("T", &m, &n, &alpha, A, &lda, x, &inc, &beta, y, &inc);
}

template <>
inline void xgemv<_Transpose>(int m, int n, const dcomplex& alpha, const dcomplex* A, int lda,
			     const dcomplex* x, const dcomplex& beta, dcomplex* y)
{
  int inc = 1;
  zgemv_("N", &m, &n, &alpha, A, &lda, x, &inc, &beta, y, &inc);
}

// template <bool Transpose>
// inline void xgemv(int m, int n, const dcomplex& alpha, const dcomplex* A, int lda,
// 		  const dcomplex* x, const dcomplex& beta, dcomplex* y)

// {
//   std::cerr << "Not implemented yet! " << std::endl;
// }

// template <>
// inline void xgemv<false>(int m, int n, const dcomplex& alpha, const dcomplex* A, int lda,
// 			  const dcomplex* x, const dcomplex& beta, dcomplex* y)
// {
//   int inc = 1;
//   zgemv_("T", &m, &n, &alpha, A, &lda, x, &inc, &beta, y, &inc);
// }

// template <>
// inline void xgemv<true>(int m, int n, const dcomplex& alpha, const dcomplex* A, int lda,
// 			     const dcomplex* x, const dcomplex& beta, dcomplex* y)
// {
//   int inc = 1;
//   zgemv_("N", &m, &n, &alpha, A, &lda, x, &inc, &beta, y, &inc);
// }

template <enum TypeOfMatrix TN>
inline void xgemv(int m, int n, double alpha, const double* A, int lda,
		  const double* x, const double beta, double* y)
{
  std::cerr << "Not implemented yet! " << std::endl;
}

template <>
inline void xgemv<_Normal>(int m, int n, double alpha, const double* A, int lda,
			  const double* x, const double beta, double* y)
{
  int inc = 1;
  dgemv_("T", &m, &n, &alpha, A, &lda, x, &inc, &beta, y, &inc);
}

template <>
inline void xgemv<_Transpose>(int m, int n, double alpha, const double* A, int lda,
			     const double* x, const double beta, double* y)
{
  int inc = 1;
  dgemv_("N", &m, &n, &alpha, A, &lda, x, &inc, &beta, y, &inc);
}

// template <bool Transpose>
// inline void xgemv(int m, int n, double alpha, const double* A, int lda,
// 		  const double* x, const double beta, double* y)
// {
//   std::cerr << "Not implemented yet! " << std::endl;
// }
  
// template <>
// inline void xgemv<false>(int m, int n, double alpha, const double* A, int lda,
// 			  const double* x, const double beta, double* y)
// {
//   int inc = 1;
//   dgemv_("T", &m, &n, &alpha, A, &lda, x, &inc, &beta, y, &inc);
// }

// template <>
// inline void xgemv<true>(int m, int n, double alpha, const double* A, int lda,
// 			     const double* x, const double beta, double* y)
// {
//   int inc = 1;
//   dgemv_("N", &m, &n, &alpha, A, &lda, x, &inc, &beta, y, &inc);
// }

inline double xnrm(int N, const dcomplex* x, int incx = 1)
{
  return dznrm2_(&N, x, &incx);
}

inline double xnrm(int N, const double* x, int incx = 1)
{
  return dnrm2_(&N, x, &incx);
}

inline void xdscal(int N, double alpha, dcomplex* zx, int incx = 1)
{
  zdscal_(&N, &alpha, zx, &incx);
}

inline void xdscal(int N, double alpha, double* zx, int incx = 1)
{
  dscal_(&N, &alpha, zx, &incx);
}

inline void xsytrf(const char* uplo, int N, dcomplex* A, int lda, int* ipiv, dcomplex* work, int lwork, int& info)
{
  zsytrf_(uplo, &N, A, &lda, ipiv, work, &lwork, &info);
}
inline void xsytrf(const char* uplo, int N, double* A, int lda, int* ipiv, double* work, int lwork, int& info)
{
  dsytrf_(uplo, &N, A, &lda, ipiv, work, &lwork, &info);
}
inline void xsytri(const char* uplo, int N, dcomplex* A, int lda, int* ipiv, dcomplex* work, int& info)
{
  zsytri_(uplo, &N, A, &lda, ipiv, work, &info);
}
inline void xsytri(const char* uplo, int N, double* A, int lda, int* ipiv, double* work, int& info)
{
  dsytri_(uplo, &N, A, &lda, ipiv, work, &info);
}

#endif

