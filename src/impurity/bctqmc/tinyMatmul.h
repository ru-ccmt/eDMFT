// tiny_matmul.hpp — Blitz++ × Eigen micro-GEMM with optional transposition
// (sizes 1…10 unrolled; dynamic fallback for larger).
//   tiny::matmul(C, A, B [, transA, transB]);
// Copyright © 2025 — MIT licence.

#ifndef TINY_MATMUL_HPP
#define TINY_MATMUL_HPP

#include <blitz/array.h>
#include "Eigen/Dense"
#include <cassert>

namespace tiny {
  //namespace detail {

// ---------------------------------------------------------------------------
// Fixed‑size kernels  (fully unrolled) --------------------------------------
// ---------------------------------------------------------------------------
/*
template<int N>
inline void mult_fixed(const double* A, const double* B, double* C, bool transA=false, bool transB=false) noexcept
{
  using M  = Eigen::Map<const Eigen::Matrix<double,N,N,Eigen::RowMajor>>;
  using MC = Eigen::Map<Eigen::Matrix<double,N,N,Eigen::RowMajor>>;

  MC Cv(C);

  if(!transA && !transB)      Cv.noalias() = M(A)            * M(B);
  else if(transA && !transB)  Cv.noalias() = M(A).transpose() * M(B);
  else if(!transA && transB)  Cv.noalias() = M(A)            * M(B).transpose();
  else                        Cv.noalias() = M(A).transpose() * M(B).transpose();
}
*/
// ---------------------------------------------------------------------------
// Dynamic kernel for larger M,N,K -------------------------------------------
// ---------------------------------------------------------------------------
/*
inline void mult_tiny(const double* A, int m, int k, bool transA,
		      const double* B, int n,        bool transB,
		      double*       C               // must point to an m×n block
		      ) noexcept
{
  using MapA = Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>;
  using MapB = MapA;
  using MapC = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>;

  int aRows = transA ? k : m;
  int aCols = transA ? m : k;
  int bRows = transB ? n : k;
  int bCols = transB ? k : n;

  MapA MA(A, aRows, aCols);
  MapB MB(B, bRows, bCols);
  MapC MC(C, m, n);

  if (!transA && !transB)        MC.noalias() = MA             * MB;
  else if ( transA && !transB)   MC.noalias() = MA.transpose() * MB;
  else if (!transA &&  transB)   MC.noalias() = MA             * MB.transpose();
  else                           MC.noalias() = MA.transpose() * MB.transpose();
}
*/
namespace detail {
// ---------------------------------------------------------------------------
// Fixed kernel for small N,M,K ----------------------------------------------
// ---------------------------------------------------------------------------
/*
template<int M,int K,int N>
inline void mult_fixed(const double* A,       // pointer to an M×K block
		       const double* B,       // pointer to a K×N block
		       double*       C,       // pointer to an M×N block
		       bool transA = false,
		       bool transB = false) noexcept
{
  // Eigen Maps with fixed compile-time sizes:
  using MapA = Eigen::Map< const Eigen::Matrix<double, M, K, Eigen::RowMajor> >;
  using MapB = Eigen::Map< const Eigen::Matrix<double, K, N, Eigen::RowMajor> >;
  using MapC = Eigen::Map<       Eigen::Matrix<double, M, N, Eigen::RowMajor> >;

  MapA Amap(A);
  MapB Bmap(B);
  MapC Cmap(C);
  
  if (!transA && !transB)      Cmap.noalias() = Amap             * Bmap;
  else if ( transA && !transB) Cmap.noalias() = Amap.transpose() * Bmap;
  else if (!transA &&  transB) Cmap.noalias() = Amap             * Bmap.transpose();
  else                         Cmap.noalias() = Amap.transpose() * Bmap.transpose();
}
*/  
template<int M,int K,int N>
inline void mult_fixed(const double *A, int lda,
		       const double *B, int ldb,
		       double *C, int ldc,
		       bool transA = false,
		       bool transB = false) noexcept
{
  // A dynamic Stride type so we can pass runtime lda,ldb,ldc
  using StrideType = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;
  // Fixed‐size mat maps, but with dynamic strides:
  using MapA = Eigen::Map<const Eigen::Matrix<double, M, K, Eigen::RowMajor>,0,StrideType>;
  using MapB = Eigen::Map<const Eigen::Matrix<double, K, N, Eigen::RowMajor>,0,StrideType>;
  using MapC = Eigen::Map<      Eigen::Matrix<double, M, N, Eigen::RowMajor>,0,StrideType>;

  StrideType strideA( lda, 1 );
  StrideType strideB( ldb, 1 );
  StrideType strideC( ldc, 1 );

  MapA Amap(A, M, K, strideA);
  MapB Bmap(B, K, N, strideB);
  MapC Cmap(C, M, N, strideC);

  if (!transA && !transB)      Cmap.noalias() = Amap             * Bmap;
  else if ( transA && !transB) Cmap.noalias() = Amap.transpose() * Bmap;
  else if (!transA &&  transB) Cmap.noalias() = Amap             * Bmap.transpose();
  else                         Cmap.noalias() = Amap.transpose() * Bmap.transpose();
}

  
// A fully‐strided tiny‐block multiply
inline void mult_tiny(const double* A, int m, int k, bool transA, int lda,
		      const double* B, int n,        bool transB, int ldb,
		      double*       C,                            int ldc) noexcept
{
  // Tell Eigen about runtime strides:
  using StridedA = Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>,0,Eigen::OuterStride<>>;
  using StridedB = StridedA;
  using StridedC = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>,0,Eigen::OuterStride<>>;
  
  int aRows = transA ? k : m, aCols = transA ? m : k;
  int bRows = transB ? n : k, bCols = transB ? k : n;

  StridedA MA(A, aRows, aCols, Eigen::OuterStride<>(lda));
  StridedB MB(B, bRows, bCols, Eigen::OuterStride<>(ldb));
  StridedC MC(C, m,    n,      Eigen::OuterStride<>(ldc));

  if (!transA && !transB)        MC.noalias() = MA             * MB;
  else if ( transA && !transB)   MC.noalias() = MA.transpose() * MB;
  else if (!transA &&  transB)   MC.noalias() = MA             * MB.transpose();
  else                           MC.noalias() = MA.transpose() * MB.transpose();
}
inline void mult_mixed_real_complex(const double* A,               int m, int k, bool transA, int lda,
                                    const std::complex<double>* B, int n,        bool transB, int ldb,
                                    std::complex<double>* C,                                  int ldc) noexcept
{
  using CD = std::complex<double>;
  using StridedR = Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>,0,Eigen::OuterStride<>>;
  using StridedC = Eigen::Map<const Eigen::Matrix<CD,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>,0,Eigen::OuterStride<>>;
  using StridedO = Eigen::Map<Eigen::Matrix<CD,Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>,0, Eigen::OuterStride<>>;
  
  const int aRows = transA ? k : m; const int aCols = transA ? m : k;
  const int bRows = transB ? n : k; const int bCols = transB ? k : n;

  StridedR MA(A, aRows, aCols, Eigen::OuterStride<>(lda));
  StridedC MB(B, bRows, bCols, Eigen::OuterStride<>(ldb));
  StridedO MC(C, m,     n,     Eigen::OuterStride<>(ldc));

  if (!transA && !transB)        MC.noalias() = MA             * MB;
  else if ( transA && !transB)   MC.noalias() = MA.transpose() * MB;
  else if (!transA &&  transB)   MC.noalias() = MA             * MB.transpose();
  else                           MC.noalias() = MA.transpose() * MB.transpose();
}
inline void mult_mixed_complex_real(const std::complex<double>* A, int m, int k, bool transA, int lda,
                                    const double* B,               int n,        bool transB, int ldb,
                                    std::complex<double>* C,                                  int ldc) noexcept
{
  using CD = std::complex<double>;
  using StridedC = Eigen::Map<const Eigen::Matrix<CD,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>,0,Eigen::OuterStride<>>;
  using StridedR = Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>,0,Eigen::OuterStride<>>;
  using StridedO = Eigen::Map<Eigen::Matrix<CD,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>,0,Eigen::OuterStride<>>;
  
  const int aRows = transA ? k : m; const int aCols = transA ? m : k;
  const int bRows = transB ? n : k; const int bCols = transB ? k : n;

  StridedC MA(A, aRows, aCols, Eigen::OuterStride<>(lda));
  StridedR MB(B, bRows, bCols, Eigen::OuterStride<>(ldb));
  StridedO MC(C, m,     n,     Eigen::OuterStride<>(ldc));
  
  if (!transA && !transB)        MC.noalias() = MA             * MB;
  else if ( transA && !transB)   MC.noalias() = MA.transpose() * MB;          // NOTE: transpose (not adjoint)
  else if (!transA &&  transB)   MC.noalias() = MA             * MB.transpose();
  else                           MC.noalias() = MA.transpose() * MB.transpose();
}
// a small‐block, fully‐strided complex GEMM
inline void mult_tiny_complex(const std::complex<double>* A, int m, int k, bool transA, int lda,
			      const std::complex<double>* B, int n,        bool transB, int ldb,
			      std::complex<double>*       C,                            int ldc) noexcept
{
  // Tell Eigen about runtime strides:
  using StridedA = Eigen::Map<const Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>,0,Eigen::OuterStride<>>;
  using StridedB = StridedA;
  using StridedC = Eigen::Map<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>,0,Eigen::OuterStride<>>;
  
  int aRows = transA ? k : m, aCols = transA ? m : k;
  int bRows = transB ? n : k, bCols = transB ? k : n;

  StridedA MA(A, aRows, aCols, Eigen::OuterStride<>(lda));
  StridedB MB(B, bRows, bCols, Eigen::OuterStride<>(ldb));
  StridedC MC(C, m,    n,      Eigen::OuterStride<>(ldc));

  if (!transA && !transB)        MC.noalias() = MA             * MB;
  else if ( transA && !transB)   MC.noalias() = MA.transpose() * MB;
  else if (!transA &&  transB)   MC.noalias() = MA             * MB.transpose();
  else                           MC.noalias() = MA.transpose() * MB.transpose();
}

inline void gemm_dyn(const double *A, const double *B, double *C, int m, int n, int k, int lda, int ldb, int ldc, bool transA=false, bool transB=false)
{
  // Simple triple-loop, block k for LLC if m and n are smallish
  constexpr int KBLOCK = 32;  // tune on your CPU
  for (int j = 0; j < n; ++j) {
    for (int ii = 0; ii < m; ++ii)
      C[ii*ldc + j] = 0.0;

    for (int kk = 0; kk < k; kk += KBLOCK) {
      const int kend = std::min(k, kk + KBLOCK);
      for (int i = 0; i < m; ++i) {
	double s = 0.0;
	for (int p = kk; p < kend; ++p) {
	  const double a = transA ? A[p*lda + i] : A[i*lda + p]; // A(p,i) : A(i,p) 
	  const double b = transB ? B[j*ldb + p] : B[p*ldb + j]; // B(j,p) : B(p,j)
  	  s += a * b;
	}
	C[i*ldc + j] += s;
      }
    }
  }
}
} // namespace detail

  
inline void dot(blitz::Array<double,2>& C,
		const blitz::Array<double,2>& A, const blitz::Array<double,2>& B, bool transA=false, bool transB=false)
{
  // 1) Compute the “logical” dimensions
  int m  = transA ? A.extent(1) : A.extent(0);
  int kA = transA ? A.extent(0) : A.extent(1);
  int kB = transB ? B.extent(1) : B.extent(0);
  int n  = transB ? B.extent(0) : B.extent(1);
  assert(kA == kB && "ERROR in dot: inner dims must match k_A!=k_B!");
  int k = kA;
  // 2) Sanity‐check that C’s backing store is big enough
  if (C.extent(0) < m || C.extent(1) < n)
    std::cout<<"ERROR in dot C.shape="<<C.shape()<<" expecting at least ("<<m<<","<<n<<") size"<<std::endl;
  // 3) Re‐bind C’s view to the top‐left m×n block (no allocation)
  if (C.extent(0)!=m or C.extent(1)!=n) C.reference(C(Range(0,m-1),Range(0,n-1)));
  
  if (m==1 and k==1 and n==1){
    C(0,0) = A(0,0) * B(0,0);
    return;
  }
  // 4) Grab pointers 
  const double* Ap = A.dataFirst();
  const double* Bp = B.dataFirst();
  double*       Cp = C.dataFirst();
  // 5) Grab leading dims
  int lda = A.stride(0);//transA ? A.stride(1) : A.stride(0);
  int ldb = B.stride(0);//transB ? B.stride(1) : B.stride(0);
  int ldc = C.stride(0);  // C is never transposed
  // 6) Dispatch: tiny “Eigen” path vs. generic block‐GEMM
  detail::mult_tiny(Ap, m, k, transA, lda, Bp, n, transB, ldb, Cp, ldc);
}

inline void dot(Array<std::complex<double>,2>& C,
		const Array<std::complex<double>,2>& A, const Array<std::complex<double>,2>& B, bool transA = false, bool transB = false)
{// C <- C * C
  // 1) logical dims
  int m  = transA ? A.extent(1) : A.extent(0);
  int kA = transA ? A.extent(0) : A.extent(1);
  int kB = transB ? B.extent(1) : B.extent(0);
  int n  = transB ? B.extent(0) : B.extent(1);
  assert(kA == kB && "dot: inner dims must match");
  int k = kA;
  // 2) Sanity‐check that C’s backing store is big enough
  if (C.extent(0) < m || C.extent(1) < n)
    std::cout<<"ERROR in dot C.shape="<<C.shape()<<" expecting at least ("<<m<<","<<n<<") size"<<std::endl;
  // 3) re‐bind C to the top‐left m×n block
  if (C.extent(0)!=m or C.extent(1)!=n) C.reference( C(Range(0,m-1), Range(0,n-1)));
  
  if (m==1 and k==1 and n==1){
    C(0,0) = A(0,0) * B(0,0);
    return;
  }
  // 4) raw pointers
  const std::complex<double>* Ap = A.dataFirst();
  const std::complex<double>* Bp = B.dataFirst();
  std::complex<double>*       Cp = C.dataFirst();
  // 5) leading dims (always row-stride)
  int lda = A.stride(0);
  int ldb = B.stride(0);
  int ldc = C.stride(0);
  // 6) dispatch
  detail::mult_tiny_complex(Ap, m, k, transA, lda, Bp, n, transB, ldb, Cp, ldc);
}
inline void dot(blitz::Array<std::complex<double>,2>& C,
                const blitz::Array<std::complex<double>,2>& A, const blitz::Array<double,2>& B, bool transA = false, bool transB = false)
{// C <- C * R   (single Eigen GEMM, no real/imag splits)
  int m  = transA ? A.extent(1) : A.extent(0);
  int kA = transA ? A.extent(0) : A.extent(1);
  int kB = transB ? B.extent(1) : B.extent(0);
  int n  = transB ? B.extent(0) : B.extent(1);
  assert(kA == kB && "dot: inner dims must match");
  int k = kA;
  if (C.extent(0) < m || C.extent(1) < n) 
    std::cout << "ERROR in dot C.shape=" << C.shape() << " expecting at least (" << m << "," << n << ") size\n";
  // Bind C to the top-left m×n block
  if (C.extent(0)!=m or C.extent(1)!=n) C.reference( C(blitz::Range(0,m-1), blitz::Range(0,n-1)) );
  
  if (m==1 && k==1 && n==1) { C(0,0) = A(0,0) * B(0,0); return; }
  
  const std::complex<double>* Ap = A.dataFirst();
  const double*               Bp = B.dataFirst();
  std::complex<double>*       Cp = C.dataFirst();

  int lda = A.stride(0);   // Blitz row-major: leading dim = row stride
  int ldb = B.stride(0);
  int ldc = C.stride(0);
  // Assumes no aliasing between C and A/B.
  detail::mult_mixed_complex_real(Ap, m, k, transA, lda, Bp, n, transB, ldb, Cp, ldc);
}
inline void dot(blitz::Array<std::complex<double>,2>& C,
                const blitz::Array<double,2>& A, const blitz::Array<std::complex<double>,2>& B, bool transA = false, bool transB = false)
{// C <- R * C  (single Eigen GEMM)
  int m  = transA ? A.extent(1) : A.extent(0);
  int kA = transA ? A.extent(0) : A.extent(1);
  int kB = transB ? B.extent(1) : B.extent(0);
  int n  = transB ? B.extent(0) : B.extent(1);
  assert(kA == kB && "dot: inner dims must match");
  int k = kA;
  if (C.extent(0) < m || C.extent(1) < n)
    std::cout<<"ERROR in dot C.shape="<<C.shape()<<" expecting at least ("<<m<<","<<n<<") size\n";
  // bind C to m×n (top-left) and call one mixed kernel
  if (C.extent(0)!=m or C.extent(1)!=n) C.reference( C(blitz::Range(0,m-1), blitz::Range(0,n-1)) );
  
  if (m==1 && k==1 && n==1) { C(0,0) = A(0,0) * B(0,0); return; }

  const double*               Ap = A.dataFirst();
  const std::complex<double>* Bp = B.dataFirst();
  std::complex<double>*       Cp = C.dataFirst();

  int lda = A.stride(0);
  int ldb = B.stride(0);
  int ldc = C.stride(0);

  // Contract: C must not alias A or B. (Your pipeline already enforces this.)
  detail::mult_mixed_real_complex(Ap, m, k, transA, lda, Bp, n, transB, ldb, Cp, ldc);
}
/*
inline void dot(Array<std::complex<double>,2>& C,
		const Array<std::complex<double>,2>& A, const Array<double,2>& B, bool transA = false, bool transB = false)
{
  // 1) logical dims
  int m  = transA ? A.extent(1) : A.extent(0);
  int kA = transA ? A.extent(0) : A.extent(1);
  int kB = transB ? B.extent(1) : B.extent(0);
  int n  = transB ? B.extent(0) : B.extent(1);
  assert(kA == kB && "dot: inner dims must match");
  int k = kA;
  if (m==1 and k==1 and n==1){
    C(0,0) = A(0,0) * B(0,0);
    return;
  }
  // 2) Sanity‐check that C’s backing store is big enough
  if (C.extent(0) < m || C.extent(1) < n){
    std::cout<<"ERROR in dot C.shape="<<C.shape()<<" expecting at least ("<<m<<","<<n<<") size"<<std::endl;
  }

  // temporary storage for real matrices
  Array<double,2> Ar(A.shape()), Cr(m,n), Ci(m,n);
  Ar = real(A);
  // 4) raw pointers
  double* Ap1 = Ar.dataFirst();
  const double* Bp = B.dataFirst();
  double*       Cpr = Cr.dataFirst();
  int lda = Ar.stride(0);
  int ldb = B.stride(0);
  int ldc = Cr.stride(0);
  detail::mult_tiny(Ap1, m, k, transA, lda, Bp, n, transB, ldb, Cpr, ldc);

  //cout<<"Cr="<< Cr << endl;
  
  Ar = imag(A);
  Ap1 = Ar.dataFirst();
  double*       Cpi = Ci.dataFirst();
  detail::mult_tiny(Ap1, m, k, transA, lda, Bp, n, transB, ldb, Cpi, ldc);

  //cout<<"Ci="<< Ci << endl;
  
  // 3) re‐bind C to the top‐left m×n block
  C.reference( C(Range(0,m-1), Range(0,n-1)));
  for (int i=0; i<m; ++i)
    for (int j=0; j<n; ++j)
      C(i,j) = std::complex<double>(Cr(i,j),Ci(i,j));
}
*/
void minus_hermitian_conjugate(Array<std::complex<double>, 2>& A)
{// produces A -> -A^\dagger
  int N = A.extent(0);
  assert(A.extent(1) == N && "Matrix must be square for in-place Hermitian conjugation");
  using std::conj;
  for (int i = 0; i<N; ++i) {
    // Diagonal: just conjugate
    A(i, i) = -conj(A(i, i));
    for (int j=i+1; j<N; ++j) {
      std::complex<double> temp = conj(A(i, j));
      A(i, j) = -conj(A(j, i));
      A(j, i) = -temp;
    }
  }
}
// ---------------------------------------------------------------------------
// Public wrapper ------------------------------------------------------------
// ---------------------------------------------------------------------------
/*
inline void matmul(blitz::Array<double,2>& C,
		   const blitz::Array<double,2>& A, const blitz::Array<double,2>& B, bool transA = false, bool transB = false)
{
  // After transposition, matrices must be square and compatible.
  const int n = static_cast<int>(transA ? A.extent(1) : A.extent(0));
  assert((transA?A.extent(0):A.extent(1)) == n);
  assert((transB?B.extent(1):B.extent(0)) == n);
  assert((transB?B.extent(0):B.extent(1)) == n);
  if (C.extent(0)!=n || C.extent(1)!=n) C.resize(n,n);
  const double* Ap = A.dataFirst();
  const double* Bp = B.dataFirst();
  double*       Cp = C.dataFirst();
  switch(n) {
  case 1:  detail::mult_fixed<1 >(Ap,Bp,Cp,transA,transB); break;
  case 2:  detail::mult_fixed<2 >(Ap,Bp,Cp,transA,transB); break;
  case 3:  detail::mult_fixed<3 >(Ap,Bp,Cp,transA,transB); break;
  case 4:  detail::mult_fixed<4 >(Ap,Bp,Cp,transA,transB); break;
  case 5:  detail::mult_fixed<5 >(Ap,Bp,Cp,transA,transB); break;
  case 6:  detail::mult_fixed<6 >(Ap,Bp,Cp,transA,transB); break;
  case 7:  detail::mult_fixed<7 >(Ap,Bp,Cp,transA,transB); break;
  case 8:  detail::mult_fixed<8 >(Ap,Bp,Cp,transA,transB); break;
  case 9:  detail::mult_fixed<9 >(Ap,Bp,Cp,transA,transB); break;
  case 10: detail::mult_fixed<10>(Ap,Bp,Cp,transA,transB); break;
  default: detail::mult_dyn       (Ap,Bp,Cp,n,transA,transB);
  }
  
}
*/

/*
enum class Transpose { No, Yes };
inline void dot_T(blitz::Array<double,2>& C,
		  const blitz::Array<double,2>& A, Transpose ta,
		  const blitz::Array<double,2>& B, Transpose tb)
{
  const int Am = ta == Transpose::No ? A.extent(0) : A.extent(1);
  const int Ak = ta == Transpose::No ? A.extent(1) : A.extent(0);
  const int Bk = tb == Transpose::No ? B.extent(0) : B.extent(1);
  const int Bn = tb == Transpose::No ? B.extent(1) : B.extent(0);
  assert(Ak == Bk && "Inner dimensions must match");
  if (C.extent(0) != Am || C.extent(1) != Bn) C.resize(Am, Bn);
  const double* Ap = A.dataFirst();
  const double* Bp = B.dataFirst();
  double*       Cp = C.dataFirst();
  detail::gemm_dyn(Ap, Bp, Cp, Am, Bn, Ak, ta, tb);
}
*/
inline void dot_diag(blitz::Array<std::complex<double>,2>& C, const blitz::Array<double,1>& b)
{// This multiplies matrix with diagonal matrix stored as a vector
 // In Python the equivalent code would be
 //     C = C @ diag(b)
 //
  const int m = C.extent(0);
  const int n = C.extent(1);
  assert(b.extent(0) == n && "Length(b) must equal #cols(A)");
  
  const double*     __restrict bp = b.dataFirst();
  complex<double>*  __restrict Cp = C.dataFirst();
  const int lda = C.stride(0);;                // leading dimension
  for (int i = 0; i < m; ++i) {
    complex<double>*       __restrict Ci = Cp + i*lda;
    //#pragma omp simd                       // auto-vectorises the inner loop
    for (int j = 0; j < n; ++j) Ci[j] *= bp[j];
  }
}
inline void dot_diag(blitz::Array<double,2>& C, const blitz::Array<double,1>& b)
{// This multiplies matrix with diagonal matrix stored as a vector
 // In Python the equivalent code would be
 //     C = C @ diag(b)
 //
  const int m = C.extent(0);
  const int n = C.extent(1);
  assert(b.extent(0) == n && "Length(b) must equal #cols(A)");
  
  const double* __restrict bp = b.dataFirst();
  double*       __restrict Cp = C.dataFirst();
  const int lda = C.stride(0);;                // leading dimension
  for (int i = 0; i < m; ++i) {
    double*     __restrict Ci = Cp + i*lda;
    //#pragma omp simd                       // auto-vectorises the inner loop
    for (int j = 0; j < n; ++j) Ci[j] *= bp[j];
  }
}
inline void dot_diag(const blitz::Array<double,1>& b, blitz::Array<std::complex<double>,2>& C)
{//  This multiplies diagonal matrix with another matrix. 
 //  In Python the equivalent code would be 
 //       C = diag(b) @ C
 //
  const int m = C.extent(0);
  const int n = C.extent(1);
  if (b.extent(0)!=m){
    cout<<"b.N="<<b.extent(0)<<" m="<<m<<endl;
  }
  assert(b.extent(0) == m && "Length(b) must equal #rows(A)");
  const double*    __restrict bp = b.dataFirst();      // row scalars
  complex<double>* __restrict Cp = C.dataFirst();
  const int ldc = C.stride(0);
  for (int i = 0; i < m; ++i) {
    const double  scale = bp[i];          // scalar for this row
    //const double* __restrict Ai = Ap + i*lda;
    complex<double>*       __restrict Ci = Cp + i*ldc;
    /* contiguous inner loop → vectorises automatically            */
    /* ‘simd’ is a hint; leave it out if you don’t use OpenMP.     */
    //#pragma omp simd
    for (int j = 0; j < n; ++j) Ci[j] *= scale;
  }
}

inline void dot_diag(blitz::Array<double,2>& C,
                     const blitz::Array<double,2>& A, const blitz::Array<double,1>& b)
{// This multiplies matrix with diagonal matrix stored as a vector
 // In Python the equivalent code would be
 //     C = A @ diag(b)
 //
  const int m = A.extent(0);
  const int n = A.extent(1);
  assert(b.extent(0) == n && "Length(b) must equal #cols(A)");
  
  if (C.extent(0) != m || C.extent(1) != n) C.resize(m, n);
  
  const double* __restrict bp = b.dataFirst();
  const double* __restrict Ap = A.dataFirst();   // row-major
  double*       __restrict Cp = C.dataFirst();
  
  const int lda = A.stride(0);                // leading dimension
  const int ldc = C.stride(0);
  for (int i = 0; i < m; ++i) {
    const double* __restrict Ai = Ap + i*lda;
    double*       __restrict Ci = Cp + i*ldc;
    //#pragma omp simd                       // auto-vectorises the inner loop
    for (int j = 0; j < n; ++j) Ci[j] = Ai[j] * bp[j];
  }
}
inline void dot_diag(blitz::Array<double,2>& C,
		     const blitz::Array<double,1>& b, const blitz::Array<double,2>& A)
{//  This multiplies diagonal matrix with another matrix. 
 //  In Python the equivalent code would be 
 //       C = diag(b) @ A
 //
  const int m = A.extent(0);
  const int n = A.extent(1);
  assert(b.extent(0) == m && "Length(b) must equal #rows(A)");

  /* resize C if needed, but keep alias-safe: C may be A */
  if (C.extent(0) != m || C.extent(1) != n) C.resize(m, n);

  const double* __restrict bp = b.dataFirst();      // row scalars
  const double* __restrict Ap = A.dataFirst();      // Blitz is row-major
  double*       __restrict Cp = C.dataFirst();

  const int lda = A.stride(0);        // leading dimension (= #cols)
  const int ldc = C.stride(0);
  for (int i = 0; i < m; ++i) {
    const double  scale = bp[i];          // scalar for this row
    const double* __restrict Ai = Ap + i*lda;
    double*       __restrict Ci = Cp + i*ldc;
    /* contiguous inner loop → vectorises automatically            */
    /* ‘simd’ is a hint; leave it out if you don’t use OpenMP.     */
    //#pragma omp simd
    for (int j = 0; j < n; ++j) Ci[j] = scale * Ai[j];
  }
}

template<typename T>
inline T trace_dot(const blitz::Array<T,2>& A, const blitz::Array<T,2>& B)
{// This function computes the trace of the product of two matrices A and B, i.e.,
  // Tr(A*B)
  const int m = A.extent(0);
  const int k = A.extent(1);
  // sanity checks: AB must be square m×m
  assert(B.extent(0) == k && B.extent(1) == m);
  const T* __restrict Ap = A.dataFirst();      // row-major
  const T* __restrict Bp = B.dataFirst();
  const int lda = A.stride(0);                                // leading dims
  const int ldb = B.stride(0);
  T sum = 0.0;
  for (int i = 0; i < m; ++i) {
    const T* __restrict Ai = Ap + i*lda;     // row -i of A
    /*  Inner loop walks the *columns* of A (contiguous)
	and the *rows* of B (strided by m), which is fine for k ≤ 10.    */
    //#pragma omp simd reduction(+:sum)
    for (int j = 0; j < k; ++j)
      sum += Ai[j] * Bp[j*ldb + i];             // A_{ij} * B_{ji}
  }
  return sum;
}
inline complex<double> trace_dot(const blitz::Array<double,2>& A, const blitz::Array<complex<double>,2>& B)
{// This function computes the trace of the product of two matrices A and B, i.e.,
  // Tr(A*B)
  const int m = A.extent(0);
  const int k = A.extent(1);
  // sanity checks: AB must be square m×m
  assert(B.extent(0) == k && B.extent(1) == m);
  const double* __restrict Ap = A.dataFirst();      // row-major
  const complex<double>* __restrict Bp = B.dataFirst();
  const int lda = A.stride(0);                                // leading dims
  const int ldb = B.stride(0);
  complex<double> sum = 0.0;
  for (int i = 0; i < m; ++i) {
    const double* __restrict Ai = Ap + i*lda;     // row -i of A
    /*  Inner loop walks the *columns* of A (contiguous)
	and the *rows* of B (strided by m), which is fine for k ≤ 10.    */
    //#pragma omp simd reduction(+:sum)
    for (int j = 0; j < k; ++j)
      sum += Ai[j] * Bp[j*ldb + i];             // A_{ij} * B_{ji}
  }
  return sum;
}
inline complex<double> trace_dot(const blitz::Array<complex<double>,2>& A, const blitz::Array<double,2>& B)
{// This function computes the trace of the product of two matrices A and B, i.e.,
  // Tr(A*B)
  const int m = A.extent(0);
  const int k = A.extent(1);
  // sanity checks: AB must be square m×m
  assert(B.extent(0) == k && B.extent(1) == m);
  const complex<double>* __restrict Ap = A.dataFirst();      // row-major
  const double* __restrict Bp = B.dataFirst();
  const int lda = A.stride(0);                                // leading dims
  const int ldb = B.stride(0);
  complex<double> sum = 0.0;
  for (int i = 0; i < m; ++i) {
    const complex<double>* __restrict Ai = Ap + i*lda;     // row -i of A
    /*  Inner loop walks the *columns* of A (contiguous)
	and the *rows* of B (strided by m), which is fine for k ≤ 10.    */
    //#pragma omp simd reduction(+:sum)
    for (int j = 0; j < k; ++j)
      sum += Ai[j] * Bp[j*ldb + i];             // A_{ij} * B_{ji}
  }
  return sum;
}
complex<double> trace(const Array<complex<double>,2>& A)
{
  assert(A.extent(0) == A.extent(1) && "tiny::trace: ERROR Matrix must be squared.");
  complex<double> csum=0;
  for (int i=0; i<A.extent(0); ++i) csum += A(i,i);
  return csum;
}

inline void cdiag(Array<complex<double>,2>& GGi, const Array<double,1>& Fi)
{
  const int n = GGi.extent(0);
  const int m = GGi.extent(1);
  const int k = Fi.extent(0);
  assert(n == m && n == k && "tiny::diag dimensions of GG and Fi dont match.");
  GGi = 0;
  for (int i=0; i<n; ++i) GGi(i,i) = Fi(i);
}
} // namespace tiny

#endif // TINY_MATMUL_HPP


