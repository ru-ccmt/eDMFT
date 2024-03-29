// @Copyright 2007 Kristjan Haule
// 
#include "sfunction.h"

void Multiply(function2D<double>& C, const function2D<double>& A, const function2D<double>& B)
{
  int N = A.size_N();
  int K = A.size_Nd();
  int M = B.size_Nd();
  C.N = N; C.Nd = M;
  if (K==1 && N==1 && M==1){
      C(0,0) = A(0,0)*B(0,0);
      return;
  }
  if ( (N+M+K)/3. < 12){
    const double* __restrict__ a = A.MemPt();
    const double* __restrict__ b = B.MemPt();
    double* __restrict__ c = C.MemPt();
    int lda = A.lda();
    int ldb = B.lda();
    int ldc = C.lda();
    
    double acc00, acc01, acc10, acc11;

    for (int j=0; j < M; j += 2){
      for(int i = 0; i < N; i += 2 ){
	acc00 = acc01 = acc10 = acc11 = 0;
	if (j+1<M && i+1<N){
	  for (int k = 0; k < K; k++){
	    acc00 += a[(i+0)*lda+k] * b[k*ldb+j+0];
	    acc01 += a[(i+0)*lda+k] * b[k*ldb+j+1];
	    acc10 += a[(i+1)*lda+k] * b[k*ldb+j+0];
	    acc11 += a[(i+1)*lda+k] * b[k*ldb+j+1];
	  }
	  c[(i+0)*ldc+j+0] = acc00;
	  c[(i+0)*ldc+j+1] = acc01;
	  c[(i+1)*ldc+j+0] = acc10;
	  c[(i+1)*ldc+j+1] = acc11;
	}else if (j+1<M && i+1==N){
	  for (int k = 0; k < K; k++){
	    acc00 += a[(i+0)*lda+k] * b[k*ldb+j+0];
	    acc01 += a[(i+0)*lda+k] * b[k*ldb+j+1];
	  }
	  c[(i+0)*ldc+j+0] = acc00;
	  c[(i+0)*ldc+j+1] = acc01;
	}else if (i+1<N && j+1==M){
	  for (int k = 0; k < K; k++){
	    acc00 += a[(i+0)*lda+k] * b[k*ldb+j+0];
	    acc10 += a[(i+1)*lda+k] * b[k*ldb+j+0];
	  }
	  c[(i+0)*ldc+j+0] = acc00;
	  c[(i+1)*ldc+j+0] = acc10;
	}else if (i+1==N && j+1==M){
	  for (int k = 0; k < K; k++) acc00 += a[i*lda+k] * b[k*ldb+j];
	  c[i*ldc+j] = acc00;
	}else{
	  std::cout<<"It should not happen!"<<std::endl;
	}
      }
    }
    /*
    const double* a = A.MemPt();
    const double* b = B.MemPt();
    double* c = C.MemPt();
    int lda = A.lda();
    int ldb = B.lda();
    int ldc = C.lda();
    for (int i=0; i<N; i++){
      for (int j=0; j<M; j++){
  	double sum=0;
  	for (int k=0; k<K; k++) sum += a[i*lda+k]*b[k*ldb+j];
  	c[i*ldc+j] = sum;
      }
    }
    */
    /*
    for (int i=0; i<N; i++){
      for (int j=0; j<M; j++){
   	double sum=0;
   	for (int k=0; k<K; k++) sum += A(i,k)*B(k,j);
   	C(i,j) = sum;
      }
    }
    */
  }else{
    C.MProduct(A,B); // Call to BLAS
  }
}




void Multiply_N_T(function2D<double>& C, const function2D<double>& A, const function2D<double>& B)
{
  int N = A.size_N();
  int K = A.size_Nd();
  int M = B.size_N();
  if (K!=B.size_Nd()){cerr<<"Multiply_N_T failed. Wrong sizes!"<<endl;}
  C.N = N; C.Nd = M;
  if (K==1 && N==1 && M==1){
      C(0,0) = A(0,0)*B(0,0);
      return;
  }
  if ( (N+M+K)/3. < 12){
    const double* __restrict__ a = A.MemPt();
    const double* __restrict__ b = B.MemPt();
    double* __restrict__ c = C.MemPt();
    int lda = A.lda();
    int ldb = B.lda();
    int ldc = C.lda();
    
    double acc00, acc01, acc10, acc11;

    for (int j=0; j < M; j += 2){
      for(int i = 0; i < N; i += 2 ){
	acc00 = acc01 = acc10 = acc11 = 0;
	if (j+1<M && i+1<N){
	  for (int k = 0; k < K; k++){
	    acc00 += a[(i+0)*lda+k] * b[(j+0)*ldb+k];
	    acc01 += a[(i+0)*lda+k] * b[(j+1)*ldb+k];
	    acc10 += a[(i+1)*lda+k] * b[(j+0)*ldb+k];
	    acc11 += a[(i+1)*lda+k] * b[(j+1)*ldb+k];
	  }
	  c[(i+0)*ldc+j+0] = acc00;
	  c[(i+0)*ldc+j+1] = acc01;
	  c[(i+1)*ldc+j+0] = acc10;
	  c[(i+1)*ldc+j+1] = acc11;
	}else if (j+1<M && i+1==N){
	  for (int k = 0; k < K; k++){
	    acc00 += a[(i+0)*lda+k] * b[(j+0)*ldb+k];
	    acc01 += a[(i+0)*lda+k] * b[(j+1)*ldb+k];
	  }
	  c[(i+0)*ldc+j+0] = acc00;
	  c[(i+0)*ldc+j+1] = acc01;
	}else if (i+1<N && j+1==M){
	  for (int k = 0; k < K; k++){
	    acc00 += a[(i+0)*lda+k] * b[(j+0)*ldb+k];
	    acc10 += a[(i+1)*lda+k] * b[(j+0)*ldb+k];
	  }
	  c[(i+0)*ldc+j+0] = acc00;
	  c[(i+1)*ldc+j+0] = acc10;
	}else if (i+1==N && j+1==M){
	  for (int k = 0; k < K; k++) acc00 += a[i*lda+k] * b[j*ldb+k];
	  c[i*ldc+j] = acc00;
	}else{
	  std::cout<<"It should not happen!"<<std::endl;
	}
      }
    }
  }else{
    C.Product(A,B); // Call to BLAS
  }
}
