// --------------------------------------------------------------------
// traceInv_Eigen.cpp
// --------------------------------------------------------------------
#include "Eigen/Dense"
#include <iostream>

double traceOfInverse(const Eigen::MatrixXd& A)
{
  const int n = A.rows();
  // LU factorisation with partial pivoting (good for general dense matrices)
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(A);
  
  double traceInv = 0.0;
  
  for (int i = 0; i < n; ++i){
    // e_i : i-th canonical basis vector
    Eigen::VectorXd ei = Eigen::VectorXd::Unit(n, i);
    // Solve A x = e_i  ->  x = A^{-1} e_i
    Eigen::VectorXd xi = lu.solve(ei);
    traceInv += xi(i);             // diagonal element (A^{-1})_ii
  }
  return traceInv;
}

// --------------------------------------------------------------------
/*
int main(){
  // Example: random 7×7 matrix with a positive diagonal shift so it's invertible
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(7,7);
  A += 7.0 * Eigen::MatrixXd::Identity(7,7);
  
  std::cout << "A =\n" << A << "\n\n";
  
  double trInv = traceOfInverse(A);
  std::cout << "Tr(A^{-1}) from LU-solves : " << trInv << '\n';
  
  // Cross-check: form the full inverse once (for small n this is fine)
  double trInv_check = A.inverse().trace();
  std::cout << "Tr(A^{-1}) from full inverse : " << trInv_check << '\n';
}
*/
