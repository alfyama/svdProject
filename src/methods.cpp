#include "methods.h"
#include "StdTypes.h"
#include "matrix.h"
#include <cstdlib>
#include <iostream>

void Housholder_svdcmp(MatrixD &matrix, MatrixD &U ,MatrixD &B, MatrixD &V) {

  /* Householder method:
   * Input:
   * matrix: Original matrix of size m x n
   * U: Identity matrix of size m x n
   * B: Empty matrix, size m x n
   * V: Identity matrix, size n x n
   * Output: Matrices U, B and V s.t. matrix = U*B*V^T
   * U, V: products of Householder matrices
   * B: upper diagonal matrix
   * */

  U = IdentityD(U);
  V = IdentityD(matrix.num_cols());

  // Housholder reduction to bidiagonal form.
  int m = matrix.num_rows();
  int n = matrix.num_cols();

  const double epsilon = std::numeric_limits<double>::epsilon();
  double norm = 0;

  MatrixD A = MatrixD(matrix);
  VectorD v(n);
  Double scale = 0.0;





}
