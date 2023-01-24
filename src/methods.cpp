#include "methods.h"
#include "StdTypes.h"
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

  // TODO Check size of matrix, if needed resize
  // if ((matrix.num_rows() != B.size()) && (matrix.num_rows() != V.num_rows()) &&
  //     (matrix.num_rows() != V.num_cols())) {
  //   std::cerr << "Error: size of matrix (matrices) or vector are incorrect."
  //             << std::endl;
  //   exit(EXIT_FAILURE);
  // }

  if(!U.is_identity() || U.num_cols() != matrix.num_cols() || )

  // Housholder reduction to bidiagonal form.
  int m = matrix.num_rows();
  int n = matrix.num_cols();

  const double epsilon = std::numeric_limits<double>::epsilon();
  double norm = 0;

  MatrixD A = MatrixD(matrix);
  VectorD v(n);
  Double scale = 0.0;




  // Accumulation of right-hand transformations.
  // Accumulation of left-hand transformations.
  // Diagonalization of the bidiagonal form:
}
