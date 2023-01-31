#include "methods.h"
#include "StdTypes.h"
#include "matrix.h"
#include <cstdlib>
#include <iostream>

/*
** NOTE the following Householder are based on the following notes
** https://www.cs.utexas.edu/users/flame/laff/alaff/chapter11-reduction-to-bidirectional-form.html
*/

void HouseholderTransformCol(MatrixD &A, int i, VectorD &w, VectorD &hv) {

  /*
  ** Computes the Householder vector and stores this vector
  ** on the zeros introduces by the Householder operation ?
  */

  int n, m;
  n = A.num_cols();
  m = A.num_rows();

  double sgn, norm_col, scale;
  double diagonal_entry, p, s, h;
  int k, l, j;

  norm_col = 0.0;

  // Scale
  for (l = i; l < m; l++) {
    scale += fabs(A(l, i));
  }

  // We need to scale the values
  // to get the unitary vector that is needed
  // for the Householder transform
  if (scale != 0.0) {
    for (l = i; l < m; l++) {
      A(l, i) /= scale;
      norm_col += A(l, i) * A(l, i);
    }
    diagonal_entry = A(i, i);
    sgn = diagonal_entry >= 0 ? 1 : -1;

    p = sqrt(norm_col) * sgn;
    // Note that we norm_col2 is not actually
    // the norm of the column that we find on the pseudo code
    // this is done so that we dont calculate sqrt that might introduce
    // rounding errors?

    h = diagonal_entry * p - norm_col;

    // Now the entry on the main diagonal
    // holds the norm of the first column of A
    A(i, i) = diagonal_entry - p;

    // Update the matrix
    // This performs the matix * matrix
    // of the Houselhoder Transform and the
    // original matrix A
    for (j = i + 1; j < n; j++) {
      s = 0.0;
      for (k = i; k < m; k++) {
        s += A(k, i) * A(k, j);
      }
      diagonal_entry = s / h;
      for (k = i; k < m; k++) {
        A(k, j) += diagonal_entry * A(k, i);
      }
    }
    // Undo scaling, this thogether with
    for (l = i; l < m; l++) {
      A(l, i) *= scale;
    }
  }
  w[i] = p * scale;

}

void HouseholderTransformRow(MatrixD &A, int i, VectorD &hv) {
  int n, m;
  n = A.num_cols();
  m = A.num_rows();
  double sgn, norm_col, scale;
  double uperdiag_entry, p, s, h;
  int k, l, j;

  norm_col = 0.0;

  for (l = i + 1; l < n; l++) {
    scale += fabs(A(i, l));
  }

  if (scale != 0.0) {
    // Scale the rows
    for (l = i + 1; l < n; l++) {
      A(i, l) /= scale;
      s += A(i, l) * A(i, l);
    }
    uperdiag_entry = A(i, i + 1);
    sgn = uperdiag_entry >= 0 ? 1 : -1;

    p = sqrt(norm_col) * sgn;
    h = uperdiag_entry * p - s;
    A(i, i + 1) = uperdiag_entry - p;

    for (j = i + 1; j < n; j++) {
      hv[j] = A(i, j) / h;
    }

    // Update the matrix
    for (l = i + 1; l < m; l++) {
      s = 0.0;
      for (k = i + 1; k < n; k++) {
        s += A(l, k) * A(i, k);
      }
      for (k = i + 1; k < n; k++) {
        A(l, k) += s * hv[k];
      }
    }
    for (k = i + 1; k < n; k++) {
      A(i, k) *= scale;
    }
  }
}

void HouseholderReductionToBidiagonal(MatrixD &A, VectorD &w_) {

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

  IdentityD U = IdentityD(A.num_rows());
  IdentityD V = IdentityD(A.num_cols());

  // Housholder reduction to bidiagonal form.
  int m = A.num_rows();
  int n = A.num_cols();

  const double epsilon = std::numeric_limits<double>::epsilon();
  double norm = 0;

  VectorD w(n);
  VectorD hv(n);

#if DEBUG
  std::cout << "Original Matrix" << std::endl;
  A.display();
  std::cout << std::endl;
#endif

  int i;
  for (i = 0; i < n; i++) {
    // Transformations for the cols
    if (i < m) {
      HouseholderTransformCol(A, i, w, hv);

#if DEBUG
  std::cout << "Transform col" << std::endl;
  A.display();
  std::cout << std::endl;
#endif

    }
    // Transformations for the rows
    if (i < n - 2) {
      HouseholderTransformRow(A, i, hv);
#if DEBUG
  std::cout << "Transform row" << std::endl;
  A.display();
  std::cout << std::endl;
#endif

    }
  }

#if DEBUG
  std::cout << "Bidiagonal matrix" << std::endl;
  A.display();
  std::cout << std::endl;
#endif
}

void Householder_svdcmp(MatrixD &A, VectorD &W) {

  int i, n, m;

  m = A.num_rows();
  n = A.num_cols();

  HouseholderReductionToBidiagonal(A, W);

  for (i = 0; i < n; i++) {
    int l = i;
  }
}
