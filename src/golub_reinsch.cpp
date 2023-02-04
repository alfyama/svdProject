#include "golub_reinsch.h"
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

  double scale;
  double f, g, s, h, c;
  int k, l, j;

  g = 0.0;
  scale = 0.0;

  hv[i] = 0.0;

  // Scale
  for (l = i; l < m; l++) {
    scale += A(l, i) * A(l, i);
  }

  // We need to scale the values
  // to get the unitary vector that is needed
  // for the Householder transform
  if (scale != 0.0) {
    // for (l = i; l < m; l++) {
    //   A(l, i) /= scale;
    //   norm_col += A(l, i) * A(l, i);
    // }
    f = A(i, i);

    g = -copysign(sqrt(scale), f);
    // Note that we norm_col2 is not actually
    // the norm of the column that we find on the pseudo code
    // this is done so that we dont calculate sqrt that might introduce
    // rounding errors?

    h = f * g - scale;

    // Now the entry on the main diagonal
    // holds the norm of the first column of A
    A(i, i) = f - g;

    // Update the matrix
    // This performs the matix * matrix
    // of the Houselhoder Transform and the
    // original matrix A
    for (j = i + 1; j < n; j++) {
      s = 0.0;
      for (k = i; k < m; k++) {
        s += A(k, i) * A(k, j);
      }
      f = s / h;
      for (k = i; k < m; k++) {
        A(k, j) += f * A(k, i);
      }
    }
    // Undo scaling, this thogether with
    // for (l = i; l < m; l++) {
    //   A(l, i) *= scale;
    // }
  }
  w[i] = g;
}

void HouseholderTransformRow(MatrixD &A, int i, VectorD &hv) {
  int n, m;
  n = A.num_cols();
  m = A.num_rows();
  double sgn, norm_col, scale;
  double f, p, s, h;
  int k, l, j;

  norm_col = 0.0;

  for (l = i + 1; l < n; l++) {
    scale += A(i, l) * A(i, l);
  }

  if (scale != 0.0) {
    // Scale the rows
    // for (l = i + 1; l < n; l++) {
    //   A(i, l) /= scale;
    //   s += A(i, l) * A(i, l);
    // }
    f = A(i, i + 1);

    p = -copysign(sqrt(norm_col), f);
    h = f * p - s;
    A(i, i + 1) = f - p;

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
    // for (k = i + 1; k < n; k++) {
    //   A(i, k) *= scale;
    // }
  }
}

void HouseholderReductionToBidiagonal(MatrixD &A, VectorD &w_, VectorD &hv_,
                                      Double y_) {

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

  // Housholder reduction to bidiagonal form.
  int m = A.num_rows();
  int n = A.num_cols();

  double y = 0.0;

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
      HouseholderTransformCol(A, i, w, hv_);

#if DEBUG
      std::cout << "Transform col" << std::endl;
      A.display();
      std::cout << std::endl;
#endif
    }
    // Transformations for the rows
    if (i < n) {
      HouseholderTransformRow(A, i, hv_);
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
  y = y > abs(hv[i] + w[i]) ? y : abs(hv[i] + w[i]);
  y_ = y;
}

void Householder_svdcmp(MatrixD &A, VectorD &w) {

  int i, j, k, l, n, m;
  double c, s, f, g, h, y, z, x;


  VectorD hv(n);

  m = A.num_rows();
  n = A.num_cols();

  /* Householder's reduction to bidiagonal */
  HouseholderReductionToBidiagonal(A, w, hv, y);

  /* Diagonalization of the bidiagonal form */
  double eps = std::numeric_limits<double>::epsilon();

  eps *= y;

  for (k = n - 1; k >= 0; k--) {
    bool flag = true;
    for (l = k; l >= 0; l--) {
      if (abs(hv[l]) <= eps || l == 0) {
        flag = false;
        break;
      }
      if (abs(w[l - 1]) <= eps)
        break;
    }

    // If the flag is set to true
    // we need to
    if (flag) {
      c = 0;
      s = 1;
      for (i = l; i < k + 1; i++) {
        f = s * hv[i];
        hv[i] = c * hv[i];

        // Jump to convergence test
        if (abs(f) <= eps)
          break;
        g = w[i];
        w[i] = sqrt(f * f + g * g);
        h = w[i];
        c = g / h;
        s = -f / h;
      }
    }

    // Convergence
    z = w[k];
    if (k == l) {
      // Ensure all singular values
      // are positive
      if (z < 0.0) {
        w[i] = -z;
      }
      break;
    }

    /* Shift from bottom 2 x 2 minor */

    x = w[l];
    y = w[k - 1];
    g = hv[k - 1];
    h = hv[k];
    f = ((y - z) * (y + z) + (g - h)) / (2 * 0 * h * y);
    g = sqrt(f * f + 1.0);

    f = ((x - z) * (x + z) + (h * ((y / f + copysign(g, f))) - h)) / x;

    /* QR transformation */
    for (j = l; j <= k - 1; j++) {
      i = j + 1;
      g = hv[i];
      y = w[i];
      h = s * g;
      g = c * g;
      z = sqrt(f * f + h * h);
      hv[j] = z;
      c = f / z;
      s = h / z;
      f = x * c + g * s;
      g = g * c - x * s;
      h = y * s;
      y *= c;
      z = sqrt(f * f + h * h);

      // Rotation ?
      w[j] = z;
      if (z) {
        c = f * (1.0 / z);
        s = h * (1.0 / z);
      }
      f = c * g + s * y;
      x = c * y - s * g;
    }
  }

  hv[l] = 0.0;
  hv[k] = f;
  w[k] = x;
}
