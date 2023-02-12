#ifndef GR_METHOD
#define GR_METHOD

#include "StdTypes.h"
#include "matrix.h"
#include "vector.h"
#include <limits>

/* L2 norm that avoids underflow and overflow */
template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
T L2Norm(const T x, const T y) {

  auto xabs = abs(x);
  auto yabs = abs(y);

  auto xsqr = [](auto x) { return x == 0.0 ? 0.0 : (x) * (x); };

  return (xabs > yabs   ? xabs * sqrt(1.0 + xsqr(yabs / xabs))
          : yabs == 0.0 ? 0.0
                        : yabs * sqrt(1.0 + xsqr(xabs / yabs)));
}

template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
void HouseholderReductionToBidiagonal(Matrix<T> &A, Vector<T> &w,
                                       Vector<T> &hv, T &y_) {
  /* Householder method:
   * Input:
   * matrix: Original matrix A of size m x n
   * */

  // Housholder reduction to bidiagonal form.
  int m = A.num_rows();
  int n = A.num_cols();

  T y;
  T scale;
  T f, g, s, h;
  int k, l, j;

  y = 0.0;
  g = 0.0;

  int i;
  for (i = 0; i < n; ++i) {

    hv[i] = g;
    scale = 0.0;
    // Transformations for the cols
    if (i < m) {
      // HouseholderTransformColF(A, i, w_, hv_);
      // Scale
      for (l = i; l < m; l++) {
        scale += A(l, i) * A(l, i);
      }
      // We need to scale the values
      // to get the unitary vector that is needed
      // for the Householder transform
      if (scale == 0.0)
        g = 0.0;
      else {
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
      }
    }
    w[i] = g;

    g = 0.0;
    scale = 0.0;

    // Transformations for the rows
    if (i + 1 != n && i + 1 <= m) {
      // HouseholderTransformRowF(A, i, hv);
      for (l = i + 1; l < n; l++) {
        scale += A(i, l) * A(i, l);
      }

      if (scale == 0.0)
        g = 0.0;

      else {

        f = A(i, i + 1);

        g = -copysign(sqrt(scale), f);
        h = f * g - scale;
        A(i, i + 1) = f - g;

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
      }
    }
    y = y > abs(hv[i] + w[i]) ? y : abs(hv[i] + w[i]);
  }
  y_ = y;
}

template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
void GolubReinsch_svd(Matrix<T> &A, Vector<T> &w) {
  int i, j, k, l, n, iter;
  T c, s, f, g, h, y, z, x;

  y = 0.0;
  // m = A.num_rows();
  n = A.num_cols();

  Vector<T> hv(n);

  /* Householder's reduction to bidiagonal */
  HouseholderReductionToBidiagonal(A, w, hv, y);
#if DEBUG
  std::cout << "Householder Reduction To Bidiagonal" << std::endl;
  A.display();
  std::cout << std::endl;
#endif

  /* Diagonalization of the bidiagonal form */
  T eps = std::numeric_limits<T>::epsilon();

  eps *= y;

  for (k = n - 1; k >= 0; k--) {
    for (iter = 0; iter < 50; iter++) {
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
      // we do cancellation of hv[l]
      if (flag) {
        c = 0.0;
        s = 1.0;
        for (i = l; i < k + 1; i++) {
          f = s * hv[i];
          hv[i] = c * hv[i];

          // Jump to convergence test
          if (abs(f) <= eps)
            break;
          g = w[i];
          w[i] = L2Norm(f, g);
          // w[i] = sqrt(f * f + g * g);
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
        if (z < 0.0)
          w[i] = -z;

        break;
      }
      if (iter == 49) {
        std::cerr << "No convergence after 49 iterations" << std::endl;
        exit(EXIT_FAILURE);
      }

      /* Shift from bottom 2 x 2 minor */
      x = w[l];
      y = w[k - 1];
      g = hv[k - 1];
      h = hv[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = L2Norm(f, (T)1.0);
      // g = sqrt(f * f + 1.0);

      f = ((x - z) * (x + z) + h * ((y / (f + copysign(g, f))) - h)) / x;

      c = 1.0;
      s = 1.0;
      /* QR transformation */
      for (j = l; j <= k - 1; j++) {
        i = j + 1;
        g = hv[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = L2Norm(f, h);
        // z = sqrt(f * f + h * h);
        hv[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y *= c;
        z = L2Norm(f, h);
        // z = sqrt(f * f + h * h);

        // Rotation ?
        w[j] = z;
        if (z) {
          c = f * (1.0 / z);
          s = h * (1.0 / z);
        }
        f = c * g + s * y;
        x = c * y - s * g;
      }

      hv[l] = 0.0;
      hv[k] = f;
      w[k] = x;
    }
  }
}

#endif
