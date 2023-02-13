#ifndef JACOBIROT_H_
#define JACOBIROT_H_

#include "matrix.h"
#include "vector.h"

template <class T> void jacobiRotation(Matrix<T> &B, Vector<T> &w) {

  int n, m;

  m = B.num_rows();
  n = B.num_cols();

  int vs = 1;
  T epsilon = 1e-12;
  T Norm;
  T alpha;
  Matrix<T> u(m, vs);
  Matrix<T> v(m, vs);
  Matrix<T> H(m, m);
  Matrix<T> I(m, m);

  I.Identity();
  H.Identity();
  for (int i = 0; i < n; i++) {
    u.zero();
    v.zero();

    Norm = 0.0;
    for (int j = i; j < m; j++) {
      u(j, 0) = B(j, i);         // B.data_[j * num_cols_ + i];
      Norm += u(j, 0) * u(j, 0); // u.data_[j] * u.data_[j];
    }

    Norm = sqrt(Norm);

    alpha = u(i, 0) < 0 ? Norm : -Norm;

    Norm = 0.0;
    for (int j = i; j < m; j++) {
      v(j, 0) = j == i ? u(j, 0) + alpha : u(j, 0); // u.data_[j] + alpha : u.data_[j];
      Norm += v(j, 0) * v(j, 0);  // v.data_[j] * v.data_[j];
    }
    Norm = sqrt(Norm);

    if (Norm < epsilon)
      continue;

    for (int j = i; j < m; j++)
      v(j, 0) /= Norm;

    H = I - (v * v.transpose()) * 2.0;

    B = H * B;
  }

  T N;
  T sin;
  T cos;
  T r;
  T B_[4];
  T R_[4];

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      N += B(i, j) * B(i, j); // B.data_[i * n + j] * B.data_[i * n + j];
  N = sqrt(N);

  T s = 0.0;
  while (sqrt(s) <= (epsilon * epsilon) * (N * N)) {
    for (int i = 0; i < n - 1; i++) {
      for (int j = i + 1; j < n; j++) {
        s += B(i, j) * B(j, i);

        r = sqrt((B(i, i) * B(i, i)) + (B(j, i) * B(j, i)));

        cos = B(i, i) / r;
        sin = B(j, i) / r;
        // X = B.data_[i * num_cols_ + i]
        // Y = B.data_[j * num_cols_ + i]
        B_[0] = cos * B(i, i) + sin * B(j, i);
        // cos * B.data_[i * num_cols_ + i] + sin * B.data_[j * num_cols_ + i];
        B_[1] = cos * B(i, j) + sin * B(j, j);
        // cos * B.data_[i * num_cols_ + j] + sin * B.data_[j * num_cols_ + j];
        B_[2] = -sin * B(i, i) + cos * B(j, i);
        //-sin * B.data_[i * num_cols_ + i] + cos * B.data_[j * num_cols_ + i];
        B_[3] = -sin * B(i, j) + cos * B(j, i);
        //-sin * B.data_[i * num_cols_ + j] +
        // cos * B.data_[j * num_cols_ + j];
        // X = B.data_[i * num_cols_ + i]
        // Y = B.data_[j * num_cols_ + i]

        cos = B_[0] / sqrt((B_[0] * B_[0]) + (B_[1] * B_[1]));
        sin = B_[1] / sqrt((B_[0] * B_[0]) + (B_[1] * B_[1]));

        R_[0] = cos * B_[0] + sin * B_[1];
        R_[1] = -sin * B_[0] + cos * B_[1];
        R_[2] = cos * B_[2] + sin * B_[3];
        R_[3] = -sin * B_[2] + cos * B_[3];

        B(i, i) = R_[0];
        B(j, j) = R_[3];

        // U = U * R;
        // V = V * R;
      }
    }
  }
  // for(int i = 0; i < n; i++){
  //     w[i] = B(i,i);
  // }
  // w.reorder();
  // return singular values
}

#endif // JACOBIROT_H_
