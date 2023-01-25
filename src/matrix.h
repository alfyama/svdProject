#ifndef SVD_MATRIX
#define SVD_MATRIX

#include <cassert>
#include <cstddef>
#include <iostream>
#include <omp.h>
#include <string>

#include "StdTypes.h"
#include "vector.h"

template <class Mtype> class Matrix {
public:
  // Constructors
  Matrix();
  explicit Matrix(int m, int n)
      : num_rows_(m), num_cols_(n), data_(new Mtype[n * m]) {}

  // Copy constructor
  Matrix(const Matrix<Mtype> &other) {
    num_cols_ = other.num_cols_;
    num_rows_ = other.num_rows_;
  }

    Matrix(int m, int n, const Mtype *data) : num_rows_(m), num_cols_(n){
      data_ = new Mtype[m * n];
      int i;
      for(i=0; i < m*n; i++){
        data_[i] = data[i];
      }
    }

  // Destructors
  ~Matrix() { delete[] data_; }

  // Dimensions methods
  inline int num_rows() const { return num_rows_; }

  inline int num_cols() const { return num_cols_; }

  // NOTE This is only for debugging
  void display() {
    int i, j;
    for (i = 0; i < num_rows_; i++) {
      for (j = 0; j < num_cols_; j++) {
        std::cout << data_[i * num_cols_ + j] << ' ';
      }
      std::cout << "\n";
    }
  }

  // get/set element by index
  Mtype &operator()(int i, int j) {
    checkRange(i, j);
    return data_[i * num_cols_ + j];
  }

  const Mtype &operator()(int i, int j) const {
    checkRange(i, j);
    return data_[i * num_cols_ + j];
  }

  // Row swap function
  // swaps row at index row1 and row2
  void row_swap(int row1, int row2) {
    checkRowRange(row1);
    checkRowRange(row2);
    std::swap(data_[row1 * num_cols_], data_[row2 * num_cols_], num_cols_);
  }

  void col_swap(int col1, int col2) {
    checkColumnRange(col1);
    checkColumnRange(col2);
    std::swap(data_[col1 * num_rows_], data_[col2 * num_rows_], num_rows_);
  }

  Matrix transpose() const {
    Matrix trp(num_cols_, num_rows_);
    int i, j;
    for (i = 0; i < num_rows_; i++) {
      for (j = 0; j < num_cols_; j++) {
        trp(j, i) = (*this)(i, j);
      }
    }
  }

  bool resize(int n, int m) {
    num_rows_ = n;
    num_cols_ = m;
    int n_elements = n * m;
    delete[] data_;
    data_ = new Mtype[n_elements];
    if (data_ != nullptr) {
      int i;
      for (i = 0; i < n_elements; i++)
        data_[i] = 0.0;
      return true;
    }
    return false;
  }

  bool is_identity() {
    int i, j;
    // TODO Do parrallel code to check if matrix is identity

    for (i = 0; i < num_rows_; i++) {
      for (j = 0; j < num_cols_; j++) {
        if (i != j && data_[i * num_cols_ + j] != 0)
          return false;
        if (i == j && data_[i * num_cols_ + j] != 1)
          return false;
      }
    }
    return true;
  }

  // Ovelaod operators
  template<class T> friend Matrix<T> operator+ (const Matrix<T>& lhs, const Matrix<T>& rhs);
  template<class T> friend Matrix<T> operator+ (const T& lhs, const Matrix<T> rhs);
  template<class T> friend Matrix<T> operator+ (const Matrix<T> lhs, const T& rhs);

  template<class T> friend Matrix<T> operator- (const Matrix<T>& lhs, const Matrix<T>& rhs);
  template<class T> friend Matrix<T> operator- (const T& lhs, const Matrix<T> rhs);
  template<class T> friend Matrix<T> operator- (const Matrix<T> lhs, const T& rhs);

  template<class T> friend Matrix<T> operator* (const Matrix<T>& lhs, const Matrix<T>& rhs);
  template<class T> friend Matrix<T> operator* (const T& lhs, const Matrix<T> rhs);
  template<class T> friend Matrix<T> operator* (const Matrix<T> lhs, const T& rhs);
  template<class T> friend Matrix<T> operator* (const Vector<T>& lhs, const Matrix<T>& rhs);
  template<class T> friend Vector<T> operator* (const Matrix<T>& lhs, const Vector<T>& rhs);


private:
  int num_rows_;
  int num_cols_;
  Mtype *data_;

  inline void checkRange(int i, int j) const {
    assert(0 <= i && i < num_rows_ && 0 <= j && j < num_cols_);
  }

  inline void checkColumnRange(int j) const {
    assert(0 <= j && j <= num_cols_);
  }

  inline void checkRowRange(int i) const { assert(0 <= i && i <= num_cols_); }
};


template <class Mtype> class Identity : public Matrix<Mtype> {
public:
  explicit Identity(int n) : Matrix<Mtype>::Matrix(n, n) {
    int i, j;
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        (*this)(i, j) = (i == j) ? 1 : 0;
      }
    }
  }

  explicit Identity(int m, int n) : Matrix<Mtype>::Matrix(n, n) {
    int i, j;
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
        (*this)(i, j) = (i == j) ? 1 : 0;
      }
    }
  }

    // This allows us to have rectangular matrices, but with a diagonal of ones
  Identity(const Matrix<Mtype> &M)
      : Matrix<Mtype>::Matrix(M.num_rows(), M.num_cols()) {
    (*this).resize(M.num_rows(), M.num_rows());
    int i, j;
    for (i = 0; i < M.num_rows(); i++) {
      for (j = 0; j < M.num_cols(); j++) {
        (*this)(i, j) = (i == j) ? 1 : 0;
      }
    }
  }
};

// TODO
/* Overload * operator */
// template<class T>
// Vector<T> operator* (const Vector<T>& lhs, const Matrix<T>& rhs){

//   Matrix<T> return_mat(rhs.num_rows_m, rhs.num_cols_);

// }


// template<class T>
// Vector<T> operator* (const Matrix<T>& lhs, const Vector<T>& rhs){

// }

template<class T>
Matrix<T> operator* (const Matrix<T>& lhs, const Matrix<T>& rhs){
  int rhs_rows = rhs.num_rows_;
  int rhs_cols = rhs.num_cols_;
  int lhs_rows = lhs.num_rows_;
  int lhs_cols = lhs.num_cols_;

  if (lhs_cols == rhs_rows){
    T *tempdata = new T[lhs_rows * rhs_cols];

    int lr, rc;
    // loop over  each row of LHS
    for(lr = 0; lr < lhs_rows; lr++){
      // loop over each column of RHS
      for(rc=0; rc < rhs_cols; rc++){
        T element = 0.0;

        int i;
        for(i = 0; i < lhs_cols; i++){
          element += lhs(lr, i) * rhs(i, rc);
        }
        tempdata[(lr * rhs_cols) + rc] = element;
      }
    }
    Matrix<T> result(lhs_rows, rhs_cols, tempdata);
    delete [] tempdata;
    return result;
  }
  else {
    Matrix<T> result(1,1);
    return result;
  }
}





// typedefs so that the user does not create
// matrices with types like string or char
// limits the class of matrices to numeric types
// without the need of creating an implementation for
// each type

typedef Matrix<Double> MatrixD;
typedef const Matrix<Double> cMatrixD;

typedef Identity<Double> IdentityD;
typedef const Identity<Double> cIdentityD;

typedef Matrix<LDouble> MatrixLD;
typedef const Matrix<LDouble> cMatrixLD;

#endif
