#ifndef SVD_MATRIX
#define SVD_MATRIX

#include <iostream>
#include <cassert>
#include <cstddef>
#include <string>
#include <omp.h>

#include "StdTypes.h"

template <class Mtype> class Matrix {
public:
  // Constructors
  Matrix();
  explicit Matrix(int n, int m)
      : num_rows_(n), num_cols_(m), data_(new Mtype[n * m]) {}

  // Copy constructor
  Matrix(const Matrix<Mtype> &other) {
    num_cols_ = other.num_cols_;
    num_rows_ = other.num_rows_;
  }

  // Destructors
  ~Matrix() { delete[] data_; }

  // Dimensions methods
  inline int num_rows() const { return num_rows_; }

  inline int num_cols() const { return num_cols_; }


  // NOTE This is only for debugging
    void display() {
      int i,j;
      for(i=0; i<num_rows_;i++){
        for(j=0; j<num_cols_;j++){
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
    int n_elements = n*m;
    delete [] data_;
    data_ = new Mtype[n_elements];
    if(data_ != nullptr){
      int i;
      for(i=0; i < n_elements; i++)
        data_[i] = 0.0;
      return true;
    }
    return false;
  }

  bool is_identity() {
    int i, j;
    // TODO Do parrallel code to check if matrix is identity

    for(i=0; i < num_rows_; i++){
      for(j=0; j< num_cols_; j++){
        if( i != j  && data_[i*num_cols_ + j] != 0) return false;
        if ( i == j && data_[i*num_cols_ + j] != 1) return false;
      }
    }
    return true;
  }



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


template<class Mtype>
class Identity : public Matrix<Mtype> {

public:
    explicit Identity(int n, int m) : Matrix<Mtype>::Matrix(n, m) {
      int i;
      for(i = 0; i < n; i++){
        (*this)(i,i) = 1;
      }
    }

};

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
