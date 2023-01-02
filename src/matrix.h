#ifndef SVD_MATRIX
#define SVD_MATRIX

#include <cassert>
#include <cstddef>

template <class Mtype> class Matrix {
public:
  // Constructors
  Matrix();
  Matrix(int n, int m) : num_rows_(n), num_cols_(m), data_(new Mtype[n * m]) {}

  // Destructors
  ~Matrix() { delete[] data_; }

  // Dimensions methods
  size_t num_rows() const { return num_rows_; }

  size_t num_cols() const { return num_cols_; }

  // get/set element by index
  Mtype& operator()(int i, int j) {
      checkRange(i, j);
      return data_[i * num_cols_ + j];
  }
  const Mtype& operator()(int i, int j) const {
      checkRange(i, j);
      return data_[i*num_cols_ + j];
  }

  // TODO Row/Col operations

private:
  int num_rows_;
  int num_cols_;
  Mtype *data_;

  void checkRange(int i, int j) const {
      assert(0 <= i && i < num_rows_ &&
             0 <= j && j < num_cols_);
  }

};

#endif
