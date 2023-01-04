#ifndef SVD_MATRIX
#define SVD_MATRIX

#include <cassert>
#include <cstddef>
#include <string>

template <class Mtype> class Matrix {
public:
  // Constructors
  Matrix();
  Matrix(int n, int m) : num_rows_(n), num_cols_(m), data_(new Mtype[n * m]) {}

  // Destructors
  ~Matrix() { delete[] data_; }

  // Dimensions methods
  inline size_t num_rows() const { return num_rows_; }

  inline size_t num_cols() const { return num_cols_; }

  // get/set element by index
  Mtype& operator()(int i, int j) {
      checkRange(i, j);
      return data_[i * num_cols_ + j];
  }
  const Mtype& operator()(int i, int j) const {
      checkRange(i, j);
      return data_[i*num_cols_ + j];
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

  // TODO Resize matrix
  // Eliminates rows/columns with only zeros
  bool resize(){
      return true;
  }

private:
  int num_rows_;
  int num_cols_;
  Mtype *data_;

  inline void checkRange(int i, int j) const {
      assert(0 <= i && i < num_rows_ &&
             0 <= j && j < num_cols_);
  }

  inline void checkColumnRange(int j) const {
      assert(0 <= j && j <= num_cols_);
  }

  inline void checkRowRange(int i) const {
      assert(0 <= i && i <= num_cols_);
  }
};

#endif
