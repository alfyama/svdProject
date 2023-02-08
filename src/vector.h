#ifndef SVD_VECTOR
#define SVD_VECTOR

#include <cmath>

#include "StdTypes.h"

template <class Vtype> class Vector {
public:
  Vector();
  explicit Vector(int n) : num_elements_(n), data_(new Vtype[n]){};

  // Destructor
  ~Vector() { delete[] data_; }

  // Assigment operator
  Vector &operator=(const Vector &rhs);

  inline int size() const { return num_elements_; };

  // get/set
  inline Vtype &operator[](const int index) {
    checkRange(index);
    return data_[index];
  }
  inline const Vtype &operator[](const int index) const {
    checkRange(index);
    return data_[index];
  }

  Vtype *begin() { return data_; }
  Vtype *end() { return data_ + num_elements_; }

  bool resize() { return true; }

  // Vtype norm(){
  //   Vtype val;
  //   int i;
  //   for(i = 0; i < num_elements_; i++){
  //     val = data_[i] * data_[i];
  //   }
  //   return sqrt(val);
  // }

private:
  int num_elements_;
  Vtype *data_;

  inline void checkRange(int index) const {
    assert(0 <= index && index < num_elements_);
  }
};

typedef Vector<Float> VectorF;

typedef Vector<Double> VectorD;
typedef const Vector<Double> cVectorD;

typedef Vector<LDouble> VectorLD;
typedef const Vector<LDouble> cVectorLD;

#endif
