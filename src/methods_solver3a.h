#ifndef SVD_MEHTODS
#define SVD_METHODS

#include "StdTypes.h"
#include "matrix.h"
#include "vector.h"
#include <limits>


//methods for solver3a
void thresholdcheck(VectorD supdiagonalvector);
int reversecountzeros(VectorD countingvector, int start, int end);
int reversecountnonzeros(VectorD countingvector, int start, int end);
VectorD calculatesingularvalues(VectorD diagonalvector, int B_size);
void segmentandprocess(VectorD diagvector, VectorD supdiagvector, int B_size, int B2_size, int B3_size);
VectorD Solver3_methods(VectorD diagvector, VectorD supdiagvector, int B_size);


#endif