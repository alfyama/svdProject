#ifndef SVD_MEHTODS
#define SVD_METHODS

#include "StdTypes.h"
#include "matrix.h"
#include "vector.h"
#include <limits>


void GolubReinsch_svd(MatrixF &A, VectorF &w);

void GolubReinsch_svd(MatrixD &A, VectorD &w);

void GolubReinsch_svd(MatrixLD &A, VectorLD &w);


#endif
