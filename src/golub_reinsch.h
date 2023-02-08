#ifndef GR_METHOD
#define GR_METHOD

#include "StdTypes.h"
#include "matrix.h"
#include "vector.h"
#include <limits>


void GolubReinsch_svdF(MatrixF &A, VectorF &w);

void GolubReinsch_svdD(MatrixD &A, VectorD &w);

void GolubReinsch_svdLD(MatrixLD &A, VectorLD &w);


#endif
