#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <fstream>
#include <sstream>

#include "matrix.h"


void readMatrixCsv(std::string fileName, std::vector<float>& mat, int &m, int &n);

void readMatrixCsv(std::string fileName, std::vector<double>& mat, int &m, int &n);

void readMatrixCsv(std::string fileName, std::vector<long double>& mat, int &m, int &n);


#endif // UTILS_H_
