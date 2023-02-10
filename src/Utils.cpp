#include "Utils.h"

void readMatrixCsvF(std::string fileName, std::vector<float> &mat, int &m,
                    int &n) {
  std::ifstream file(fileName);
  std::string line;

  int rows = 0;
  int cols = 0;

  while (std::getline(file, line)) {
    std::istringstream iss(line);
    float num;
    int current_col = 0;

    rows++;
    while (iss >> num) {
      mat.push_back(num);
      if (iss.peek() == ',')
        iss.ignore();
      current_col++;
    }

    if (rows == 1) {
      cols = current_col;
    } else if (cols != current_col) {
      std::cerr << "Error reading matrix from input file. Inconsistent number "
                   "of columns"
                << std::endl;
    }
  }
  m = rows;
  n = cols;
}

void readMatrixCsvD(std::string fileName, std::vector<double> &mat, int &m,
                    int &n) {
  std::ifstream file(fileName);
  std::string line;

  int rows = 0;
  int cols = 0;

  while (std::getline(file, line)) {
    std::istringstream iss(line);
    double num;
    int current_col = 0;

    rows++;
    while (iss >> num) {
      mat.push_back(num);
      if (iss.peek() == ',')
        iss.ignore();
      current_col++;
    }
    if (rows == 1)
      cols = current_col;
  }
  m = rows;
  n = cols;
}

void readMatrixCsvLD(std::string fileName, std::vector<long double> &mat,
                     int &m, int &n) {
  std::ifstream file(fileName);
  std::string line;

  int rows = 0;
  int cols = 0;

  while (std::getline(file, line)) {
    std::istringstream iss(line);
    long double num;

    while (iss >> num) {
      mat.push_back(num);
      cols++;
    }
    rows++;
    if (rows == 1)
      cols /= rows;
  }
  m = rows;
  n = cols;
}
