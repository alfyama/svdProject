#ifndef UTILS_H_
#define UTILS_H_

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "matrix.h"

#define SOLUTIONS_PATH "results"

template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
void readMatrixCsv(std::string fileName, std::vector<T> &mat, int &m, int &n) {
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

template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
void writeSolutionToCsv(Vector<T> &sol, std::string filenName, auto starttime,
                        auto endtime, int m, int n) {
  std::filesystem::create_directory(SOLUTIONS_PATH);
  std::string filePath = SOLUTIONS_PATH "/" + filenName;
  std::ofstream file(filePath);

  std::chrono::duration<double> elapsed_seconds = endtime - starttime;

  if (!file.is_open()) {
    std::cerr << "Failed to open file " << filePath << "\n";
    exit(EXIT_FAILURE);
  }

  file << "Time,Rows,Cols\n";

  file << elapsed_seconds.count() << ",";
  file << m << "," << n <<"\n";

  std::cout << std::endl;
  int i;
  for (i = 0; i < sol.size(); i++) {
    file << std::fixed << std::setprecision(16) <<sol[i];
    if (i != sol.size() - 1)
      file << ",";
  }
  file.close();
}

#endif // UTILS_H_
