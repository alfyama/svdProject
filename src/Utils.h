#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <fstream>
#include <sstream>
#include <filesystem>

#include "matrix.h"

#define SOLUTIONS_PATH "results"

template<class T, class=std::enable_if_t<std::is_arithmetic<T>::value>>
void readMatrixCsv(std::string fileName, std::vector<T>& mat, int &m, int &n){
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


template<class T, class=std::enable_if_t<std::is_arithmetic<T>::value>>
void writeSolutionToCsv(Vector<T> &sol, std::string filenName) {
  std::filesystem::create_directory(SOLUTIONS_PATH);
  std::string filePath = SOLUTIONS_PATH "/" + filenName;
  std::ofstream file(filePath);

  if (!file.is_open()) {
    std::cerr << "Failed to open file " << filePath << "\n";
    exit(EXIT_FAILURE);
  }

  int i;
  for (i = 0; i < sol.size(); i++) {
    file << sol[i];
    if (i != sol.size() - 1)
      file << ",";
  }

  file << "\n";
  file.close();
}

template<class T>
std::string createResultFileName(std::string filename, std::string type,
                                 std::string method) {

  int start = filename.find_last_of("/") + 1;
  int end = filename.find("matrix.csv");
  std::string sol = filename.substr(start, end - start);
  if (type == "float")
    sol += "F";
  else if (type == "double")
    sol += "D";
  else if (type == "long double")
    sol += "LD";
  else
    ;

  sol += method;
  sol += ".csv";
  return sol;
}

#endif // UTILS_H_
