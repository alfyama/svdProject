#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <chrono>

#include "Utils.h"
#include "golub_reinsch.h"
#include "matrix.h"
#include "dqds.h"
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

int main(int argc, char *argv[]) {

  std::string flagMethod;
  std::string flagType;
  std::string fileName;

  for (int i = 1; i < argc; ++i) {
    std::string argument(argv[i]);
    if (argument.substr(0, 8) == "-method=") {
      flagMethod = argument.substr(8);
      std::transform(flagMethod.begin(), flagMethod.end(), flagMethod.begin(),
                     ::tolower);
    } else if (argument.substr(0, 6) == "-type=") {
      flagType = argument.substr(6);
      std::transform(flagType.begin(), flagType.end(), flagType.begin(),
                     ::tolower);
    } else if (fileName.empty()) {
      fileName = argument;
    } else {
      std::cout << "Invalid argument: " << argument << std::endl;
    }
  }

  if (flagMethod.empty() || flagType.empty() || fileName.empty()) {
    std::cerr << "Usage: svd -method=<METHOD> -type=<TYPE> file" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Check type flag
  if (flagType == "float") {
    std::vector<float> data;
    int m, n;
    readMatrixCsv(fileName, data, m, n);
    MatrixF A(m, n, data.data());
    VectorF w(n);

#if DEBUG
    std::cout << "Original Matrix: " << std::endl;
    A.display();
    std::cout << std::endl;
#endif

    if (flagMethod == "gr") {
      std::cout << "Golub Reinsch method " << std::endl;
      // start method timer
      auto start = std::chrono::high_resolution_clock::now();
      GolubReinsch_svd(A, w);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Method finished " << std::endl;
      std::string solFileName =
          createResultFileName(fileName, flagType, flagMethod);
      writeSolutionToCsv(w, solFileName, start, end);

    } else if (flagMethod == "dqds") {
      std::cout << "Differential quotient-difference method " << std::endl;
      // start method timer
      auto start = std::chrono::high_resolution_clock::now();
      dqds_main(A, w);
      auto end = std::chrono::high_resolution_clock::now();
      std::string solFileName =
         createResultFileName(fileName, flagType, flagMethod);
      writeSolutionToCsv(w, solFileName, start, end);

    } else if (flagMethod == "gdg") {

    } else {
      std::cout << "Error -method=<METHOD>" << std::endl;
    }

  } else if (flagType == "double") {
    std::vector<double> data;
    int m, n;
    readMatrixCsv(fileName, data, m, n);
    MatrixD A(m, n, data.data());
    VectorD w(n);
    if (flagMethod == "gr") {
      std::cout << "Golub Reinsch method " << std::endl;
      // start method timer
      auto start = std::chrono::high_resolution_clock::now();
      GolubReinsch_svd(A, w);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Method finished " << std::endl;
      std::string solFileName =
          createResultFileName(fileName, flagType, flagMethod);
      writeSolutionToCsv(w, solFileName, start, end);

    } else if (flagMethod == "dqds") {
      std::cout << "Differential quotient-difference method " << std::endl;
      // start method timer
      auto start = std::chrono::high_resolution_clock::now();
      dqds_main(A, w);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Method finished " << std::endl;
      std::string solFileName =
          createResultFileName(fileName, flagType, flagMethod);
      writeSolutionToCsv(w, solFileName, start, end);

    } else if (flagMethod == "gdg") {

    } else {
      std::cout << "Error -method=<METHOD>" << std::endl;
    }

  } else if (flagType == "long double") {
    std::vector<long double> data;
    int m, n;
    readMatrixCsv(fileName, data, m, n);
    MatrixLD A(m, n, data.data());
    VectorLD w(n);
    if (flagMethod == "gr") {
      std::cout << "Golub Reinsch method " << std::endl;
      // start method timer
      auto start = std::chrono::high_resolution_clock::now();
      GolubReinsch_svd(A, w);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Method finished " << std::endl;
      std::string solFileName =
          createResultFileName(fileName, flagType, flagMethod);
      writeSolutionToCsv(w, solFileName, start, end);

    } else if (flagMethod == "dqds") {
      std::cout << "Differential quotient-difference method " << std::endl;
      // start method timer
      auto start = std::chrono::high_resolution_clock::now();
      dqds_main(A, w);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Method finished " << std::endl;
      std::string solFileName =
          createResultFileName(fileName, flagType, flagMethod);
      writeSolutionToCsv(w, solFileName, start, end);

    } else if (flagMethod == "gdg") {

    } else {
    std::cout << "Error -type=<TYPE>" << std::endl;
    exit(EXIT_FAILURE);
  }
  }

  exit(EXIT_SUCCESS);
}
