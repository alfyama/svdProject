#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <filesystem>


#include "Utils.h"
#include "golub_reinsch.h"
#include "matrix.h"
#include "methods_solver3a.h"

#define SOLUTIONS_PATH "results"

void writeSolutionToCsv(VectorF &sol, std::string fileName);
void writeSolutionToCsv(VectorD &sol);
void writeSolutionToCsv(VectorLD &sol);

std::string createResultFileName(std::string filename, std::string type, std::string method);


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

  // TODO Read the file containing the matrix data
  // pass the data to an array and then pass this array to the contstructor
  // for the class Matrix.
  // Then dependeding on the type flag and method flag select the corresponding
  // method

  // Check type flag
  if (flagType == "float") {
    std::vector<float> data;
    int m, n;
    readMatrixCsv(fileName, data, m, n);
    int num_singvals = std::min(m, n);
    MatrixF A(m, n, data.data());
    VectorF w(num_singvals);

    if (flagMethod == "gr") {
      // Golub Reinsch method
      GolubReinsch_svdF(A, w);
      std::string solFileName = createResultFileName(fileName, flagType, flagMethod);
      writeSolutionToCsv(w, solFileName);

    // } else if (flagMethod == "petermethod") {
    //   VectorD w = Solver3_main(A);
    //   writeSolutionToCsv(w);

    // } else if (flagMethod == "gdg") {

    // } else {
    //   std::cout << "Error -method=<METHOD>" << std::endl;
    // }

  } else if (flagType == "double") {
    std::vector<double> data;
    int m, n;
    readMatrixCsv(fileName, data, m, n);
    MatrixD A(m, n, data.data());

  } else if (flagType == "long double") {
    std::vector<long double> data;
    int m, n;
    readMatrixCsv(fileName, data, m, n);
    MatrixLD A(m, n, data.data());
  } else {
    std::cout << "Error -type=<TYPE>" << std::endl;
    exit(EXIT_FAILURE);
  }
  }
  exit(EXIT_SUCCESS);
}

void writeSolutionToCsv(VectorF &sol, std::string filenName){
  std::filesystem::create_directory(SOLUTIONS_PATH);
  std::string filePath = SOLUTIONS_PATH "/" + filenName;
  std::ofstream file(filePath);

  if(!file.is_open()) {
    std::cerr << "Failed to open file" << "\n";
    exit(EXIT_FAILURE);
  }

  int i;
  for(i=0; i < sol.size(); i++){
    file << sol[i];
    if(i != sol.size() - 1)
      file << ",";

  }
  file.close();
}

std::string createResultFileName(std::string filename, std::string type, std::string method){
  std::string sol = filename.substr(0,filename.find("matrix.csv"));
  if(type == "float") sol += "F";
  else if (type == "double") sol += "D";
  else if (type == "long double") sol += "LD";
  else;

  sol += method;
  sol += ".csv";
  return sol;
}
