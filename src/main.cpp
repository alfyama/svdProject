#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>


#include "methods_solver3a.h"
#include "golub_reinsch.h"
#include "matrix.h"


using namespace std;

struct solver {
  virtual void operator()(...) = 0;
  virtual ~solver() {}
};

struct solver1 : solver {
  virtual void operator()(...) override {}
};

// Match string from argv[1]
// to int so it can be used with switch
enum SolverOpt {
  solver1,
  // solver 2
  // ...
  INVALID_SOLVER
};

SolverOpt stringArgOpt(string s);
VectorD Solver3_main(MatrixD bidiagmatrix);




int main(int argc, char *argv[]) {

  std::string flagMethod;
  std::string flagType;
  std::string fileName;

  for (int i = 1; i < argc; ++i) {
    std::string argument(argv[i]);
    if (argument.substr(0,8) == "-method=") {
      flagMethod = argument.substr(8);
    } else if (argument.substr(0,6) == "-type=") {
      fileName = argument.substr(6);
    } else if (fileName.empty()) {
        fileName = argument;
    } else {
        std::cout << "Invalid argument: " << argument << std::endl;
    }
  }

  if (flagMethod.empty() || flagType.empty() || fileName.empty() ) {
    std::cerr << "Usage: svd -method=<METHOD> -type=<TYPE> file" << std::endl;
    return 1;
  }


  std::ifstream inpData(fileName);
  if(!inpData) {

  }



  Double data1[] = {2, 3, 4, 5, 1,  1,  1,  6,  5,  6,
                    7, 7, 8, 9, 10, 11, 12, 13, 14, 15};

  Double data2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  MatrixD A(5, 4, data1);

  MatrixD A2(3, 3, data2);

  VectorD W(3);

  Householder_svdcmp(A, W);

  // unique_ptr<solver> msolver;

  // SolverOpt opt = stringArgOpt(s);

  // switch (opt) {
  //     case solver1:
  //         msolver = make_unique<struct solver1>();
  //         break;
  //     case INVALID_SOLVER:
  //         cerr << "Invalid option " << endl;
  //         // Print functionality of the whole program instead of just cerr
  //         break;

  //     default:
  //         break; // FIXME this should be properly handled

  // }

  return 0;
}

SolverOpt stringArgOpt(string s) {

  SolverOpt opt;
  // TODO use proper method for string
  // and literal comparison
  if (s == "solver1") {
    opt = solver1;
  }
  // else if (opt == "solver2") {
  //     ...
  // }
  else {
    opt = INVALID_SOLVER;
  }

  return opt;
}

VectorD Solver3_main(MatrixD bidiagmatrix)
{
    //assert that matrix is bi-diagonal

    int n = bidiagmatrix.num_cols();
    VectorD diag(n);
    VectorD supdiag(n);
    VectorD solution(n);

    // convert and square bidiagonal matrix into two vectors: diag, and supdiag
    for (int i = 0; i < n; i++)
    {
        diag[i] = bidiagmatrix(i, i) * bidiagmatrix(i, i);
    }
    for (int i = 0; i < n - 1; i++)
    {
        supdiag[i] = bidiagmatrix(i, i + 1) * bidiagmatrix(i, i + 1);
    }

    solution = Solver3_methods(diag, supdiag, n);

    return solution;
}
