#include "matrix.h"
#include "methods.h"
#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>

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

int main(int argc, char *argv[]) {

  Double data1[] = {2, 3, 4, 5, 1,  1,  1,  6,  5,  6,
                   7, 7, 8, 9, 10, 11, 12, 13, 14, 15};

  Double data2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  MatrixD A(5, 4, data1);

  MatrixD A2(3, 3, data2);

  VectorD W(3);
  MatrixD U(3, 3);
  MatrixD V(3, 3);

  Householder_svdcmp(A, U, W, V);
  // Householder_svdcmp(A2, U, W, V);

  // exit(1);

  // TODO add functionality to check
  // for valid arguments
  // string s = argv[1];

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
