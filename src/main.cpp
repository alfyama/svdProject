#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>
#include "matrix.h"

using namespace std;

struct solver{
    virtual void operator() (...) = 0;
    virtual ~solver() {}
};

struct solver1 : solver{
    virtual void operator() (...) override { }
};

// Match string from argv[1]
// to int so it can be used with switch
enum SolverOpt{
  solver1,
  // solver 2
  // ...
  INVALID_SOLVER
};

SolverOpt stringArgOpt(string s);

int main(int argc, char *argv[]) {

   MatrixD A(2,2);

   A(0,0) = 3.2;
   A(1,1) = 1;
   A(0,1) = 10;
   A(1,0) = 0;

   A.display();

   IdentityD B(3,3);
   B.display();

   exit(1);


    // TODO add functionality to check
    // for valid arguments
    string s = argv[1];

    unique_ptr<solver> msolver;

    SolverOpt opt = stringArgOpt(s);

    switch (opt) {
        case solver1:
            msolver = make_unique<struct solver1>();
            break;
        case INVALID_SOLVER:
            cerr << "Invalid option " << endl;
            // Print functionality of the whole program instead of just cerr
            break;

        default:
            break; // FIXME this should be properly handled

    }

    
    return 0;

}

SolverOpt stringArgOpt(string s){

    SolverOpt opt;
    // TODO use proper method for string
    // and literal comparison
    if(s == "solver1"){
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
