#include <iostream>
#include <vector>
#include <stdio.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/extensions/HelperMacros.h>

#include "hmat.hpp"
#include "iterativeSolver.hpp"
#include "helperFun.hpp"


// form the sparse matrix for the laplacian operation using
//  five point finite difference scheme
Eigen::MatrixXd Laplacian(int nx, int ny);

// output vector object (for debugging purpose)
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec);



class SolverTest: public CppUnit::TestCase {
  
  /*----------------Creating a Test Suite----------------------*/
  CPPUNIT_TEST_SUITE(SolverTest);
  CPPUNIT_TEST(direct_solve);
  CPPUNIT_TEST_SUITE_END();

public:
  SolverTest(): CppUnit::TestCase("Test HMat solver") {}

  void direct_solve() {
  
    // default parameters
    int nx = 32, ny = 32;
    int numLevels = 3;
    int maxRank = 16;
    AdmissType admiss = WEAK;
    int nRhs = 2;  

    // descritize lapacian operator
    Eigen::MatrixXd A = Laplacian(nx, ny);
  
    // random right hand side
    int N = nx*ny;
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Random( N, nRhs );

    // z-order permutation matrix
    ZorderPermute perm(nx, ny, numLevels);
  
    // reordered matrix
    Eigen::MatrixXd Aperm = perm*A;
    Aperm = Aperm*perm.inverse();

    // build hierarchical tree
    HMat Ah(Aperm, maxRank, numLevels, admiss, nx, ny);

    // get accuracy and timing for the h-solver
    Eigen::MatrixXd x1 = Ah/rhs;
    
    CPPUNIT_ASSERT( (Aperm*x1 - rhs).norm() < 1e-10 );
  }

};


int main(int argc, char *argv[]) {
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(SolverTest::suite());
  runner.run();
  return 0;
}

Eigen::MatrixXd Laplacian(int nx, int ny) {
  int N = nx*ny;
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero( N, N );

  for (int x=0; x<nx; x++) {
    for (int y=0; y<ny; y++) {
      int s = x+y*nx;
      A(s,s) = 4;
      if (x > 0)
	A(s,s-1) = -1;
      if (x < nx-1)
	A(s,s+1) = -1;
      if (y > 0)
	A(s,s-nx) = -1;
      if (y < ny-1)
	A(s,s+nx) = -1;
    }
  }
  return A;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
  for (int i=0; i<vec.size(); i++)
    os << vec[i] << " ";
  os << std::endl;
  return os;
}
