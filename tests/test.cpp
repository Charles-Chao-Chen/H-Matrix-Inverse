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

class SolverTest: public CppUnit::TestCase {

private:
  int N;
  Eigen::MatrixXd A;
  Eigen::MatrixXd Aperm;
  ZorderPermute   perm;
  HMat           *Ah;
  HMat           *M;
  
  /*----------------Creating a Test Suite----------------------*/
  CPPUNIT_TEST_SUITE(SolverTest);
  CPPUNIT_TEST(direct_solve1);
  CPPUNIT_TEST(direct_solve2);
  CPPUNIT_TEST(hmat_vec_mul);
  CPPUNIT_TEST(cg_solve);
  CPPUNIT_TEST(fp_solve);
  CPPUNIT_TEST_SUITE_END();

public:
  SolverTest(): CppUnit::TestCase("Test HMat solver") {
    int nx = 32, ny = 32;
    int numLevels = 4;
    int maxRank = 16;
    AdmissType admiss = WEAK;

    // problem size
    this->N = nx*ny;
    
    // descritize lapacian operator
    this->A = Laplacian(nx, ny);
  
    // z-order permutation matrix
    this->perm = ZorderPermute(nx, ny, numLevels);
  
    // reordered matrix
    this->Aperm = perm*A*perm.inverse();

    // build hierarchical tree
    this->Ah = new HMat(Aperm, maxRank, numLevels, admiss, nx, ny);

    // reduced rank matrix as preconditioner
    int rank = maxRank - 2;
    this->M = new HMat(Aperm, rank, numLevels, admiss, nx, ny);
  }

  ~SolverTest() {
    delete Ah;
  }

  void direct_solve1() {
    std::cout << "Test direct solve ..." << std::endl;
    // test h-solver accuracy
    int nRhs = 2;  
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Random( N, nRhs );
    Eigen::MatrixXd x1 = Ah->solve( rhs );
    CPPUNIT_ASSERT( (Aperm*x1 - rhs).norm() < 1e-10 );
  }

  void direct_solve2() {
    std::cout << "Test direct solve ..." << std::endl;
    // test the accuracy for the original problem
    //  i.e. A x = b
    //      (P*A*P') (P*x) = P*b
    int nRhs = 2;  
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Random( N, nRhs );
    Eigen::MatrixXd rhsPerm = perm*rhs;
    Eigen::MatrixXd x3 = Ah->solve( rhsPerm );
    Eigen::MatrixXd xOrig = perm.inverse()*x3;
    CPPUNIT_ASSERT( (A*xOrig - rhs).norm() < 1e-10 );
  }

  void hmat_vec_mul() {
    std::cout << "Test hmat * vec ..." << std::endl;
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(N, 2);
    CPPUNIT_ASSERT( (Ah->multiply(x) - Aperm*x).norm() < 1e-10 );
  }

  void fp_solve() {
    std::cout << "Test fixed point iterative solve ..." << std::endl;
    Eigen::VectorXd b = Eigen::VectorXd::Random(N);
    FixedPoint fp;
    Eigen::VectorXd fp_x = fp.solve(Aperm, b, (*M));
    CPPUNIT_ASSERT( (Aperm*fp_x - b).norm() < 1e-10 );
  }

  void cg_solve() {
    std::cout << "Test CG iterative solve ..." << std::endl;
    Eigen::VectorXd b = Eigen::VectorXd::Random(N);
    ConjugateGradient cg;
    Eigen::VectorXd cg_x = cg.solve(Aperm, b, (*M));
    CPPUNIT_ASSERT( (Aperm*cg_x - b).norm() < 1e-10 );
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
