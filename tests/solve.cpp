#include <iostream>
#include <vector>
#include <stdio.h>

#include "hmat.hpp"
#include "zperm.hpp"
#include "timer.hpp"

// stopping criteria for iterative solve
const int    ITER_NUM = 1000;
const double ITER_TOL = 1e-10;

// form the sparse matrix for the laplacian operation using
//  five point finite difference scheme
Eigen::MatrixXd Laplacian(int nx, int ny);

// output vector object (for debugging purpose)
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec);
  
int main(int argc, char *argv[]) {

  // default grid size
  int nx = 32, ny = 32;
  int numLevels = 3;
  int maxRank = 16;
  AdmissType admiss = WEAK;
  int nRhs = 2;
  
  // parse command line arguments
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-n")) {
      nx = ny = atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "-l")) {
      numLevels = atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "-r")) {
      maxRank = atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "-rhs")) {
      nRhs = atoi(argv[++i]);
      continue;
    }
  }

  printf("\n------------ configuration ------------\n");
  printf("Problem (grid) size   : (%d x %d)\n", nx, ny);
  printf("Hierarchy level       : %d\n", numLevels);
  printf("Off-diag rank (bound) : %d\n", maxRank);
  printf("Admissibility         : %s\n", (admiss==WEAK) ? "WEAK" : "unkown");
  printf("# of right hand size  : %d\n", nRhs);
  printf("---------------------------------------\n\n");
  
  Timer t;

  // descritize lapacian operator
  t.start();
  Eigen::MatrixXd A = Laplacian(nx, ny);
  t.stop(); t.get_elapsed_time("generate A");
  
  // random right hand side
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Random( nx*ny, nRhs );

  // z-order permutation matrix
  ZorderPermute perm(nx, ny, numLevels);
  
  // reordered matrix
  t.start();
  Eigen::MatrixXd Aperm = perm*A;
  t.stop(); t.get_elapsed_time("permute A (dgemm)");

  t.start();
  Aperm = Aperm*perm.inverse();
  t.stop(); t.get_elapsed_time("permute A (dgemm)");
  
  // build hierarchical tree
  t.start();
  HMat Ah(Aperm, maxRank, numLevels, admiss, nx, ny);
  t.stop(); t.get_elapsed_time("build tree");
  
  // get accuracy and timing for the h-solver
  t.start();
  Eigen::MatrixXd x1 = Ah/rhs;
  t.stop(); t.get_elapsed_time("solver");
  std::cout << "Fast solver residule : "
	    << (Aperm*x1 - rhs).norm()
	    << std::endl;

    /*
  // get accuracy and time for the standard LU method
  t.start();
  Eigen::MatrixXd x2 = A.lu().solve( rhs );
  t.stop();  t.get_elapsed_time();
  std::cout << "Direct solver residule : "
	    << (A*x2 - rhs).norm()
	    << std::endl;

  // test the accuracy for the original problem
  //  i.e. A x = b
  //      (P*A*P') (P*x) = P*b
  Eigen::MatrixXd rhsPerm = perm*rhs;
  Eigen::MatrixXd x3 = Ah.solve( rhsPerm );
  Eigen::MatrixXd xOrig = perm.inverse()*x3;
  printf("Residule: %e\n", (A*xOrig - rhs).norm() );
  */

  int N = nx*ny;
    
  // TODO: use the solver as preconditioner for
  //  fix point iterationa and GMRES
  // (1) implement maxRank : done
  // note: the accuracy is 1.0e with less rank than 16
  // (2) implement hmat * vec
  // (3) implement fix point iteration
  

  /*
    // debugging hmat * vec
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(N, 1);
  std::cout << "correct : \n" << Aperm*x << std::endl;
  std::cout << "hmat * vec : \n" << Ah*x << std::endl;
  std::cout << "debugging hmat * vec : "
	    << (Ah*x - Aperm*x).norm()
	    << std::endl;
    */
    
  /*
  int niter = 1e5;
  Eigen::MatrixXd x_cur = Eigen::MatrixXd::Zero(N, nRhs);
  for (int i=0; i<niter; i++) {
    Eigen::MatrixXd r = rhs - Aperm * x_cur;
    Eigen::MatrixXd del = Ah/r ;
    std::cout << "residule : " << r.norm() << std::endl;
    x_cur += del;

    if (r.norm() <= ITER_TOL) {
      std::cout << "Converged!\n" << " iter # : "
		<<  i+1 << std::endl;
      break;
    }
  }
  */

  // TODO: use as preconditioner for CG
  // (1) implement CG : check
  // (2) implement PCG

  // right hand size
  Eigen::VectorXd b = Eigen::VectorXd::Random(N);
  Eigen::VectorXd x = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd r_cur = b;
  Eigen::VectorXd r_pre = b;
  Eigen::VectorXd p = r_cur;
  Eigen::VectorXd Ap;
  double alpha, beta;
  int j = 0;
  while (j            < ITER_NUM &&
	 r_cur.norm() > ITER_TOL ) {
    
    Ap = Aperm * p;
    alpha = r_pre.dot(r_pre) / Ap.dot(p);
    x += alpha * p;
    r_cur = r_pre - alpha * Ap;
    beta = r_cur.dot(r_cur) / r_pre.dot(r_pre);
    p = r_cur + beta * p;
    r_pre = r_cur;
    j++;
    std::cout << "residule : " << r_cur.norm() << std::endl;
  }
  std::cout << "Converged!\n" << " iter # : " <<  j << std::endl;

  
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

