#include <iostream>
#include <vector>
#include <stdio.h>

#include "hmat.hpp"
#include "iterativeSolver.hpp"
#include "helperFun.hpp"

// form the sparse matrix for the laplacian operation using
//  five point finite difference scheme
Eigen::MatrixXd Laplacian(int nx, int ny);

// output vector object (for debugging purpose)
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec);

int main(int argc, char *argv[]) {
  
  // default parameters
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
  t.stop(); t.show_elapsed_time("generate A");
  
  // random right hand side
  int N = nx*ny;
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Random( N, nRhs );

  // z-order permutation matrix
  ZorderPermute perm(nx, ny, numLevels);
  
  // reordered matrix
  t.start();
  Eigen::MatrixXd Aperm = perm*A;
  t.stop(); t.show_elapsed_time("permute A (dgemm)");

  t.start();
  Aperm = Aperm*perm.inverse();
  t.stop(); t.show_elapsed_time("permute A (dgemm)");
  
  // build hierarchical tree
  t.start();
  HMat Ah(Aperm, maxRank, numLevels, admiss, nx, ny);
  t.stop(); t.show_elapsed_time("build tree");

  Eigen::MatrixXd x = Ah/rhs;
  std::cout << "residule : " << (Aperm*x - rhs).norm() << std::endl;
  
  /*
  // TODO: use as preconditioner for GMRES
  // (1) implement GMRES with O(n^2) QR
  Eigen::VectorXd b = Eigen::VectorXd::Random(N);
  Eigen::VectorXd x = Eigen::VectorXd::Zero(N);
  double beta = b.norm();
  Eigen::MatrixXd Q0(N, 1);
  Q0.cols(0) = b / beta;
  for (int n=0; n<N; n++) {
  v = Aperm * Q0.rightCols(1);
  for (int j=0; j<n; j++) {
  h = Q.col(j).dot(v);
  v -= h*Q.col(j); 
  }
  h_n+1 = v.norm();
  -++++++
  }
  */

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

