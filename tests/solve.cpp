#include <iostream>
#include <vector>
#include <stdio.h>

#include "hmat.hpp"
#include "timer.hpp"

class ZorderPermute {

public:
  ZorderPermute();
  ZorderPermute(int nx, int ny, int level);
  
  // apply the permutation matrix
  Eigen::MatrixXd convert(const Eigen::MatrixXd&);
  
private:  
  void BuildMapOnQuadrant
  (int* map, int& index, int curLevel, int thisXSize, int thisYSize);
  
  // private variable
  int nx_;
  int ny_;
  int numLevel_;
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P;
};


// form the sparse matrix for the laplacian operation using
//  five point finite difference scheme
void Laplacian(Eigen::MatrixXd& A, int nx, int ny);

void BuildNaturalToHierarchicalMap
( std::vector<int>& map, int xSize, int ySize, int numLevels );

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec);
  
int main(int argc, char *argv[]) {

  int nx = 8, ny = 8;
  Eigen::MatrixXd A;
  Laplacian(A, nx, ny);

  int numLevels = 4;    
  ZorderPermute perm(nx, ny, numLevels);
  Eigen::MatrixXd Aperm = perm.convert(A);
  
  int maxRank = 10;
  AdmissType admiss = WEAK;
  HMat Ah(Aperm, maxRank, numLevels, admiss, nx, ny);

  // random right hand side
  int nRhs = 2;
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Random( nx*ny, nRhs );

  Timer t; t.start();
  Eigen::MatrixXd x1 = Ah.solve( rhs );
  t.stop(); t.get_elapsed_time();
  std::cout << "Fast solver residule : "
	    << (Aperm*x1 - rhs).norm()
	    << std::endl;
  
  t.start();
  Eigen::MatrixXd x2 = A.lu().solve( rhs );
  t.stop();  t.get_elapsed_time();
  std::cout << "Direct solver residule : "
	    << (Aperm*x1 - rhs).norm()
	    << std::endl;

  //TODO: test the accuracy for the original matrix A
  // A x = b
  // (P*A*P') (P*x) = P*b

  Eigen::MatrixXd rhsPerm = perm.convert(rhs);
  Eigen::MatrixXd x1 = Ah.solve( rhsPerm1 );
  Eigen::MatrixXd xOrig = perm.transpose_convert(x1);
  printf("Residule: %e\n", (A*xOrig - rhs).norm() );
    
  return 0;
}

void Laplacian(Eigen::MatrixXd& A, int nx, int ny) {
  int N = nx*ny;
  A = Eigen::MatrixXd::Zero( N, N );

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
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
  for (int i=0; i<vec.size(); i++)
    os << vec[i] << " ";
  os << std::endl;
  return os;
}


ZorderPermute::ZorderPermute()
  : nx_(0), ny_(0), numLevel_(0) {}

ZorderPermute::ZorderPermute(int nx, int ny, int level)
  : nx_(nx), ny_(ny), numLevel_(level) {

#ifdef DEBUG
  std::cout << "Grid : " << nx_ << " x " << ny_ << std::endl;
  std::cout << " hierarchy level : " << level << std::endl;
#endif
  int N = nx_*ny_;
  int map[N];

  // Fill the mapping from the 'natural' x-y-z ordering
  int index = 0;
  int rootLevel = 0;
  BuildMapOnQuadrant(map, index, rootLevel, nx_, ny_);
  assert(index == N);
    
  //P = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>(N);
  P.resize(N);
  for (int i=0; i<N; i++)
    P.indices()[i] = map[i];
}

Eigen::MatrixXd ZorderPermute::convert(const Eigen::MatrixXd& A) {
  assert( P.size() != 0);
  return P*A*P.transpose();
}

void ZorderPermute::BuildMapOnQuadrant
(int* map, int& index, int curLevel, int thisXSize, int thisYSize)
{
    if( curLevel == numLevel_-1 )
    {
        // Stamp these indices into the buffer
        for( int j=0; j<thisYSize; ++j )
        {
            int* thisRow = &map[j*nx_];
            for( int i=0; i<thisXSize; ++i )
                thisRow[i] = index++;
        }
    }
    else
    {
        const int leftWidth = thisXSize/2;
        const int rightWidth = thisXSize - leftWidth;
        const int bottomHeight = thisYSize/2;
        const int topHeight = thisYSize - bottomHeight;

        // Recurse on the lower-left quadrant
        BuildMapOnQuadrant
        ( &map[0], index, curLevel+1,
          leftWidth, bottomHeight );
        // Recurse on the lower-right quadrant
        BuildMapOnQuadrant
        ( &map[leftWidth], index, curLevel+1,
          rightWidth, bottomHeight );
        // Recurse on the upper-left quadrant
        BuildMapOnQuadrant
        ( &map[bottomHeight*nx_], index, curLevel+1,
          leftWidth, topHeight );
        // Recurse on the upper-right quadrant
        BuildMapOnQuadrant
        ( &map[bottomHeight*nx_+leftWidth], index, curLevel+1,
          rightWidth, topHeight );
    }
}

