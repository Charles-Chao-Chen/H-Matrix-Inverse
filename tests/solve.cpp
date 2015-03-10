#include <iostream>
#include <vector>
#include <stdio.h>

#include "hmat.hpp"
#include "timer.hpp"

class ZorderPermute {
  typedef Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> PMatrix;
public:
  // three constructors
  //  (1) defalut constructor
  //  (2) takes the grid information and compute the permutation matrix
  //  (3) takes the permutation matrix
  ZorderPermute();
  ZorderPermute(int nx, int ny, int level);
  explicit ZorderPermute(const PMatrix&);
  
  // the inverse (transpose) permutation
  ZorderPermute inverse();
  ZorderPermute transpose();
  
  // operator overloading functions
  //  so the object can be used like a permutation matrix
  friend Eigen::MatrixXd operator*(const ZorderPermute&,
				   const Eigen::MatrixXd&);
  friend Eigen::MatrixXd operator*(const Eigen::MatrixXd&,
				   const ZorderPermute&);
private:  
  void BuildMapOnQuadrant
  (int* map, int& index, int curLevel, int thisXSize, int thisYSize);
  
  // private variable
  int nx_;
  int ny_;
  int numLevel_;
  PMatrix P;
};

// form the sparse matrix for the laplacian operation using
//  five point finite difference scheme
void Laplacian(Eigen::MatrixXd& A, int nx, int ny);


// output vector object (for debugging purpose)
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec);
  
int main(int argc, char *argv[]) {

  int nx = 8, ny = 8;
  Eigen::MatrixXd A;
  Laplacian(A, nx, ny);

  int numLevels = 4;    
  ZorderPermute perm(nx, ny, numLevels);
  Eigen::MatrixXd Aperm = perm*A*perm.inverse();
  
  int maxRank = 10;
  AdmissType admiss = WEAK;
  HMat Ah(Aperm, maxRank, numLevels, admiss, nx, ny);

  // random right hand side
  int nRhs = 2;
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Random( nx*ny, nRhs );

  // get accuracy and time for the h-solver
  Timer t; t.start();
  Eigen::MatrixXd x1 = Ah.solve( rhs );
  t.stop(); t.get_elapsed_time();
  std::cout << "Fast solver residule : "
	    << (Aperm*x1 - rhs).norm()
	    << std::endl;

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

ZorderPermute::ZorderPermute(const PMatrix& Pinv)
  : nx_(0), ny_(0), numLevel_(0) {
  P = Pinv;
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

ZorderPermute ZorderPermute::inverse() {
  return ZorderPermute( P.inverse() ); 
}

ZorderPermute ZorderPermute::transpose() {
  return inverse();
}

Eigen::MatrixXd operator*(const ZorderPermute& perm,
			  const Eigen::MatrixXd& A) {
  assert( perm.P.size() != 0);
  return (perm.P)*A;
}

Eigen::MatrixXd operator*(const Eigen::MatrixXd& A,
			  const ZorderPermute& perm) {
  assert( perm.P.size() != 0);
  return A*(perm.P);
}

