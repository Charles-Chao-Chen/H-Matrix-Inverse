#include "utility.hpp"

#include <iostream>

void ComputeLowRank
(EMatrix& UMat, EMatrix& VMat, const EMatrix& A) {

#ifdef DEBUG
  std::cout << "  Form the low rank block ..." << std::endl;
#endif

  int nRow = A.rows();
  int nCol = A.cols();
  int k = nRow > nCol ? nCol : nRow; // the smaller
  EMatrix Ufull = EMatrix::Zero( nRow, k );
  EMatrix Vfull = EMatrix::Zero( k, nCol );
  int rank = 0;
  for (int c=0; c<nCol; c++) {
    for (int r=0; r<nRow; r++) {
      if (A(r,c) != 0) {
#ifdef DEBUG
	assert( rank < k );
#endif
	Ufull(r,rank) = A(r,c);
	Vfull(rank,c) = 1;
	rank++;
      }
    }
  }
  if (rank>0) {
    UMat = Ufull.leftCols(rank);
    VMat = Vfull.topRows(rank);
  }
  else {
    // zero vectors
    UMat = Ufull.leftCols(1);
    VMat = Vfull.topRows(1);
  }
#ifdef DEBUG
  std::cout << "Error for low rank approximation: "
	    << (A - UMat*VMat).norm()
	    << std::endl;
#endif
}

void ComputeLowRank_SVD
(EMatrix& UMat, EMatrix& VMat, const EMatrix& A, const int maxRank) {
  
#ifdef DEBUG
  std::cout << "  Form the low rank block ..." << std::endl;
  //std::cout << A << std::endl;
#endif

  Eigen::JacobiSVD<EMatrix> svd(A, Eigen::ComputeThinU
				|  Eigen::ComputeThinV);
  Eigen::VectorXd S = svd.singularValues();
  EMatrix U = svd.matrixU();
  EMatrix V = svd.matrixV();

  // handle (numerically) zero matrix
  if ( S(0) < SVD_ZERO_TOL ) {
    UMat = EMatrix::Zero( A.rows(), 0 );
    VMat = EMatrix::Zero( 0, A.cols() );
  }
  else {
    U.col(0) *= S(0);
    int i=1;
    for (; i<S.size() && i<maxRank; i++) {
      if (S(i)/S(i-1) > SVD_RANK_TOL)
	U.col(i) *= S(i);
      else
	break;
    }
    UMat = U.leftCols(i);
    VMat = V.leftCols(i).transpose();
  }

#ifdef DEBUG
  std::cout << "Matrix size : " << S.size()
	    << " rank : " << UMat.cols()
	    << std::endl;
  std::cout << "Error for low rank approximation: "
	    << (A - UMat*VMat).norm()
	    << std::endl;
#endif
}

// compute the partial sum of the input array
//  e.g. given : 1, 2, 3
//       ouput : 0, 1, 3, 6
// (assume 0 at the end of the array, so for the above example,
//  the actual input should be : 1, 2, 3, 0
void PartialSum( std::vector<int>& vec) {
  int sum  = 0;
  int size = vec.size();
  int i = 0;
  for (; i<size-1; i++) {
    int tmp = vec[i];
    vec[i] = sum;
    sum += tmp;
  }
  vec[i] = sum;
}

bool Admissible (Point2 source, Point2 target, AdmissType admissType) {
  if ( admissType == STRONG )
    return Point2::max( Point2::abs( source - target ) ) > 1;
  else
    return source != target;
}

