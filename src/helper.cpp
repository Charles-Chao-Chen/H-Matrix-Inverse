#include "helper.hpp"

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
(EMatrix& UMat, EMatrix& VMat, const EMatrix& A, double eps) {
  
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
  if ( S(0) < eps ) {
    UMat = EMatrix::Zero( A.rows(), 0 );
    VMat = EMatrix::Zero( 0, A.cols() );
  }
  else {
    U.col(0) *= S(0);
    int i=1;
    for (; i<S.size(); i++) {
      if (S(i)/S(i-1) > eps)
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

// The low rank prepresentation for 2x2 low rank blocks
//    -----------------------
//    |          |          |
//    |  u0*v0T  |  u1*v1T  |
//    |          |          |
//    -----------------------
//    |          |          |
//    |  u2*v2T  |  u3*v3T  |
//    |          |          |
//    -----------------------
//  can be either of the following two forms:
//
//    (1)
//    |    |    |    |    |
//    |    |    |    |    |
//    | u0 |    | u1 |    |
//    |    |    |    |    |    |     v0T     |             |
//    |    |    |    |    |    |     v2T     |             |
//    --------------------- *
//    |    |    |    |    |    |             |     v1T     |
//    |    |    |    |    |    |             |     v3T     |
//    |    | u2 |    | u3 |
//    |    |    |    |    |
//    |    |    |    |    |
//
//    (2)
//    |    |    |    |    |
//    |    |    |    |    |
//    | u0 | u1 |    |    |
//    |    |    |    |    |    |     v0T     |             |
//    |    |    |    |    |    |             |     v1T     |
//    --------------------- *
//    |    |    |    |    |    |     v2T     |             |
//    |    |    |    |    |    |             |     v3T     |
//    |    |    | u2 | u3 |
//    |    |    |    |    |
//    |    |    |    |    |
//
// Here we adopt the first case for the following two functions.

EMatrix FormUfrom2x2
(const EMatrix& U0, const EMatrix& U1,
 const EMatrix& U2, const EMatrix& U3) {

#ifdef DEBUG
  assert( U0.rows() == U1.rows() );
  assert( U2.rows() == U3.rows() );
#endif
  long rows[] = {U0.rows(), U2.rows(), 0};
  long cols[] = {U0.cols(), U1.cols(), U2.cols(), U3.cols(), 0};
  PartialSum( rows, 3 );
  PartialSum( cols, 5 );
#ifdef DEBUG
  assert( rows[2] == U0.rows()+U2.rows() );
  assert( cols[4] == U0.cols()+U1.cols()+U2.cols()+U3.cols() );
#endif
  
  EMatrix U = EMatrix::Zero( rows[2], cols[4] );
  U.block(rows[0], cols[0], U0.rows(), U0.cols()) = U0;
  U.block(rows[0], cols[2], U1.rows(), U1.cols()) = U1;
  U.block(rows[1], cols[1], U2.rows(), U2.cols()) = U2;
  U.block(rows[1], cols[3], U3.rows(), U3.cols()) = U3;  
  return U;
}

EMatrix FormVfrom2x2
(const EMatrix& V0, const EMatrix& V1,
 const EMatrix& V2, const EMatrix& V3) {

#ifdef DEBUG
  assert( V0.cols() == V2.cols() );
  assert( V1.cols() == V3.cols() );
#endif
  long cols[] = {V0.cols(), V1.cols(), 0};
  long rows[] = {V0.rows(), V1.rows(), V2.rows(), V3.rows(), 0};
  PartialSum( cols, 3);
  PartialSum( rows, 5);
#ifdef DEBUG
  assert( cols[2] == V0.cols()+V2.cols() );
  assert( rows[4] == V0.rows()+V1.rows()+V2.rows()+V3.rows() );
#endif
  
  EMatrix V = EMatrix::Zero( rows[4], cols[2] );
  V.block(rows[0], cols[0], V0.rows(), V0.cols()) = V0;
  V.block(rows[1], cols[0], V2.rows(), V2.cols()) = V2;
  V.block(rows[2], cols[1], V1.rows(), V1.cols()) = V1;
  V.block(rows[3], cols[1], V3.rows(), V3.cols()) = V3;  
  return V;
}

bool Admissible (Point2 source, Point2 target, AdmissType admissType) {
  if ( admissType == STRONG )
    return Point2::max( Point2::abs( source - target ) ) > 1;
  else
    return source != target;
}

