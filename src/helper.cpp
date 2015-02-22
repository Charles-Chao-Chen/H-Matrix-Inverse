#include "helper.hpp"

#include <iostream>

bool Admissible( Dim2 source, Dim2 target, AdmissType admissType ) {
  if ( admissType == STRONG )
    return Max( Abs( source - target ) ) > 1;
  else
    return source != target;
}


void ComputeLowRank
(EMatrix& UMat, EMatrix& VMat,
 const EMatrix& A) {

#ifdef DEBUG
  std::cout << "  Form the low rank block ..." << std::endl;
#endif
  //JacobiSVD<MatrixXd> svd( A, ComputeThinU | ComputeThinV );

  int nRow = A.rows();
  int nCol = A.cols();
  int k = nRow > nCol ? nCol : nRow; // the smaller
  EMatrix Ufull = EMatrix::Zero( nRow, k );
  EMatrix Vfull = EMatrix::Zero( k, nCol );
  int rank = 0;
  for (int c=0; c<nCol; c++) {
    for (int r=0; r<nRow; r++) {
      if (A(r,c) != 0) {
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


void Copy
(EMatrix& DMat, const EMatrix& A) {
#ifdef DEBUG
  std::cout << "  Form the dense block..." << std::endl;
#endif
}


Dim2 ZorderIdx( int idx ) {
  return Dim2( idx>>1, idx&1 );
}
