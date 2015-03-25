#ifndef HELPER_H
#define HELPER_H

#include "dim2.hpp"   // for Admissible()
#include "macros.hpp" // for Admissible() and EMatrix type

#include <vector> // for PartialSum()

// the simplist way to get low rank factorization from a sparse
//  matrix : A = sum_k A(i_k, j_k) e_{i_k} e_{j_k}^T
void ComputeLowRank
(EMatrix& UMat, EMatrix& VMat, const EMatrix& A);

void ComputeLowRank_SVD
(EMatrix& UMat, EMatrix& VMat, const EMatrix& A, const int maxRank);

// compute the partial sum of the input array
//  e.g. given : 1, 2, 3
//       ouput : 0, 1, 3, 6
// (assume 0 at the end of the array, so for the above example,
//  the actual input should be : 1, 2, 3, 0
template< typename T >
void PartialSum(T vec[], int size) {
  T sum  = 0;
  int i = 0;
  for (; i<size-1; i++) {
    T tmp = vec[i];
    vec[i] = sum;
    sum += tmp;
  }
  vec[i] = sum;
}

void PartialSum( std::vector<int>& vec);

// check admissibility
bool Admissible(Point2 source, Point2 target, AdmissType admissType);

inline Point2 ZorderIdx( int idx ) {
  return Point2( idx>>1, idx&1 );
}

#endif // HELPER_H
