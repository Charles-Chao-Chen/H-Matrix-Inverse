#ifndef HELPER_H
#define HELPER_H

#include "dim2.hpp"
#include "macros.hpp"
#include <Eigen/Dense>

bool Admissible( Dim2 source, Dim2 target, AdmissType admissType );

Dim2 ZorderIdx( int idx );

void ComputeLowRank
(EMatrix& UMat, EMatrix& VMat,
   const EMatrix& A);

void Copy
(EMatrix& DMat, const EMatrix& A);



#endif // HELPER_H
