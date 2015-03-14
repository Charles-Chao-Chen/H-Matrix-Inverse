#ifndef MACROS_H
#define MACROS_H

#include <Eigen/Dense>
typedef Eigen::MatrixXd EMatrix;

#include <float.h>
const double SVD_ZERO_TOL = DBL_EPSILON * 1e2;

const double SVD_RANK_TOL = DBL_EPSILON * 1e3;

enum AdmissType {
  WEAK,
  STRONG,
};

#endif // MACROS_H
