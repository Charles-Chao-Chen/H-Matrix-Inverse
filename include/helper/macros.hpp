#ifndef MACROS_H
#define MACROS_H

#include <typeinfo> // for typeid

#include <Eigen/Dense>
typedef Eigen::MatrixXd EMatrix;
typedef Eigen::VectorXd EVector;

// constants used in ComputeLowRank_SVD()
#include <float.h>
const double SVD_ZERO_TOL = DBL_EPSILON * 1e2;
const double SVD_RANK_TOL = DBL_EPSILON * 1e3;

// admissibility type
enum AdmissType {
  WEAK,
  STRONG,
};

// error message
#include <cstdlib> // for EXIT_FAILURE
#include <cassert>
#define ErrorMessage(msg) {			\
    std::cerr					\
      << "Error in file : " << __FILE__		\
      << ", function : " << __func__		\
      << ", line : " << __LINE__		\
      << "\n\t" << msg << std::endl;		\
    assert(false);				\
  }
  //exit(EXIT_FAILURE)

#endif // MACROS_H
