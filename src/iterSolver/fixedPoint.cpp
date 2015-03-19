#include "fixedPoint.hpp"

#include <assert.h>
#include <iostream>

typedef FixedPoint FP;

// Gauss-Seidel
EMatrix FP::solve(const EMatrix& A, const EVector& rhs) {
  EMatrix L = A.triangularView<Eigen::Lower>();
  EMatrix U = A.triangularView<Eigen::StrictlyUpper>();
  EMatrix x = rhs;
  int j=0;
  Timer t; t.start();
  while (j        < ITER_MAX_NUM &&
	 (rhs - U * x).norm() > ITER_TOL ) {
    
    x = A.triangularView<Eigen::Lower>().solve(rhs - U * x);
    #if true
    std::cout << "residule : " << (rhs - U * x).norm() << std::endl;
#endif
  }
  t.stop();
  this->num_iter = j;
  this->residule = (rhs - U * x).norm();
  this->time     = t.get_elapsed_time();
  std::cout << "\n========================" << std::endl;
  if (num_iter < ITER_MAX_NUM) {
    std::cout << "Fixed point iteration has converged!"  << std::endl;
  }
  else {
    std::cout << "Fixed point iteration has NOT converge!" << std::endl;
  }
  std::cout << " iter #   : " << num_iter     << std::endl
	    << " residule : " << residule     << std::endl
	    << " time     : " << time << " s" << std::endl;
  std::cout << "========================\n" << std::endl;
  return x;  
}

EMatrix FP::solve(const EMatrix& A, const EVector& rhs, const Pcond& M) {
  assert( rhs.cols() == 1 );
  EVector x = EVector::Zero( rhs.rows() );
  EVector r = rhs;
  EVector del;
  int j=0;
  Timer t; t.start();
  while (j        < ITER_MAX_NUM &&
	 r.norm() > ITER_TOL ) {

    del = M / r ;
    x  += del;
    r   = rhs - A * x;
    j  ++;
#ifdef DEBUG
    std::cout << "residule : " << r.norm() << std::endl;
#endif
  }
  t.stop();
  this->num_iter = j;
  this->residule = r.norm();
  this->time     = t.get_elapsed_time();
  std::cout << "\n========================" << std::endl;
  if (num_iter < ITER_MAX_NUM) {
    std::cout << "Fixed point iteration has converged!"  << std::endl;
  }
  else {
    std::cout << "Fixed point iteration has NOT converge!" << std::endl;
  }
  std::cout << " iter #   : " << num_iter     << std::endl
	    << " residule : " << residule     << std::endl
	    << " time     : " << time << " s" << std::endl;
  std::cout << "========================\n" << std::endl;
  return x;
}
