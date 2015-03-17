#include "cg.hpp"
#include "timer.hpp"

#include <assert.h>
#include <iostream>

EMatrix ConjugateGradient::solve(const EMatrix& A, const EVector& rhs) {
  assert( rhs.cols() == 1 );
  EVector x = EVector::Zero( rhs.rows() );
  EVector r0 = rhs;
  EVector r1 = rhs;
  EVector p = r0;
  EVector Ap;
  double a, b;
  int j = 0;
  Timer t; t.start();
  while (j         < ITER_MAX_NUM &&
	 r0.norm() > ITER_TOL ) {
    
    Ap = A * p;
    a  = r1.dot(r1) / Ap.dot(p);
    x += a * p;
    r0 = r1 - a * Ap;
    b  = r0.dot(r0) / r1.dot(r1);
    p  = r0 + b * p;
    r1 = r0;
    j++;
#ifdef DEBUG
    std::cout << "residule : " << r0.norm() << std::endl;
#endif
  }
  t.stop();
  this->num_iter = j;
  this->residule = r1.norm();
  this->time     = t.get_elapsed_time();
  std::cout << "\n========================" << std::endl;
  if (num_iter < ITER_MAX_NUM) {
    std::cout << "CG has converged!"  << std::endl;
  }
  else {
    std::cout << "CG does NOT converge!" << std::endl;
  }
  std::cout << " iter #   : " << num_iter     << std::endl
	    << " residule : " << residule     << std::endl
	    << " time     : " << time << " s" << std::endl;
  std::cout << "========================\n" << std::endl;
  return x;
}
