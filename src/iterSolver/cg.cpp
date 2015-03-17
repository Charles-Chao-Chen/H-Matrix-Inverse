#include "cg.hpp"

#include <assert.h>
#include <iostream>

typedef ConjugateGradient CG;

EMatrix CG::solve(const EMatrix& A, const EVector& rhs) {
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
    a  = r0.dot(r0) / Ap.dot(p);
    x += a * p;
    r1 = r0 - a * Ap;
    b  = r1.dot(r1) / r0.dot(r0);
    p  = r1 + b * p;
    r0 = r1;
    j++;
#ifdef DEBUG
    std::cout << "residule : " << r1.norm() << std::endl;
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
    std::cout << "CG has NOT converge!" << std::endl;
  }
  std::cout << " iter #   : " << num_iter     << std::endl
	    << " residule : " << residule     << std::endl
	    << " time     : " << time << " s" << std::endl;
  std::cout << "========================\n" << std::endl;
  return x;
}

EMatrix CG::solve(const EMatrix& A, const EVector& rhs, const Pcond& M) {
  assert( rhs.cols() == 1 );
  EVector x = EVector::Zero( rhs.rows() );
  EVector r0 = rhs;
  EVector z0 = M / rhs;
  EVector r1 = r0;
  EVector z1 = z0;
  EVector p = z1;
  EVector Ap;
  double a, b;
  int j = 0;
  Timer t; t.start();
  while (j         < ITER_MAX_NUM &&
	 r0.norm() > ITER_TOL ) {
    
    Ap = A * p;
    a  = z0.dot(r0) / Ap.dot(p);
    x += a * p;
    r1 = r0 - a * Ap;
    z1 = M / r1;
    b  = z1.dot(r1) / z0.dot(r0);
    p  = z1 + b * p;
    r0 = r1;
    z0 = z1;
    j++;
#ifdef DEBUG
    std::cout << "residule : " << r1.norm() << std::endl;
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
