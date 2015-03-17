#ifndef iter_solver_base_hpp
#define iter_solver_base_hpp

#include "macros.hpp"

// this interface, especially the solve() funtion is inspired by
//  matlab syntax and probably Eigen iterative solver interface
class IterSolverBase {
public:
  // constructor
  IterSolverBase();
  
  virtual EMatrix solve(const EMatrix&, const EVector&) = 0;

  /*
  //EMatrix solve(const EMatrix&, const EMatrix&, );
  EMatrix solve(const EMatrix&, const EMatrix&, const double);
  EMatrix solve(const EMatrix&, const EMatrix&, const double,
		const int);
  EMatrix solve(const EMatrix&, const EMatrix&, const double,
		const int, );
  */
public:
  // stopping criteria for iterative solve
  static int    ITER_MAX_NUM;
  static double ITER_TOL;

protected:
  int num_iter;
  double time;
  double residule;
};

#endif
