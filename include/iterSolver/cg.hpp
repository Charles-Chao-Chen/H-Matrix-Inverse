#ifndef cg_hpp
#define cg_hpp

#include "iterSolverBase.hpp"

class ConjugateGradient : public IterSolverBase {
public:
  EMatrix solve(const EMatrix&, const EVector&);
};

/*
class CG : pulbic iterSolverBase {
public:
  CG();
  EMatrix solve(const EMatrix&, const EMatrix&);
  EMatrix solve(const EMatrix&, const EMatrix&, );
  EMatrix solve(const EMatrix&, const EMatrix&, const double);
  EMatrix solve(const EMatrix&, const EMatrix&, const double,
		const int);
  EMatrix solve(const EMatrix&, const EMatrix&, const double,
		const int, );
private:
  static int    ITER_MAX_NUM;
  static double ITER_TOL;
  int num_iter;
  double time;
  double error;
};
*/

#endif
