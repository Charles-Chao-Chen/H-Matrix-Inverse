#ifndef iter_solver_base_hpp
#define iter_solver_base_hpp

// this interface, especially the solve() funtion is inspired by
// matlab syntax and probably Eigen iterative solver interface
class IterSolverBase {
  
  virtual EMatrix solve(const EMatrix&, const EMatrix&) = 0;

  /*
  //EMatrix solve(const EMatrix&, const EMatrix&, );
  EMatrix solve(const EMatrix&, const EMatrix&, const double);
  EMatrix solve(const EMatrix&, const EMatrix&, const double,
		const int);
  EMatrix solve(const EMatrix&, const EMatrix&, const double,
		const int, );
  */
private:
  static int    ITER_MAX_NUM;
  static double ITER_TOL;
  
  int num_iter;
  double time;
  double error;
};

int    IterSolverBase::ITER_MAX_NUM = 1000;
double IterSolverBase::ITER_TOL     = 1e-10;

#endif
