#include "iterSolverBase.hpp"

IterSolverBase::IterSolverBase()
  : num_iter(-1), time(-1), residule(-1) {}

int    IterSolverBase::ITER_MAX_NUM = 1000;
double IterSolverBase::ITER_TOL     = 1e-10;
