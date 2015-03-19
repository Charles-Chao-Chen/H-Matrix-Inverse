#ifndef fixed_point_hpp
#define fixed_point_hpp

#include "solverBase.hpp"

class FixedPoint : public IterSolverBase {
public:
  // Gauss-Seidel
  EMatrix solve(const EMatrix&, const EVector&);
  EMatrix solve(const EMatrix&, const EVector&, const Pcond&);
};

#endif
