#ifndef fixed_point_hpp
#define fixed_point_hpp

#include "iterSolverBase.hpp"

class FixedPoint : public IterSolverBase {
public:
  EMatrix solve(const EMatrix&, const EVector&, const Pcond&);
};

#endif
