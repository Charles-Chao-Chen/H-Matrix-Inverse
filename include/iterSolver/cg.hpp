#ifndef cg_hpp
#define cg_hpp

#include "iterSolverBase.hpp"

class ConjugateGradient : public IterSolverBase {
public:
  EMatrix solve(const EMatrix&, const EVector&);
  EMatrix solve(const EMatrix&, const EVector&, const Pcond&);
};

#endif
