#ifndef HMAT_HPP
#define HMAT_HPP

#include "node.hpp"


class HMat {
  
public:
  HMat();
  HMat(const EMatrix&, int, int, AdmissType, int, int);
  ~HMat();

  EMatrix solve(const EMatrix&);
  
private:
  
  // the following three are helper functions for the solver
  EMatrix solve
  (const EMatrix&, const Node*);

  EMatrix solve_2x2
  (const EMatrix&, const Node*, int);

  EMatrix RecoverSolution
  (const EMatrix& x0, const EMatrix& x1, int,
   const EMatrix& V0, const EMatrix& V1);

  // helper for deconstructor
  void DestroyNode(Node* node);
  
  // global information of the tree
  int maxRank_;
  int numLevels_;
  AdmissType admissType_;

  // root of the tree
  Node* treeRoot_;
};




#endif // HMAT_HPP
