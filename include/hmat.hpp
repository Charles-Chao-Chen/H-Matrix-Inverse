#ifndef HMAT_HPP
#define HMAT_HPP

#include "node.hpp"

class HMat {
  
public:
  HMat();
  
  // This constructor takes a matrix with z-ordering
  //  and construct the hierarchical tree structure.
  //  Usage:
  //   HMat (A, maxRank, numLevels, admiss, xSize, ySize)
  //  where A is the matrix and xSize and ySize are the grid sizes
  HMat(const EMatrix&, int, int, AdmissType, int, int);

  // destroy the hierarchical tree
  ~HMat();

  // solve the linear system with the right hand side as the input
  //  it also accepts a single right hand size as an eigen vector
  EMatrix solve(const EMatrix&);
  
private:
  
  // the following three are helper functions for the solver
  EMatrix solve
  (const EMatrix&, const Node*);

  EMatrix solve_2x2
  (const EMatrix&, const Node*, int);

  // helper function for the destructor
  void DestroyNode(Node* node);
  
  // global information of the tree
  int maxRank_;
  int numLevels_;
  AdmissType admissType_;

  // root of the hierarchical tree
  Node* treeRoot_;
};


#endif // HMAT_HPP
