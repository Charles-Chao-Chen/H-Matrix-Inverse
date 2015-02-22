#ifndef HMAT_HPP
#define HMAT_HPP

#include "node.hpp"


class HMat {
  
public:
  HMat();
  HMat(const EMatrix&, int, int, AdmissType, int, int);
  ~HMat();

  void solve(EMatrix&);
  
private:

  //void solve(EMatrix&, Node*);

  EMatrix solve(EMatrix&, Node*);

  void solve_2x2(Node*, EMatrix&);

  // helper for deconstructor
  void DestroyNode(Node* node);
  
  // global information of the tree
  int maxRank_;
  int numLevels_;
  AdmissType admissType_;

  
  Node* treeRoot_;
};




#endif // HMAT_HPP
