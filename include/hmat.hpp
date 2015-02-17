#ifndef HMAT_HPP
#define HMAT_HPP

#include "node.hpp"


class HMat {
  
public:
  HMat();
  HMat(Eigen::MatrixXd&, int, int, AdmissType, int, int);

  ~HMat();
  
private:

  // global information of the tree
  int maxRank_;
  int numLevels_;
  AdmissType admissType_;

  
  Node* treeRoot;
};




#endif // HMAT_HPP
