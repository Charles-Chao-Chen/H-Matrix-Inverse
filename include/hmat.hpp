#ifndef HMAT_HPP
#define HMAT_HPP

#include "node.hpp"


class HMat {
  
public:
  HMat();
  HMat(const Eigen::MatrixXd&, int, int, AdmissType, int, int);
  ~HMat();

  void solve(Eigen::MatrixXd&);
  
private:

  // global information of the tree
  int maxRank_;
  int numLevels_;
  AdmissType admissType_;

  
  Node* treeRoot_;
};




#endif // HMAT_HPP
