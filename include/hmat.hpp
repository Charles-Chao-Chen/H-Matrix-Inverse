#ifndef HMAT_HPP
#define HMAT_HPP

#include "node.hpp"


class HMat {
  
public:
  HMat();
  HMat(Eigen::MatrixXd&, int, int, AdmissType, int, int);
  
private:
  /*
  // global information of the h-matrix
  int maxRank_;
  int numLevels_;
  AdmissType admiss_;
  */
  
  Node* root;
};




#endif // HMAT_HPP
