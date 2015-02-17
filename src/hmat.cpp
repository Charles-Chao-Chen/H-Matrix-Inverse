#include "hmat.hpp"
#include <iostream>

HMat::HMat() {}

HMat::HMat
(Eigen::MatrixXd& S, int numLevels, int maxRank, AdmissType admiss,
 int xSize, int ySize)
//: maxRank_(maxRank), numLevels_(numLevels), admiss_(admiss)
{

  Node::set_max_rank( maxRank );
  Node::set_num_levels( numLevels );
  Node::set_admissibility( admiss );

#ifdef DEBUG
  std::cout << "Constructing tree ..." << std::endl;
  std::cout << "max rank : " << Node::get_max_rank() << std::endl;
  std::cout << "# of levels : " << Node::get_num_levels() << std::endl;
  std::cout << "admissibility : " << Node::get_admissibility() << std::endl;
#endif

  root = new Node(S, numLevels, xSize, ySize);
}


