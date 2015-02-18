#include "hmat.hpp"
#include <iostream>

HMat::HMat() {}

HMat::HMat
(Eigen::MatrixXd& A, int maxRank, int numLevels, AdmissType admiss,
 int xSize, int ySize)
  : maxRank_(maxRank), numLevels_(numLevels), admissType_(admiss)
{

#ifdef DEBUG
  std::cout << "Constructing tree ..." << std::endl;
  std::cout << "max rank : " << maxRank_ << std::endl;
  std::cout << "# of levels : " << numLevels_ << std::endl;
  std::cout << "admissibility : "
	    << (admissType_==STRONG ? "STRONG" : "WEAK")
	    << std::endl;
#endif

  int  rootLevel = 0;   // the root starts from level 0
  Vec2 rootSource(0,0); // only one source at root level
  Vec2 rootTarget(0,0); // only one target at root level
  Vec2 rootOffset(0,0); // no offset for the root
  Vec2 rootSize(xSize, ySize);
  treeRoot_ = new Node(A, rootOffset, rootSize, admissType_,
		      rootSource, rootTarget, rootLevel, numLevels_);
}


HMat::~HMat() {
  DestroyNode( treeRoot_ );
}



