#include "hmat.hpp"
#include <iostream>

HMat::HMat() {}

HMat::HMat
(const Eigen::MatrixXd& A, int maxRank, int numLevels, AdmissType admiss,
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

  int  rootLevel =   0;          // the root starts from level 0
  Dim2 rootSource(0,0);          // only one source at root level
  Dim2 rootTarget(0,0);          // only one target at root level
  Dim2 rootSrcSize(xSize, ySize);
  Dim2 rootTgtSize(xSize, ySize);
  
  //Dim2 rootOffset(0,0);          // no offset for the root
  //Dim2 rootSize  (xSize, ySize); // the whole matrix for the root
  //treeRoot_ = new Node(A, rootOffset, rootSize, admissType_,
  //		      rootSource, rootTarget, rootLevel, numLevels_);
  treeRoot_ = new Node(A,
		       rootSource,  rootTarget,
		       rootSrcSize, rootTgtSize,
		       admissType_, rootLevel, numLevels_);
}


HMat::~HMat() {
  DestroyNode( treeRoot_ );
}



