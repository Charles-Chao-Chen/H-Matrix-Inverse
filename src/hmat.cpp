#include "hmat.hpp"
#include <iostream>

HMat::HMat() {}

HMat::HMat
(const Eigen::MatrixXd& A, int maxRank, int numLevels, AdmissType admiss,
 int xSize, int ySize)
  : maxRank_(maxRank), numLevels_(numLevels), admissType_(admiss)
{

#ifdef DEBUG
  std::cout << "Begin building h-tree ..." << std::endl;
  std::cout << "  max rank : " << maxRank_ << std::endl;
  std::cout << "  # of levels : " << numLevels_ << std::endl;
  std::cout << "  admissibility : "
	    << (admissType_==STRONG ? "STRONG" : "WEAK")
	    << "\n\n";
#endif

  int  rootLevel  =  0;          // the root starts from level 0
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

void HMat::solve( Eigen::MatrixXd& rhs ) {
  std::cout << "Starting fast solver ..." << std::endl;
  solve( rhs, root );
}

void HMat::solve( Eigen::MatrixXd& rhs, Node* node ) {
  if (node->is_leaf()) {
    return node->DMat.PartialPivLU().solve( rhs );
  }

  // solve two 2x2 sub-problems

  Eigen::MatrixXd rhsSub0 = [rhs_first_half,  U0];
  Eigen::MatrixXd rhsSub1 = [rhs_second_half, U1];
  
  solve_2x2( node_first_2x2, rhsSub0 );
  solve_2x2( node_last_2x2,  rhsSub1 );

  // extract u0, u1 and d0, d1
  Eigen::MatrixXd d0 = rhsSub0(:, 1); // assume one column in rhs
  Eigen::MatrixXd d1 = rhsSub0(:, 1); // assume one column in rhs
  Eigen::MatrixXd u0 = rhsSub0(:, 2:end); // assume one column in rhs
  Eigen::MatrixXd u1 = rhsSub0(:, 2:end); // assume one column in rhs
  
  // recover the solution of the original problem

  // solve the small system
  int r0 = u0.cols();
  int r1 = u1.cols();
  Eigen::MatrixXd S = Eigen::MatrixXd::Idendity( r0+r1 );
  S.topRight( r0, r0 ) = V1^T * u1;
  S.bottomLeft( r0, r0 ) = V0^T * u0;

  Eigen::MatrixXd Srhs = [V1^T*d1, V0^T*d0];
  Eigen::MatrixXd [eta0, eta1] = S.lu().solve( Srhs );

  return [d0 - u0*eta0; d1 - u1*eta1];
}

// standard HODLR solver
void solve_2x2( Node* node, Eigen::MatrixXd rhs ) {

  solve_2x2( node_first_2x2, rhsSub0 );
  solve_2x2( node_last_2x2,  rhsSub1 );

  // extract u0, u1 and d0, d1
  Eigen::MatrixXd d0 = rhsSub0(:, 1); // assume one column in rhs
  Eigen::MatrixXd d1 = rhsSub0(:, 1); // assume one column in rhs
  Eigen::MatrixXd u0 = rhsSub0(:, 2:end); // assume one column in rhs
  Eigen::MatrixXd u1 = rhsSub0(:, 2:end); // assume one column in rhs
  
  // recover the solution of the original problem

  // solve the small system
  int r0 = u0.cols();
  int r1 = u1.cols();
  Eigen::MatrixXd S = Eigen::MatrixXd::Idendity( r0+r1 );
  S.topRight( r0, r0 ) = V1^T * u1;
  S.bottomLeft( r0, r0 ) = V0^T * u0;

  Eigen::MatrixXd Srhs = [V1^T*d1, V0^T*d0];
  Eigen::MatrixXd [eta0, eta1] = S.lu().solve( Srhs );

  return [d0 - u0*eta0; d1 - u1*eta1];  
}
