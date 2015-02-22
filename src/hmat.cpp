#include "hmat.hpp"
#include <Eigen/LU>

#include <iostream>

HMat::HMat() {}

HMat::HMat
(const EMatrix& A, int maxRank, int numLevels, AdmissType admiss,
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

void HMat::DestroyNode(Node* node) {
  if ( node->is_leaf() ) {
    delete node;
    return;
  }
  else {
    for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
	DestroyNode( node->child(i,j) );
    
    delete node;
  }
}

// this is in-place solve
void HMat::solve( EMatrix& rhs ) {
  std::cout << "Starting fast solver ..." << std::endl;
  solve( rhs, treeRoot_ );
}

EMatrix HMat::solve( EMatrix& rhs, Node* node ) {
  if (node->is_leaf()) {
    const EMatrix& denseBlock = node->get_dense_matrix();
    //return denseBlock.PartialPivLU().solveInPlace( rhs );
    return denseBlock.lu().solve( rhs );
  }

  return EMatrix::Identity(2,2);
  
  // solve two 2x2 sub-problems
  /*
  EMatrix rhsSub0 = [rhs_first_half,  U0];
  EMatrix rhsSub1 = [rhs_second_half, U1];
  
  solve_2x2( node_first_2x2, rhsSub0 );
  solve_2x2( node_last_2x2,  rhsSub1 );

  // extract u0, u1 and d0, d1
  EMatrix d0 = rhsSub0(:, 1); // assume one column in rhs
  EMatrix d1 = rhsSub0(:, 1); // assume one column in rhs
  EMatrix u0 = rhsSub0(:, 2:end); // assume one column in rhs
  EMatrix u1 = rhsSub0(:, 2:end); // assume one column in rhs
  
  // recover the solution of the original problem

  // solve the small system
  int r0 = u0.cols();
  int r1 = u1.cols();
  EMatrix S = EMatrix::Idendity( r0+r1 );
  S.topRight( r0, r0 ) = V1^T * u1;
  S.bottomLeft( r0, r0 ) = V0^T * u0;

  EMatrix Srhs = [V1^T*d1, V0^T*d0];
  EMatrix [eta0, eta1] = S.lu().solve( Srhs );

  return [d0 - u0*eta0; d1 - u1*eta1];
  */
}

// standard HODLR solver
void HMat::solve_2x2( Node* node, EMatrix& rhs ) {
  /*
  solve_2x2( node_first_2x2, rhsSub0 );
  solve_2x2( node_last_2x2,  rhsSub1 );

  // extract u0, u1 and d0, d1
  EMatrix d0 = rhsSub0(:, 1); // assume one column in rhs
  EMatrix d1 = rhsSub0(:, 1); // assume one column in rhs
  EMatrix u0 = rhsSub0(:, 2:end); // assume one column in rhs
  EMatrix u1 = rhsSub0(:, 2:end); // assume one column in rhs
  
  // recover the solution of the original problem

  // solve the small system
  int r0 = u0.cols();
  int r1 = u1.cols();
  EMatrix S = EMatrix::Idendity( r0+r1 );
  S.topRight( r0, r0 ) = V1^T * u1;
  S.bottomLeft( r0, r0 ) = V0^T * u0;

  EMatrix Srhs = [V1^T*d1, V0^T*d0];
  EMatrix [eta0, eta1] = S.lu().solve( Srhs );

  return [d0 - u0*eta0; d1 - u1*eta1];
  */
}
