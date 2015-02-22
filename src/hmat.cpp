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

EMatrix HMat::solve( const EMatrix& rhs ) {
  std::cout << "Starting fast solver ..." << std::endl;
  return solve( rhs, treeRoot_ );
}

EMatrix HMat::solve( const EMatrix& rhs, const Node* node ) {
  if (node->is_leaf()) {
    const EMatrix& denseBlock = node->get_dmat();
    return denseBlock.lu().solve( rhs );
  }
  
  // solve two sub-problems,
  //  which correspond to the first and the last 2x2 matrix blocks

  const EMatrix U0 = node->get_topU();
  const EMatrix U1 = node->get_botU();
  const EMatrix V0 = node->get_topV();
  const EMatrix V1 = node->get_botV();
  
#ifdef DEBUG
  assert( U0.rows()+U1.rows() == rhs.rows() );
#endif
  EMatrix rhsSub0( U0.rows(), rhs.cols() + U0.cols() );
  EMatrix rhsSub1( U1.rows(), rhs.cols() + U1.cols() );
  rhsSub0.leftCols( rhs.cols() ) = rhs.topRows(    U0.rows() );
  rhsSub1.leftCols( rhs.cols() ) = rhs.bottomRows( U1.rows() );
  rhsSub0.rightCols( U0.cols() ) = U0;
  rhsSub1.rightCols( U1.cols() ) = U1;
  
  // solve the two diagonal 2x2 blocks
  const EMatrix x0 = solve_2x2( rhsSub0, node, 0 );
  const EMatrix x1 = solve_2x2( rhsSub1, node, 2 );

  // extract u0, u1 and d0, d1
  const EMatrix& d0 = x0.leftCols( rhs.cols() );
  const EMatrix& d1 = x1.leftCols( rhs.cols() );
  const EMatrix& u0 = x0.rightCols( U0.cols() );
  const EMatrix& u1 = x1.rightCols( U1.cols() );

  return RecoverSolution( x0, x1, rhs.cols(), V0, V1 );
  
    /*
  EMatrix rhsSub1 = [rhs_second_half, U1];
  

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
//  refering to : http://arxiv.org/abs/1403.5337
EMatrix HMat::solve_2x2
(const EMatrix& rhs, const Node* node, int first) {

  // -----------------------
  // |          |          |
  // |  dense   |  lowrank |
  // |          |          |
  // -----------------------
  // |          |          |
  // | lowrank  |   dense  |
  // |          |          |
  // -----------------------

  int second = first + 1;
  const EMatrix& U0 = node->child( first, second )->get_umat();
  const EMatrix& U1 = node->child( second, first )->get_umat();
  const EMatrix& V0 = node->child( second, first )->get_vmat();
  const EMatrix& V1 = node->child( first, second )->get_vmat();

#ifdef DEBUG
  assert( U0.rows()+U1.rows() == rhs.rows() );
#endif
  EMatrix rhsSub0( U0.rows(), rhs.cols() + U0.cols() );
  EMatrix rhsSub1( U1.rows(), rhs.cols() + U1.cols() );
  rhsSub0.leftCols( rhs.cols() ) = rhs.topRows(    U0.rows() );
  rhsSub1.leftCols( rhs.cols() ) = rhs.bottomRows( U1.rows() );
  rhsSub0.rightCols( U0.cols() ) = U0;
  rhsSub1.rightCols( U1.cols() ) = U1;

  // recursively solve
  const EMatrix x0 = solve( rhsSub0, node->child( first,  first  ) );
  const EMatrix x1 = solve( rhsSub1, node->child( second, second ) );

  return RecoverSolution( x0, x1, rhs.cols(), V0, V1 );
    
  /*
  // extract u0, u1 and d0, d1
  const EMatrix& d0 = x0.leftCols( rhs.cols() );
  const EMatrix& d1 = x1.leftCols( rhs.cols() );
  const EMatrix& u0 = x0.rightCols( U0.cols() );
  const EMatrix& u1 = x1.rightCols( U1.cols() );


  // solve the small system
  int rank0 = u0.cols();
  int rank1 = u1.cols();
  int r = rank0 + rank1;
  EMatrix S = EMatrix::Identity( r, r );
  S.topRightCorner(   rank0, rank1 ) = V1 * u1;
  S.bottomLeftCorner( rank1, rank0 ) = V0 * u0;

  EMatrix Srhs( r, rhs.cols() );
  Srhs.topRows(    rank1 ) = V1*d1;
  Srhs.bottomRows( rank0 ) = V0*d0;

  const EMatrix eta = S.lu().solve( Srhs );
  EMatrix result( rhs.rows(), rhs.cols() );
  result.topRows(    d0.rows() ) = d0 - u0*eta.topRows(    rank1 );
  result.bottomRows( d1.rows() ) = d1 - u1*eta.bottomRows( rank0 );
  
  return result;
  */
}


EMatrix HMat::RecoverSolution
(const EMatrix& x0, const EMatrix& x1, int rhs_cols,
 const EMatrix& V0, const EMatrix& V1) {

  int U0_cols = x0.cols() - rhs_cols;
  int U1_cols = x1.cols() - rhs_cols;
    
  // extract u0, u1 and d0, d1
  const EMatrix& d0 = x0.leftCols( rhs_cols );
  const EMatrix& d1 = x1.leftCols( rhs_cols );
  const EMatrix& u0 = x0.rightCols( U0_cols );
  const EMatrix& u1 = x1.rightCols( U1_cols );

  // solve the small system
  int rank0 = u0.cols();
  int rank1 = u1.cols();
  int r = rank0 + rank1;
  EMatrix S = EMatrix::Identity( r, r );
  S.topRightCorner(   rank0, rank1 ) = V1 * u1;
  S.bottomLeftCorner( rank1, rank0 ) = V0 * u0;

  EMatrix Srhs( r, rhs_cols );
  Srhs.topRows(    rank1 ) = V1*d1;
  Srhs.bottomRows( rank0 ) = V0*d0;

  const EMatrix eta = S.lu().solve( Srhs );
  EMatrix result( V0.cols()+V1.cols(), rhs_cols );
  result.topRows(    d0.rows() ) = d0 - u0*eta.topRows(    rank1 );
  result.bottomRows( d1.rows() ) = d1 - u1*eta.bottomRows( rank0 );

  return result;
}
