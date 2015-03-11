#include "hmat.hpp"

#include <iostream>

// the enum type and NewRHS() function aim at forming the
//  right hand side for the sub-problem with
//  part of the original right hand side and the U matrix
enum SubProblem {
  TOP,
  BOTTOM,
};

template <SubProblem type>
static EMatrix NewRHS(const EMatrix& oldRHS, const EMatrix& U);

// assembles the solution for the original problem from the
//  solutions of two sub-problems
static EMatrix AssembleSolution
(const EMatrix& x0, const EMatrix& x1, int,
 const EMatrix& V0, const EMatrix& V1);

HMat::HMat()
  : maxRank_(0),       numLevels_(0),
    admissType_(WEAK), treeRoot_(NULL)
{}

HMat::HMat
(const EMatrix& A, int maxRank, int numLevels, AdmissType admiss,
 int xSize, int ySize)
  : maxRank_(maxRank),   numLevels_(numLevels),
    admissType_(admiss), treeRoot_(NULL)
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
  Point2 rootSource(0,0);          // only one source at root level
  Point2 rootTarget(0,0);          // only one target at root level
  Rect2 rootSrcSize(xSize, ySize);
  Rect2 rootTgtSize(xSize, ySize);

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
    const EMatrix& denseBlock = node->dmat();
    return denseBlock.lu().solve( rhs );
  }
  
  // solve two sub-problems,
  //  which correspond to the first and the last 2x2 matrix blocks
  // Note: the following four matrices seems to be redunant copics
  const EMatrix U0 = node->get_topU();
  const EMatrix U1 = node->get_botU();
  const EMatrix V0 = node->get_botV();
  const EMatrix V1 = node->get_topV();
  
#ifdef DEBUG
  assert( U0.rows()+U1.rows() == rhs.rows() );
#endif
  // new right hand side for two sub-problems :
  //   rhsSub0 = [ rhsTop, U0 ]
  //   rhsSub1 = [ rhsBot, U1 ]
  EMatrix rhsSub0 = NewRHS<TOP>   (rhs, U0);
  EMatrix rhsSub1 = NewRHS<BOTTOM>(rhs, U1);
    
  // solve the two diagonal 2x2 blocks
  EMatrix x0 = solve_2x2( rhsSub0, node, 0 );
  EMatrix x1 = solve_2x2( rhsSub1, node, 2 );
  return  AssembleSolution( x0, x1, rhs.cols(), V0, V1 );
}

// standard HODLR solver for the following matrix structure
//  -----------------------
//  |          |          |
//  |  dense   |  lowrank |
//  |          |          |
//  -----------------------
//  |          |          |
//  | lowrank  |   dense  |
//  |          |          |
//  -----------------------
//  refering to : http://arxiv.org/abs/1403.5337
EMatrix HMat::solve_2x2
(const EMatrix& rhs, const Node* node, int first) {

  int second = first + 1;
  const EMatrix& U0 = node->child( first, second )->umat();
  const EMatrix& U1 = node->child( second, first )->umat();
  const EMatrix& V0 = node->child( second, first )->vmat();
  const EMatrix& V1 = node->child( first, second )->vmat();

#ifdef DEBUG
  assert( U0.rows()+U1.rows() == rhs.rows() );
#endif
  // new right hand side for two sub-problems :
  //   rhsSub0 = [ rhsTop, U0 ]
  //   rhsSub1 = [ rhsBot, U1 ]
  EMatrix rhsSub0 = NewRHS<TOP>   (rhs, U0);
  EMatrix rhsSub1 = NewRHS<BOTTOM>(rhs, U1);
    
  // recursively solve
  const EMatrix x0 = solve( rhsSub0, node->child( first,  first  ) );
  const EMatrix x1 = solve( rhsSub1, node->child( second, second ) );
  return AssembleSolution( x0, x1, rhs.cols(), V0, V1 );
}

/*static*/ EMatrix  AssembleSolution
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
  int rank0 = U0_cols;
  int rank1 = U1_cols;
  int r = rank0 + rank1;
  EMatrix S = EMatrix::Identity( r, r );
  S.topRightCorner(   rank0, rank1 ) = V1 * u1;
  S.bottomLeftCorner( rank1, rank0 ) = V0 * u0;

  // Srhs = [ V1*d1 ;
  //          V0*d0 ]
  EMatrix Srhs( r, rhs_cols );
  Srhs << V1*d1, V0*d0;

  // TODO: optimize this solve
  EMatrix eta = S.lu().solve( Srhs );

  // result = [ d0-u0*eta0 ;
  //            d1-u1*eta1 ]
  EMatrix result( x0.rows()+x1.rows(), rhs_cols );
  result << d0 - u0*eta.topRows(    rank0 ),
            d1 - u1*eta.bottomRows( rank1 );
    
  return result;
}

template <SubProblem type>
/*static */ EMatrix NewRHS(const EMatrix& oldRHS, const EMatrix& U) {
#ifdef DEBUG
  assert( type==TOP || type==BOTTOM );
#endif
  EMatrix rhs( U.rows(), oldRHS.cols() + U.cols() );
  if (type == TOP) {
    // rhs = [ rhsTop, U ]
    rhs << oldRHS.topRows( U.rows() ), U;
  }
  else {
    // rhs = [ rhsBot, U ]
    rhs << oldRHS.bottomRows( U.rows() ), U;
  }
  return rhs;
}

