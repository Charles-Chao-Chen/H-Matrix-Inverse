#include "node.hpp"
#include "utility.hpp"

#include <assert.h>
#include <stdio.h>

EMatrix FormUfrom2x2
(const EMatrix& U0, const EMatrix& U1,
 const EMatrix& U2, const EMatrix& U3);

EMatrix FormVfrom2x2
(const EMatrix& V0, const EMatrix& V1,
 const EMatrix& V2, const EMatrix& V3);

//===================================================
//             Node class
//===================================================

// empty constructor for debugging purpose
Node::Node() {}

Node::Node(int r, int c) : mRows(r), mCols(c) {}

Node::~Node() {}

int Node::rows() const {return mRows;}

int Node::cols() const {return mCols;}

const EMatrix& Node::dmat() const {
  ErrorMessage("Illegal node type");
  EMatrix *dummy = new EMatrix;
  return  *dummy;
}
  
const EMatrix& Node::umat() const {
  ErrorMessage("Illegal node type");
  EMatrix *dummy = new EMatrix;
  return  *dummy;
}

const EMatrix& Node::vmat() const {
  ErrorMessage("Illegal node type");
  EMatrix *dummy = new EMatrix;
  return  *dummy;
}
  
const Node* Node::child(int, int) const {
  ErrorMessage("Illegal node type");
  return NULL;
}

EMatrix Node::get_topU() const {
  ErrorMessage("Illegal node type");
  return EMatrix::Zero(1, 1);
}

EMatrix Node::get_botU() const {
  ErrorMessage("Illegal node type");
  return EMatrix::Zero(1, 1);
}

EMatrix Node::get_topV() const {
  ErrorMessage("Illegal node type");
  return EMatrix::Zero(1, 1);
}

EMatrix Node::get_botV() const {
  ErrorMessage("Illegal node type");
  return EMatrix::Zero(1, 1);
}

//===================================================
//              HierNode class
//===================================================

HierNode::HierNode
(const EMatrix& A,
 const Point2& source, const Point2& target,
 const Rect2& srcSize, const Rect2& tgtSize,
 const int curLevel,   const int numLevels,
 const int maxRank,    const AdmissType admissType)
  
  : Node( tgtSize.area(), srcSize.area() )
{

#ifdef DEBUG
  printf("Build h-tree at level : %d\n", curLevel);
#endif

  // hierarchical block,
  //  sub-divide the grid into 2 x 2 blocks and
  //  the matrix will be divided into 16 pieces
  //  else {
  //nodeType = HIERARCHY;

  // Morton order is described in the following wiki page:
  //  http://en.wikipedia.org/wiki/Z-order_curve
  // Here the (grid) order looks like :
  // -----------
  // | 2  |  3 |
  // -----------
  // | 0  |  1 |
  // -----------
  // or in binary form:
  //     -------------
  // y=1 | 10  |  11 |
  //     -------------
  // y=0 | 00  |  01 |
  //     -------------
  //       x=0    x=1
  // (which happens to be the same as the natural order)

  // size of the children block along every dimension
  Rect2 xTgtSize = tgtSize.x_bisection();
  Rect2 yTgtSize = tgtSize.y_bisection();
  Rect2 xSrcSize = srcSize.x_bisection();
  Rect2 ySrcSize = srcSize.y_bisection();
      
  // offset in the matrix, starting from row 0 and column 0
  for (int t=0, startRow=0; t<4; t++) {
    Point2 tIdx = ZorderIdx( t );
    Point2 tgtChild( 2*target + tIdx );
    Rect2  tgtSizeChild( xTgtSize[ tIdx.x_ ], yTgtSize[ tIdx.y_ ] );
	
    for (int s=0, startCol=0; s<4; s++) {
      Point2 sIdx = ZorderIdx( s );
      Point2 srcChild( 2*source + sIdx );
      Rect2  srcSizeChild( xSrcSize[ tIdx.x_ ], ySrcSize[ tIdx.y_ ] );

      int     lChild = curLevel+1;	
      EMatrix Achild = A.block( startRow, startCol,
				tgtSizeChild.area(),
				srcSizeChild.area() );

      if (  Admissible( srcChild, tgtChild, admissType ) ) {
	children[t][s] = new LowrankNode(Achild, maxRank);
      }
      else if ( lChild == numLevels-1 ) {
	children[t][s] = new DenseNode(Achild);
      }
      else { 
	children[t][s] = new HierNode(Achild,
				      srcChild,     tgtChild,
				      srcSizeChild, tgtSizeChild,
				      lChild,       numLevels,
				      maxRank,      admissType);
      }
      // update column offset
      startCol += srcSizeChild.area();
    }
    // update row offset
    startRow += tgtSizeChild.area();
  }
}

bool HierNode::is_leaf() const { return false; }

const Node* HierNode::child( int i, int j ) const {
  return children[i][j];
}

EMatrix HierNode::get_topU() const {
  
  const EMatrix& U0 = children[0][2]->umat();
  const EMatrix& U1 = children[0][3]->umat();
  const EMatrix& U2 = children[1][2]->umat();
  const EMatrix& U3 = children[1][3]->umat();
  return FormUfrom2x2( U0, U1, U2, U3 );
}

EMatrix HierNode::get_botU() const {

  const EMatrix& U0 = children[2][0]->umat();
  const EMatrix& U1 = children[2][1]->umat();
  const EMatrix& U2 = children[3][0]->umat();
  const EMatrix& U3 = children[3][1]->umat();
  return FormUfrom2x2( U0, U1, U2, U3 );
}

EMatrix HierNode::get_topV() const {

  const EMatrix& V0 = children[0][2]->vmat();
  const EMatrix& V1 = children[0][3]->vmat();
  const EMatrix& V2 = children[1][2]->vmat();
  const EMatrix& V3 = children[1][3]->vmat();
  return FormVfrom2x2( V0, V1, V2, V3 );
}

EMatrix HierNode::get_botV() const {

  const EMatrix& V0 = children[2][0]->vmat();
  const EMatrix& V1 = children[2][1]->vmat();
  const EMatrix& V2 = children[3][0]->vmat();
  const EMatrix& V3 = children[3][1]->vmat();
  return FormVfrom2x2( V0, V1, V2, V3 );
}

EMatrix HierNode::multiply(const EMatrix& x) const {

  assert( this->cols() == x.rows() );
  EMatrix y = EMatrix::Zero( this->rows(), x.cols() );
  Node* const* child = &(children[0][0]);
  for (int r=0, rowBeg=0; r<4; r++) {
    int rowSize = (*child)->rows();
    for (int c=0, colBeg=0; c<4; c++) {
      int colSize = (*child)->cols();
      EMatrix tmp=(*child)->multiply(x.block(colBeg,0,colSize,x.cols()));
      y.block(rowBeg,0,rowSize,x.cols()) += tmp;
      colBeg += colSize;
      child  ++;
    }
    rowBeg += rowSize;
  }
  return y;
}

//===================================================
//              DenseNode class
//===================================================

DenseNode::DenseNode(const EMatrix& A)
  : Node(A.rows(), A.cols()) {
  this->DMat = A;
}

bool DenseNode::is_leaf() const { return true; }

const EMatrix& DenseNode::dmat() const {
  return DMat;
}

EMatrix DenseNode::multiply(const EMatrix& x) const {
  return DMat*x;
}

//===================================================
//              LowRankNode class
//===================================================

LowrankNode::LowrankNode(const EMatrix& A, const int maxRank)
  : Node(A.rows(), A.cols()) {
  ComputeLowRank_SVD( this->UMat, this->VMat, A, maxRank );
}

bool LowrankNode::is_leaf() const { return true; }

const EMatrix& LowrankNode::umat() const {
  return UMat;
}

const EMatrix& LowrankNode::vmat() const {
  return VMat;
}

EMatrix LowrankNode::multiply(const EMatrix& x) const {
  return UMat*(VMat*x);
}

//===================================================
//             helper functions
//===================================================

// The low rank prepresentation for 2x2 low rank blocks
//    -----------------------
//    |          |          |
//    |  u0*v0T  |  u1*v1T  |
//    |          |          |
//    -----------------------
//    |          |          |
//    |  u2*v2T  |  u3*v3T  |
//    |          |          |
//    -----------------------
//  can be either of the following two forms:
//
//    (1)
//    |    |    |    |    |
//    |    |    |    |    |
//    | u0 |    | u1 |    |
//    |    |    |    |    |    |     v0T     |             |
//    |    |    |    |    |    |     v2T     |             |
//    --------------------- *
//    |    |    |    |    |    |             |     v1T     |
//    |    |    |    |    |    |             |     v3T     |
//    |    | u2 |    | u3 |
//    |    |    |    |    |
//    |    |    |    |    |
//
//    (2)
//    |    |    |    |    |
//    |    |    |    |    |
//    | u0 | u1 |    |    |
//    |    |    |    |    |    |     v0T     |             |
//    |    |    |    |    |    |             |     v1T     |
//    --------------------- *
//    |    |    |    |    |    |     v2T     |             |
//    |    |    |    |    |    |             |     v3T     |
//    |    |    | u2 | u3 |
//    |    |    |    |    |
//    |    |    |    |    |
//
// Here we adopt the first case for the following two functions.

EMatrix FormUfrom2x2
(const EMatrix& U0, const EMatrix& U1,
 const EMatrix& U2, const EMatrix& U3) {

#ifdef DEBUG
  assert( U0.rows() == U1.rows() );
  assert( U2.rows() == U3.rows() );
#endif
  long rows[] = {U0.rows(), U2.rows(), 0};
  long cols[] = {U0.cols(), U1.cols(), U2.cols(), U3.cols(), 0};
  PartialSum( rows, 3 );
  PartialSum( cols, 5 );
#ifdef DEBUG
  assert( rows[2] == U0.rows()+U2.rows() );
  assert( cols[4] == U0.cols()+U1.cols()+U2.cols()+U3.cols() );
#endif
  
  EMatrix U = EMatrix::Zero( rows[2], cols[4] );
  U.block(rows[0], cols[0], U0.rows(), U0.cols()) = U0;
  U.block(rows[0], cols[2], U1.rows(), U1.cols()) = U1;
  U.block(rows[1], cols[1], U2.rows(), U2.cols()) = U2;
  U.block(rows[1], cols[3], U3.rows(), U3.cols()) = U3;  
  return U;
}

EMatrix FormVfrom2x2
(const EMatrix& V0, const EMatrix& V1,
 const EMatrix& V2, const EMatrix& V3) {

#ifdef DEBUG
  assert( V0.cols() == V2.cols() );
  assert( V1.cols() == V3.cols() );
#endif
  long cols[] = {V0.cols(), V1.cols(), 0};
  long rows[] = {V0.rows(), V1.rows(), V2.rows(), V3.rows(), 0};
  PartialSum( cols, 3);
  PartialSum( rows, 5);
#ifdef DEBUG
  assert( cols[2] == V0.cols()+V2.cols() );
  assert( rows[4] == V0.rows()+V1.rows()+V2.rows()+V3.rows() );
#endif
  
  EMatrix V = EMatrix::Zero( rows[4], cols[2] );
  V.block(rows[0], cols[0], V0.rows(), V0.cols()) = V0;
  V.block(rows[1], cols[0], V2.rows(), V2.cols()) = V2;
  V.block(rows[2], cols[1], V1.rows(), V1.cols()) = V1;
  V.block(rows[3], cols[1], V3.rows(), V3.cols()) = V3;  
  return V;
}
