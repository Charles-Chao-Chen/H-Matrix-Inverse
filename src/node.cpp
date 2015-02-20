#include "node.hpp"

#include <iostream>

// empty constructor for debugging purpose
Node::Node() {}


/*
// tree root constructor
Node::Node
(Eigen::MatrixXd& A, int numLevels, AdmissType admissType,
 int xSize, int ySize)
{

  if ( Admissible( source_, target_, admissType ) ) {
    
  }
    
  if (level_ < numLevels-1) { // sub-divide the domain

    
    blockType = HIERARCHY;

    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
	children[i][j] = new Node;
  }
  else {
    blockType = DENSE;
  }
}


Node::Node
(const Eigen::MatrixXd& A,
 const Dim2& offset,  const Dim2& size, const AdmissType& admissType,
 const Dim2& source,  const Dim2& target,
 int curLevel, int numLevels)

  : level_(curLevel),
    source_(source), target_(target),
    offset_(offset), size_(size)
{

  // low rank block
  if ( Admissible( source_, target_, admissType ) ) {
    blockType = LOWRANK;
    ComputeLowRank( UMat, VMat, A );
  }

  // dense block
  else if (curLevel == numLevels-1) {
    blockType = DENSE;
    Copy( DMat, A );
  }

  // hierarchical block, sub-divide the grid and
  //  the matrix block is divided into 16 pieces
  else {
    blockType = HIERARCHY;

    // Morton order is described in the following wiki page:
    //  http://en.wikipedia.org/wiki/Z-order_curve
    // Here the order looks like :
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

    Dim2 size0 = size/2;
    Dim2 size1 = size - size0;
    for (int t=0; t<4; t++) {
      Dim2 tIdx = ZorderDim2(t);
      //int tx = t >> 1;
      //int ty = t &  1;
      for (int s=0; s<4; s++) {
	Dim2 sIdx = ZorderDim2(s);
	//int sx = s >> 1;
	//int sy = s &  1;

	Dim2 OFchild;
	Dim2 SZchild();
	Dim2 SRCchild( 2*source_ + sIdx );
	Dim2 TGTchild( 2*target_ + tIdx );
	int  Lchild = curLevel-1;
	const Eigen::MatrixXd Achild;
      }
    }
  }
}
*/


Node::Node
(const Eigen::MatrixXd& A,
 const Dim2& source,  const Dim2& target,
 const Dim2& srcSize, const Dim2& tgtSize,
 AdmissType admissType, int curLevel, int numLevels)

  : source_ (source),  target_ (target),
    srcSize_(srcSize), tgtSize_(tgtSize)
{

  // low rank block
  if ( Admissible( source_, target_, admissType ) ) {
    blockType = LOWRANK;
    ComputeLowRank( UMat, VMat, A );
  }

  // dense block
  else if (curLevel == numLevels-1) {
    blockType = DENSE;
    Copy( DMat, A );
  }

  // hierarchical block,
  //  sub-divide the grid into 2 x 2 blocks and
  //  the matrix will be divided into 16 pieces
  else {
    blockType = HIERARCHY;

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

    Dim2 tgt0 = tgtSize_ / 2;
    Dim2 tgt1 = tgtSize_ - tgt0;
    Dim2 src0 = srcSize_ / 2;
    Dim2 src1 = srcSize_ - src0;

    Dim2 offset(0, 0);
    
    for (int t=0; t<4; t++) {
      Dim2 tIdx = ZorderDim2( t );
      Dim2 tgtChild( 2*target_ + tIdx );
      Dim2 tgtSizeChild( tgt0[ tIdx.x_ ], tgt1[ tIdx.y_ ] );
      offset.x_ += tgtSizeChild.area();
	
      for (int s=0; s<4; s++) {
	Dim2 sIdx = ZorderDim2( s );
	Dim2 srcChild( 2*source_ + sIdx );
	Dim2 srcSizeChild( src0[ sIdx.x_ ], src1[ sIdx.y_ ] );
	offset.y_ += srcSizeChild.area();

	int  lChild = curLevel-1;	
	Eigen::MatrixXd Achild
	  = A.block( offset.x_, offset.y_,
		     tgtSizeChild.area(), srcSizeChild.area() );
	children[t][s] = new Node(Achild,
				  srcChild, tgtChild,
				  srcSizeChild, tgtSizeChild,
				  admissType, lChild, numLevels);
      }
    }
  }
}


bool Node::is_leaf() {
  return blockType != HIERARCHY;
}


Node* Node::child( int i, int j ) {
  return children[i][j];
}


void DestroyNode(Node* node) {
  if ( node->is_leaf() ) {
    delete node;
    return;
  }
  else {
    for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
	DestroyNode( node->child(i,j) );
  }
}


bool Admissible( Dim2 source, Dim2 target, AdmissType admissType ) {
  if ( admissType == STRONG )
    return Max( Abs( source - target ) ) > 1;
  else
    return source != target;
}


void ComputeLowRank
(Eigen::MatrixXd& UMat, Eigen::MatrixXd& VMat,
 const Eigen::MatrixXd& A) {

  std::cout << "Build the low rank block ..." << std::endl;
}


void Copy
(Eigen::MatrixXd& DMat, const Eigen::MatrixXd& A) {

  std::cout << "Copy the dense block..." << std::endl;
}


Dim2 ZorderDim2( int idx ) {
  return Dim2( idx>>1, idx&1 );
}


/*
// initialize class Node static variables
int Node::maxRank = -1;
int Node::numLevels = -1;
AdmissType Node::admiss = WEAK;
*/

/*static*//*
void Node::set_max_rank(const int rank) {
  Node::maxRank = rank;
}*/

/*static*//*
void Node::set_num_levels(const int level) {
  Node::numLevels = level;
}*/

/*static*//*
void Node::set_admissibility(const AdmissType admiss_) {
  Node::admiss = admiss_;
}*/

/*static*//*
int Node::get_max_rank() {return Node::maxRank;}
*/
/*static*//*
int Node::get_num_levels() {return Node::numLevels;}
*/
/*static*/
/*
AdmissType Node::get_admissibility() {return Node::admiss;}

Node::Node(Eigen::MatrixXd& S, int numLevels, int xSize, int ySize)
  : {

  

}
*/  
