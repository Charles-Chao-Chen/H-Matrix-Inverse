#include "node.hpp"

#include <assert.h>
#include <stdio.h>


// empty constructor for debugging purpose
Node::Node() {}

Node::Node
(const EMatrix& A,
 const Dim2& source,  const Dim2& target,
 const Dim2& srcSize, const Dim2& tgtSize,
 AdmissType admissType, int curLevel, int numLevels)

  : source_ (source),  target_ (target),
    srcSize_(srcSize), tgtSize_(tgtSize)
{

#ifdef DEBUG
  printf("Build h-tree at level : %d\n", curLevel);
#endif
  
  // low rank block
  if ( Admissible( source_, target_, admissType ) ) {
    blockType = LOWRANK;
    ComputeLowRank( UMat, VMat, A );
  }

  // dense block
  else if (curLevel == numLevels-1) {
    blockType = DENSE;
    DMat = A;
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

    // size of the children block along every dimension
    Dim2 xTgtSize = tgtSize.x_bisection();
    Dim2 yTgtSize = tgtSize.y_bisection();
    Dim2 xSrcSize = srcSize.x_bisection();
    Dim2 ySrcSize = srcSize.y_bisection();
      
    // offset in the matrix, starting from row 0 and column 0
    for (int t=0, startRow=0; t<4; t++) {
      Dim2 tIdx = ZorderIdx( t );
      Dim2 tgtChild( 2*target_ + tIdx );
      Dim2 tgtSizeChild( xTgtSize[ tIdx.x_ ], yTgtSize[ tIdx.y_ ] );
	
      for (int s=0, startCol=0; s<4; s++) {
	Dim2 sIdx = ZorderIdx( s );
	Dim2 srcChild( 2*source_ + sIdx );
	Dim2 srcSizeChild( xSrcSize[ tIdx.x_ ], ySrcSize[ tIdx.y_ ] );

	int  lChild = curLevel+1;	
	EMatrix Achild
	  = A.block( startRow, startCol,
		     tgtSizeChild.area(), srcSizeChild.area() );
	children[t][s] = new Node(Achild,
				  srcChild, tgtChild,
				  srcSizeChild, tgtSizeChild,
				  admissType, lChild, numLevels);
	// update column offset
	startCol += srcSizeChild.area();
      }
      // update row offset
      startRow += tgtSizeChild.area();
    }
  }
}

const EMatrix Node::get_topU() const {
#ifdef DEBUG
  assert( this->blockType == HIERARCHY );
#endif
  const EMatrix& U0 = children[0][2]->UMat;
  const EMatrix& U1 = children[0][3]->UMat;
  const EMatrix& U2 = children[1][2]->UMat;
  const EMatrix& U3 = children[1][3]->UMat;
  return FormUfrom2x2( U0, U1, U2, U3 );
}

const EMatrix Node::get_botU() const {
#ifdef DEBUG
  assert( this->blockType == HIERARCHY );
#endif
  const EMatrix& U0 = children[2][0]->UMat;
  const EMatrix& U1 = children[2][1]->UMat;
  const EMatrix& U2 = children[3][0]->UMat;
  const EMatrix& U3 = children[3][1]->UMat;
  return FormUfrom2x2( U0, U1, U2, U3 );
}

const EMatrix Node::get_topV() const {
#ifdef DEBUG
  assert( this->blockType == HIERARCHY );
#endif
  const EMatrix& V0 = children[0][2]->VMat;
  const EMatrix& V1 = children[0][3]->VMat;
  const EMatrix& V2 = children[1][2]->VMat;
  const EMatrix& V3 = children[1][3]->VMat;
  return FormVfrom2x2( V0, V1, V2, V3 );
}

const EMatrix Node::get_botV() const {
#ifdef DEBUG
  assert( this->blockType == HIERARCHY );
#endif
  const EMatrix& V0 = children[2][0]->VMat;
  const EMatrix& V1 = children[2][1]->VMat;
  const EMatrix& V2 = children[3][0]->VMat;
  const EMatrix& V3 = children[3][1]->VMat;
  return FormVfrom2x2( V0, V1, V2, V3 );
}

