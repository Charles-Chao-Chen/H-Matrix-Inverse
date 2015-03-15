#include "zperm.hpp"

#include <assert.h>

ZorderPermute::ZorderPermute()
  : nx_(0), ny_(0), numLevel_(0) {}

ZorderPermute::ZorderPermute(int nx, int ny, int level)
  : nx_(nx), ny_(ny), numLevel_(level) {

#ifdef DEBUG
  std::cout << "Grid : " << nx_ << " x " << ny_ << std::endl;
  std::cout << " hierarchy level : " << level << std::endl;
#endif
  int N = nx_*ny_;
  int map[N];

  // Fill the mapping from the 'natural' x-y-z ordering
  int index = 0;
  int rootLevel = 0;
  BuildMapOnQuadrant(map, index, rootLevel, nx_, ny_);
  assert(index == N);
    
  //P = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>(N);
  P.resize(N);
  for (int i=0; i<N; i++)
    P.indices()[i] = map[i];
}

ZorderPermute::ZorderPermute(const PMatrix& Pinv)
  : nx_(0), ny_(0), numLevel_(0) {
  P = Pinv;
}

void ZorderPermute::BuildMapOnQuadrant
(int* map, int& index, int curLevel, int thisXSize, int thisYSize)
{
    if( curLevel == numLevel_-1 )
    {
        // Stamp these indices into the buffer
        for( int j=0; j<thisYSize; ++j )
        {
            int* thisRow = &map[j*nx_];
            for( int i=0; i<thisXSize; ++i )
                thisRow[i] = index++;
        }
    }
    else
    {
        const int leftWidth = thisXSize/2;
        const int rightWidth = thisXSize - leftWidth;
        const int bottomHeight = thisYSize/2;
        const int topHeight = thisYSize - bottomHeight;

        // Recurse on the lower-left quadrant
        BuildMapOnQuadrant
        ( &map[0], index, curLevel+1,
          leftWidth, bottomHeight );
        // Recurse on the lower-right quadrant
        BuildMapOnQuadrant
        ( &map[leftWidth], index, curLevel+1,
          rightWidth, bottomHeight );
        // Recurse on the upper-left quadrant
        BuildMapOnQuadrant
        ( &map[bottomHeight*nx_], index, curLevel+1,
          leftWidth, topHeight );
        // Recurse on the upper-right quadrant
        BuildMapOnQuadrant
        ( &map[bottomHeight*nx_+leftWidth], index, curLevel+1,
          rightWidth, topHeight );
    }
}

ZorderPermute ZorderPermute::inverse() {
  return ZorderPermute( P.inverse() ); 
}

ZorderPermute ZorderPermute::transpose() {
  return inverse();
}

Eigen::MatrixXd operator*(const ZorderPermute& perm,
			  const Eigen::MatrixXd& A) {
  assert( perm.P.size() != 0);
  return (perm.P)*A;
}

Eigen::MatrixXd operator*(const Eigen::MatrixXd& A,
			  const ZorderPermute& perm) {
  assert( perm.P.size() != 0);
  return A*(perm.P);
}
