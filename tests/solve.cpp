#include <iostream>
#include <vector>

#include "hmat.hpp"

void Laplacian(Eigen::MatrixXd& A, int nx, int ny);

void BuildNaturalToHierarchicalMap
( std::vector<int>& map, int xSize, int ySize, int numLevels );

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec);

  
int main(int argc, char *argv[]) {


  int nx = 4, ny = 4;
  Eigen::MatrixXd A;
  Laplacian(A, nx, ny);

  //std::cout << A() << std::endl;

  int numLevels = 2; //atoi(agrv[1]);
  std::vector<int> map;
  BuildNaturalToHierarchicalMap(map, nx, ny, numLevels);

  //std::cout << map << std::endl;

  int N = nx*ny;
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(N);
  for (int i=0; i<N; i++)
    perm.indices()[i] = map[i];
    
  Eigen::MatrixXd Aperm = perm.transpose() * A * perm;
    

  int maxRank = 10;
  AdmissType admiss = WEAK;
  HMat Ah(Aperm, maxRank, numLevels,  admiss, nx, ny);
  
  return 0;
}

void Laplacian(Eigen::MatrixXd& A, int nx, int ny) {
  int N = nx*ny;
  A = Eigen::MatrixXd::Zero( N, N );

  for (int x=0; x<nx; x++) {
    for (int y=0; y<ny; y++) {
      int s = x+y*nx;
      A(s,s) = 4;
      if (x > 0)
	A(s,s-1) = -1;
      if (x < nx-1)
	A(s,s+1) = -1;
      if (y > 0)
	A(s,s-nx) = -1;
      if (y < ny-1)
	A(s,s+nx) = -1;
    }
  }
}


void BuildMapOnQuadrant
( int* map, int& index, int level, int numLevels,
  int xSize, int ySize, int thisXSize, int thisYSize )
{
    if( level == numLevels-1 )
    {
        // Stamp these indices into the buffer
        for( int j=0; j<thisYSize; ++j )
        {
            int* thisRow = &map[j*xSize];
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
        ( &map[0], index, level+1, numLevels,
          xSize, ySize, leftWidth, bottomHeight );
        // Recurse on the lower-right quadrant
        BuildMapOnQuadrant
        ( &map[leftWidth], index, level+1, numLevels,
          xSize, ySize, rightWidth, bottomHeight );
        // Recurse on the upper-left quadrant
        BuildMapOnQuadrant
        ( &map[bottomHeight*xSize], index, level+1, numLevels,
          xSize, ySize, leftWidth, topHeight );
        // Recurse on the upper-right quadrant
        BuildMapOnQuadrant
        ( &map[bottomHeight*xSize+leftWidth], index, level+1, numLevels,
          xSize, ySize, rightWidth, topHeight );
    }
}


void BuildNaturalToHierarchicalMap
( std::vector<int>& map, int xSize, int ySize, int numLevels ){
  
    map.resize( xSize*ySize );

    // Fill the mapping from the 'natural' x-y-z ordering
    int index = 0;
    BuildMapOnQuadrant
    ( &map[0], index, 0, numLevels, xSize, ySize, xSize, ySize );

    assert(index == xSize*ySize);
}


template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
  for (int i=0; i<vec.size(); i++)
    os << vec[i] << " ";
  os << std::endl;
  return os;
}
