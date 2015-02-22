#ifndef NODE_HPP
#define NODE_HPP

#include "dim2.hpp"
#include "helper.hpp"
#include <Eigen/Dense>


class Node {

public:
  
  enum BlockType {
    DENSE,
    LOWRANK,
    HIERARCHY,
  };

  // empty constructor for debugging purpose
  Node();
  
  Node
  (const EMatrix& A,
   const Dim2& source,  const Dim2& target,
   const Dim2& srcSize, const Dim2& tgtSize,
   AdmissType admissType, int curLevel, int numLevels);
  
public:
  Node* child(int, int); // get child pointer
  bool  is_leaf();
  const EMatrix& get_dense_matrix() const;
  
private:
  int  level_;
  Dim2 source_;  // index  of source cell
  Dim2 target_;  // index  of target cell
  Dim2 srcSize_; // size of the source cell
  Dim2 tgtSize_; // size of the target cell
  //Dim2 blockSize_; // size of the matrix block, i.e. length and width

  Node* children[4][4]; // 4=2^2 i.e. sub-divide in 2d

  // Matrix data :
  //  UMat and VMat for low rank factorization
  //   note VMat is a fat matrix, i.e. VMat^T exactly
  //  DMat for dense block
  BlockType blockType;
  EMatrix   UMat;
  EMatrix   VMat;
  EMatrix   DMat;  
};


#endif // NODE_HPP
