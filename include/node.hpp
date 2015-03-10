#ifndef NODE_HPP
#define NODE_HPP

#include "dim2.hpp"
#include "helper.hpp"

#include <iostream>

class Node {

public :
  
  enum BlockType {
    DENSE,
    LOWRANK,
    HIERARCHY,
  };

  // empty constructor
  Node();
  
  Node
  (const EMatrix& A,
   const Point2& source,  const Point2& target,
   const Rect2& srcSize, const Rect2& tgtSize,
   AdmissType admissType, int curLevel, int numLevels);
  
public:
  
  bool is_leaf() const;
  
  Node* child(int, int); // get child pointer
  const Node* child(int, int) const;

  // for hierarchical nodes only
  const EMatrix  get_topU() const;
  const EMatrix  get_botU() const;
  const EMatrix  get_topV() const;
  const EMatrix  get_botV() const;

  // for dense nodes only
  const EMatrix& dmat() const;

  // for low rank nodes only
  const EMatrix& umat() const;
  const EMatrix& vmat() const;
  
private:
  int  level_;
  Point2 source_; // index  of source cell
  Point2 target_; // index  of target cell
  Rect2 srcSize_; // size of the source cell
  Rect2 tgtSize_; // size of the target cell

  Node* children[4][4]; // 4=2^2 i.e. sub-divide each dimension in 2d

  // Matrix data :
  //  UMat and VMat for low rank factorization
  //   note VMat is a fat matrix, i.e. VMat^T exactly
  //  DMat for dense block
  BlockType blockType;
  EMatrix   UMat;
  EMatrix   VMat;
  EMatrix   DMat;  
};


inline bool Node::is_leaf() const {
  return blockType != HIERARCHY;
}

inline Node* Node::child( int i, int j ) {
  return children[i][j];
}

inline const Node* Node::child( int i, int j ) const {
  return children[i][j];
}

inline const EMatrix& Node::dmat() const {
#ifdef DEBUG
  assert( this->blockType == DENSE );
#endif
  return DMat;
}

inline const EMatrix& Node::umat() const {
#ifdef DEBUG
  assert( this->blockType == LOWRANK );
#endif
  return UMat;
}

inline const EMatrix& Node::vmat() const {
#ifdef DEBUG
  assert( this->blockType == LOWRANK );
#endif
  return VMat;
}


#endif // NODE_HPP
