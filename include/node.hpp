#ifndef NODE_HPP
#define NODE_HPP

#include "dim2.hpp"
#include "helper.hpp"

#include <iostream>

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
  bool  is_leaf() const;
  
  Node* child(int, int); // get child pointer

  const Node* child(int, int) const;

  BlockType get_block_type() const;

  // for hierarchical nodes only
  const EMatrix  get_topU() const;
  const EMatrix  get_botU() const;
  const EMatrix  get_topV() const;
  const EMatrix  get_botV() const;

  // for dense nodes only
  const EMatrix& get_dmat() const;

  // for low rank nodes only
  const EMatrix& get_umat() const;
  const EMatrix& get_vmat() const;
  
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


inline bool Node::is_leaf() const {
  return blockType != HIERARCHY;
}

inline Node* Node::child( int i, int j ) {
  return children[i][j];
}

inline const Node* Node::child( int i, int j ) const {
  return children[i][j];
}

inline const EMatrix Node::get_topU() const {
#ifdef DEBUG
  assert( node->get_block_type() == HIERARCHY );
#endif  
  std::cout << "node::get_topU() to be implemented" << std::endl;
  return EMatrix::Identity(2,2);
}

inline const EMatrix Node::get_botU() const {
#ifdef DEBUG
  assert( node->get_block_type() == HIERARCHY );
#endif
  std::cout << "node::get_botU() to be implemented" << std::endl;
  return EMatrix::Identity(2,2);
}

inline const EMatrix Node::get_topV() const {
#ifdef DEBUG
  assert( node->get_block_type() == HIERARCHY );
#endif
  std::cout << "node::get_topV() to be implemented" << std::endl;
  return EMatrix::Identity(2,2);
}

inline const EMatrix Node::get_botV() const {
#ifdef DEBUG
  assert( node->get_block_type() == HIERARCHY );
#endif
  std::cout << "node::get_botV() to be implemented" << std::endl;
  return EMatrix::Identity(2,2);
}

inline const EMatrix& Node::get_dmat() const {
#ifdef DEBUG
  assert( node->get_block_type() == DENSE );
#endif
  return DMat;
}

inline const EMatrix& Node::get_umat() const {
#ifdef DEBUG
  assert( node->get_block_type() == LOWRANK );
#endif
  return UMat;
}

inline const EMatrix& Node::get_vmat() const {
#ifdef DEBUG
  assert( node->get_block_type() == LOWRANK );
#endif
  return VMat;
}


#endif // NODE_HPP
