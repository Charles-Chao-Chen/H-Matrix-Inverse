#ifndef NODE_HPP
#define NODE_HPP

#include "dim2.hpp"
#include "macros.hpp"

#include <iostream>
  
enum NodeType {
  DENSE,
  LOWRANK,
  HIERARCHY,
};

// base class for DenseNode, LowrankNode and HierNode
class Node {

public :

  // empty constructor
  Node();
  virtual bool is_leaf() const = 0;
  virtual const Node* child(int, int) const;

  // this is necessarily a virtual function,
  //  otherwise the following error occurs
  //
  //   Node* p = New DenseNode()
  //   delete p; // not calling ~DenseNode()
  virtual ~Node();
public:

  // H * vector
  //virtual EMatrix multiply(const EMatrix&) const = 0;
  
  //bool is_leaf() const;

    /*
  // get child pointer
  Node* child(int, int);
  const Node* child(int, int) const;

  // for hierarchical nodes only
  EMatrix get_topU() const;
  EMatrix get_botU() const;
  EMatrix get_topV() const;
  EMatrix get_botV() const;
*/
  
  // for dense nodes only
  //const EMatrix& dmat() const;

  // for low rank nodes only
  //const EMatrix& umat() const;
  //const EMatrix& vmat() const;
  
  //private:
  //int  level_;
  //Point2 source_; // index  of source cell
  //Point2 target_; // index  of target cell
  //Rect2 srcSize_; // size of the source cell
  //Rect2 tgtSize_; // size of the target cell

  //Node* children[4][4]; // 4=2^2 i.e. sub-divide each dimension in 2d

  // Matrix data :
  //  UMat and VMat for low rank factorization
  //   note VMat is a fat matrix, i.e. VMat^T exactly
  //  DMat for dense block
  //NodeType nodeType;
  //EMatrix   UMat;
  //EMatrix   VMat;
  //EMatrix   DMat;  
};

/*
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
*/

class DenseNode : public Node {
public:
  DenseNode();
  DenseNode(const EMatrix&);
  
  virtual bool is_leaf() const;
  
  //EMatrix multiply(const EMatrix&) const;
  const EMatrix& dmat() const;
private:
  EMatrix DMat;
};

class LowrankNode : public Node {
public:
  LowrankNode();
  LowrankNode(const EMatrix&, const int);

  virtual bool is_leaf() const;
  const EMatrix& umat() const;
  const EMatrix& vmat() const;
private:
  EMatrix UMat;
  EMatrix VMat;
};

class HierNode : public Node {
public:
  HierNode
  (const EMatrix& A,
   const Point2& source, const Point2& target,
   const Rect2& srcSize, const Rect2& tgtSize,
   const int curLevel,   const int numLevels,
   const int maxRank,    const AdmissType admissType);

  //Node* child(int, int);
  virtual const Node* child(int, int) const;

  virtual bool is_leaf() const;
  
  EMatrix get_topU() const;
  EMatrix get_botU() const;
  EMatrix get_topV() const;
  EMatrix get_botV() const;
private:
  Node* children[4][4]; // 4=2^2 i.e. sub-divide each dimension in 2d

  /*
  Point2 source_; // index  of source cell
  Point2 target_; // index  of target cell
  Rect2 srcSize_; // size of the source cell
  Rect2 tgtSize_; // size of the target cell
  */
};

#endif // NODE_HPP
