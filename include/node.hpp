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
  
  Node(int, int);
  
  // the destructor is necessarily a virtual function,
  //  otherwise the following error occurs
  //
  //   Node* p = New DenseNode()
  //   delete p; // not calling ~DenseNode()
  virtual ~Node();

  int rows() const;
  int cols() const;
  
  //------ pure virtual functions ------
  virtual bool is_leaf() const = 0;

  // H * vector
  virtual EMatrix multiply(const EMatrix&) const = 0;

  //------ overwritten in DenseNode ------
  virtual const EMatrix& dmat() const;
  
  //------ overwritten in LowrankNode ------
  virtual const EMatrix& umat() const;
  virtual const EMatrix& vmat() const;

  //------ overwritten in HierNode ------
  virtual const Node* child(int, int) const;
  virtual EMatrix get_topU() const;
  virtual EMatrix get_botU() const;
  virtual EMatrix get_topV() const;
  virtual EMatrix get_botV() const;

protected:
  int mRows;
  int mCols;
  
  //int  level_;
  //Point2 source_; // index  of source cell
  //Point2 target_; // index  of target cell
  //Rect2 srcSize_; // size of the source cell
  //Rect2 tgtSize_; // size of the target cell
};

class DenseNode : public Node {
public:
  DenseNode();
  DenseNode(const EMatrix&);
  
  virtual bool is_leaf() const;
  virtual const EMatrix& dmat() const;
  virtual EMatrix multiply(const EMatrix&) const;
private:
  EMatrix DMat;
};

class LowrankNode : public Node {
public:
  LowrankNode();
  LowrankNode(const EMatrix&, const int);

  virtual bool is_leaf() const;
  virtual const EMatrix& umat() const;
  virtual const EMatrix& vmat() const;
  virtual EMatrix multiply(const EMatrix&) const;
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

  virtual bool is_leaf() const;
  virtual const Node* child(int, int) const;
  virtual EMatrix get_topU() const;
  virtual EMatrix get_botU() const;
  virtual EMatrix get_topV() const;
  virtual EMatrix get_botV() const;
  virtual EMatrix multiply(const EMatrix&) const;
private:
  Node* children[4][4]; // 4=2^2 i.e. sub-divide each dimension in 2d
};

#endif // NODE_HPP
