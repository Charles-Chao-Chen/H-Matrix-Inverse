#ifndef NODE_HPP
#define NODE_HPP

#include "dim2.hpp"

#include "Eigen/Dense"
//#include "Eigen/Sparse"


enum AdmissType {
  WEAK,
  STRONG,
};

class Node {

public:
  
  enum BlockType {
    DENSE,
    LOWRANK,
    HIERARCHY,
  };

  // empty constructor for debugging purpose
  Node();

  /*
  Node(Eigen::MatrixXd&, int, AdmissType, int, int);

  Node(Eigen::MatrixXd& A, Dim2 offset, Dim2 Size, AdmissType admissType,
       Dim2 source, Dim2 target, int curLevel, int numLevels);

  Node
  (const Eigen::MatrixXd& A,
   const Dim2& offset,  const Dim2& size, const AdmissType& admissType,
   const Dim2& source,  const Dim2& target,
   int curLevel, int numLevels);
  */
  
  Node
  (const Eigen::MatrixXd& A,
   const Dim2& source,  const Dim2& target,
   const Dim2& srcSize, const Dim2& tgtSize,
   AdmissType admissType, int curLevel, int numLevels);
  
public:
  Node* child(int, int); // get child pointer
  bool is_leaf();

  
  /*
public:
  // global information of the tree
  static int maxRank;
  static int numLevels;
  static AdmissType admiss;
  
public:
  // functions for accessing static variables
  static void set_max_rank(const int);
  static void set_num_levels(const int);
  static void set_admissibility(const AdmissType);
  
  static int get_max_rank();
  static int get_num_levels();
  static AdmissType get_admissibility();
  */
  
private:
  int  level_;
  Dim2 source_;  // index  of source cell
  Dim2 target_;  // index  of target cell
  Dim2 srcSize_; // size of the source cell
  Dim2 tgtSize_; // size of the target cell
  
  //Dim2 offset_; // offset of the matrix block
  Dim2 blockSize_;   // size   of the matrix block, i.e. length and width

  Node* children[4][4]; // 4=2^2 i.e. sub-divide in 2d

  // Matrix data :
  //  UMat and VMat for low rank factorization
  //  DMat for dense block
  BlockType       blockType;
  Eigen::MatrixXd UMat;
  Eigen::MatrixXd VMat;
  Eigen::MatrixXd DMat;
  
};


/* --- helper functions --- */

void DestroyNode(Node* node);

bool Admissible( Dim2 source, Dim2 target, AdmissType admissType );

Dim2 ZorderDim2( int idx );


void ComputeLowRank
(Eigen::MatrixXd& UMat, Eigen::MatrixXd& VMat,
 const Eigen::MatrixXd& A);

void Copy
(Eigen::MatrixXd& DMat, const Eigen::MatrixXd& A);


#endif // NODE_HPP
