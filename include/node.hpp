#ifndef NODE_HPP
#define NODE_HPP

#include "vec2.hpp"

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
  
  // tree root constructor
  Node(Eigen::MatrixXd&, int, AdmissType, int, int);

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
  Vec2 source_; // index of source cell
  Vec2 target_; // index of target cell
  Vec2 offset_; // offset of the matrix block
  Vec2 size_;   // size of the matrix block, i.e. height and width

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

bool Admissible( Vec2 source, Vec2 target, AdmissType admissType );



#endif // NODE_HPP
