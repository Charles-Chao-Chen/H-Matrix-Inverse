#include "node.hpp"


// empty constructor for debugging purpose
Node::Node() {}


// tree root constructor
Node::Node
(Eigen::MatrixXd& A, int numLevels, AdmissType admissType,
 int xSize, int ySize)
{

  if ( Admissible( source_, target_, admissType ) ) {
    
  }
    
  if (level_ < numLevels-1) { // sub-divide the domain

    
    blockType = HIERARCHY;

    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
	children[i][j] = new Node;
  }
  else {
    blockType = DENSE;
  }
}


Node::Node
(Eigen::MatrixXd& A, Vec2 offset, Vec2 size, AdmissType admissType,
 Vec2 source, Vec2 target, int curLevel, int numLevels)
  : level_(curLevel),
    source_(source), target_(target),
    offset_(offset), size_(size)
{

  if ( Admissible( source_, target_, admissType ) ) {
    
  }

  blockType = DENSE;
  
  /*
  if (level_ < numLevels-1) { // sub-divide the domain

    
    blockType = HIERARCHY;

    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
	children[i][j] = new Node;
  }
  else {
    blockType = DENSE;
  }
  */
}


bool Node::is_leaf() {
  return blockType != HIERARCHY;
}


Node* Node::child( int i, int j ) {
  return children[i][j];
}


void DestroyNode(Node* node) {
  if ( node->is_leaf() ) {
    delete node;
    return;
  }
  else {
    for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
	DestroyNode( node->child(i,j) );
  }
}


bool Admissible( Vec2 source, Vec2 target, AdmissType admissType ) {
  if ( admissType == STRONG )
    return Max( Abs( source - target ) ) > 1;
  else
    return source != target;
}


/*
// initialize class Node static variables
int Node::maxRank = -1;
int Node::numLevels = -1;
AdmissType Node::admiss = WEAK;
*/

/*static*//*
void Node::set_max_rank(const int rank) {
  Node::maxRank = rank;
}*/

/*static*//*
void Node::set_num_levels(const int level) {
  Node::numLevels = level;
}*/

/*static*//*
void Node::set_admissibility(const AdmissType admiss_) {
  Node::admiss = admiss_;
}*/

/*static*//*
int Node::get_max_rank() {return Node::maxRank;}
*/
/*static*//*
int Node::get_num_levels() {return Node::numLevels;}
*/
/*static*/
/*
AdmissType Node::get_admissibility() {return Node::admiss;}

Node::Node(Eigen::MatrixXd& S, int numLevels, int xSize, int ySize)
  : {

  

}
*/  
