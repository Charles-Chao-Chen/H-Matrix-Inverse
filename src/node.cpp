#include "node.hpp"

// initialize class Node static variables
int Node::maxRank = -1;
int Node::numLevels = -1;
AdmissType Node::admiss = WEAK;


/*static*/
void Node::set_max_rank(const int rank) {
  Node::maxRank = rank;
}

/*static*/
void Node::set_num_levels(const int level) {
  Node::numLevels = level;
}

/*static*/
void Node::set_admissibility(const AdmissType admiss_) {
  Node::admiss = admiss_;
}

/*static*/
int Node::get_max_rank() {return Node::maxRank;}

/*static*/
int Node::get_num_levels() {return Node::numLevels;}

/*static*/
AdmissType Node::get_admissibility() {return Node::admiss;}

Node::Node(Eigen::MatrixXd& S, int numLevels, int xSize, int ySize)
  : {

  

}
  
