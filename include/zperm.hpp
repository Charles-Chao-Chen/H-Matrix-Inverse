#ifndef z_perm_hpp
#define z_perm_hpp

#include <Eigen/Dense>

class ZorderPermute {
  typedef Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> PMatrix;
public:
  // three constructors
  //  (1) defalut constructor
  //  (2) takes the grid information and compute the permutation matrix
  //  (3) takes the permutation matrix
  ZorderPermute();
  ZorderPermute(int nx, int ny, int level);
  explicit ZorderPermute(const PMatrix&);
  
  // the inverse (transpose) permutation
  ZorderPermute inverse();
  ZorderPermute transpose();
  
  // operator overloading functions
  //  so the object can be used like a permutation matrix
  friend Eigen::MatrixXd operator*(const ZorderPermute&,
				   const Eigen::MatrixXd&);
  friend Eigen::MatrixXd operator*(const Eigen::MatrixXd&,
				   const ZorderPermute&);
private:  
  void BuildMapOnQuadrant
  (int* map, int& index, int curLevel, int thisXSize, int thisYSize);
  
  // private variable
  int nx_;
  int ny_;
  int numLevel_;
  PMatrix P;
};

#endif
