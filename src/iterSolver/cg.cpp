#include "cg.hpp"

EMatrix ConjugateGradient::solve(const EMatrix& A, const EMatrix& b) {
  Eigen::VectorXd r_cur = b;
  Eigen::VectorXd r_pre = b;
  Eigen::VectorXd p = r_cur;
  Eigen::VectorXd Ap;
  double alpha, beta;
  int j = 0;
  while (j            < ITER_NUM &&
	 r_cur.norm() > ITER_TOL ) {
    
    Ap = A * p;
    alpha = r_pre.dot(r_pre) / Ap.dot(p);
    x += alpha * p;
    r_cur = r_pre - alpha * Ap;
    beta = r_cur.dot(r_cur) / r_pre.dot(r_pre);
    p = r_cur + beta * p;
    r_pre = r_cur;
    j++;
    std::cout << "residule : " << r_cur.norm() << std::endl;
  }
  std::cout << "Converged!\n" << " iter # : " <<  j << std::endl;
}
