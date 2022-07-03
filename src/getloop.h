#ifndef GETLOOP_H
#define GETLOOP_H

#include <Eigen/Dense>
#include <vector>

namespace tiger
{

void dfs_get_loop2(
  int cur, int pre, 
  std::vector<bool>& vis, 
  std::vector<std::vector<int>>& G, 
  std::vector<int>& path, 
  std::vector<std::vector<int>>& loop_lst);


void boundary_loop_by_dfs2(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& bnd);
void all_boundary_loop_by_dfs2(Eigen::MatrixXd& V, Eigen::MatrixXi& F, std::vector<Eigen::VectorXi>& bnd_vec);

};




#endif