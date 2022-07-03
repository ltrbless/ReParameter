#ifndef REPARAM_H
#define REPARAM_H

#include <Eigen/Dense>

namespace tiger
{

using vec3d = Eigen::Vector3d;

void buildTuttleParameter(Eigen::MatrixXd& V_3d, Eigen::MatrixXi& T_3d, Eigen::MatrixXd& V_uv);
void buildTuttleHoleParameter(Eigen::MatrixXd& V_3d, Eigen::MatrixXi& T_3d, Eigen::MatrixXd& V_uv, Eigen::MatrixXi& T_filled);
void buildHarmonicParameter(Eigen::MatrixXd& V_3d, Eigen::MatrixXi& T_3d, Eigen::MatrixXd& V_uv);

void removeDulplicatePoint(Eigen::MatrixXd& V, Eigen::MatrixXi& T, double eps);



} // namespace tiger


#endif