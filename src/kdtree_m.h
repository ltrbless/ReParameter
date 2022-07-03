#ifndef _KDTREE_M_H_
#define _KDTREE_M_H_

#include "kdtree.h"
#include "Eigen/Dense"
#include <vector>

typedef Eigen::Vector3d vec3d;
typedef Eigen::Vector2d vec2d;
typedef Eigen::Vector3i vec3i;
typedef Eigen::Vector2i vec2i;
typedef Eigen::MatrixXd mat;

class KdTreeM{

public:

    kdtree* kd;
    kdres* res;


    KdTreeM(int dimension);
    ~KdTreeM();

    bool Insert3DNode(vec3d& point, int id);
    std::vector<int> Query3DNodeByDistance(vec3d& point, double dis);
    
    bool Insert2DNode(vec2d& point, int id);
    std::vector<int> Query2DNodeByDistance(vec2d& point, double dis); 
};

#endif
