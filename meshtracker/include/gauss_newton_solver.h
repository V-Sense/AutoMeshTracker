#ifndef GAUSS_NEWTON_SOLVER_H
#define GAUSS_NEWTON_SOLVER_H

#include "utils.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>	

#include <Eigen/CholmodSupport>

#include <pcl/common/common.h> 
#include <pcl/PolygonMesh.h>
#include <pcl/geometry/polygon_mesh.h>

#include <iostream>
#include <vector>

// Node idx, corrspondence pairs
typedef std::vector<std::pair<
    pcl::PointXYZRGBNormal,pcl::PointXYZRGBNormal>>                   NodeCorrs;

struct GNParams
{
    GNParams(
        double alpha_rigid  = 500.0, // 1000.0
        double alpha_smooth = 500.0, //1000.0
        double alpha_point  =   0.1, //0.1 
        double alpha_plane  =   1.0) : //1.0
        m_alpha_rigid(alpha_rigid),
        m_alpha_smooth(alpha_smooth),
        m_alpha_point(alpha_point),
        m_alpha_plane(alpha_plane) 
    {}

    // weights
    double m_alpha_rigid  = 0.0;
    double m_alpha_smooth = 0.0;
    double m_alpha_point  = 0.0;
    double m_alpha_plane  = 0.0;

    void RelaxParams(void)
    {
        this->m_alpha_rigid  /= 2.0;
        this->m_alpha_smooth /= 2.0;
    }
};

class GaussNewtonSolver
{
private:
    GNParams params_;

public:
    GaussNewtonSolver(){};
    ~GaussNewtonSolver(){};

    inline void SetParams(const GNParams p)
    {
        this->params_=p;
    }

    // Solve use Gauss-Newton Iterative Minimization.
    // returns residual energy
    double solve(
        const std::vector<std::vector<int>>    &   node_node_neigh,
        const std::vector<std::vector<double>> & node_node_weights,
        const NodeCorrs                        &       constraints,
        const std::vector<std::vector<int>>    &    con_node_neigh,
        const std::vector<std::vector<double>> & con_neigh_weights,
        std::vector<Eigen::Vector3d>           &          node_pos,
        Eigen::VectorXd                        &     affine_vector);
};

#endif
