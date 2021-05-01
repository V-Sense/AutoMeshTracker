#ifndef DEFORM_GRAPH_H
#define DEFORM_GRAPH_H

#include <iterator>
#include <vector>
#include <algorithm>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <flann/flann.hpp>

#include <pcl/search/kdtree.h>
#include <pcl/common/common.h> 
#include <pcl/PolygonMesh.h>
#include <pcl/geometry/polygon_mesh.h>

#include "tri_mesh.h"
#include "mesh_processing.h"
#include "matching.h"
#include "gauss_newton_solver.h"
#include "utils.h"

typedef pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr                   CloudPtr;

// Node idx, corrspondence pairs
typedef std::vector<std::pair<
    pcl::PointXYZRGBNormal,pcl::PointXYZRGBNormal>>                   NodeCorrs;

class DeformGraph 
{

private:
    pcl::PolygonMesh graph_structure_;
    double node_sample_radius_ =  1.0; // init value only    
    double dmax_scaler_        =  1.5; // Scales the inclusion radius during vert-node weight assignment
    size_t debug_vertex_       = 2000;

public:

    std::vector<Eigen::Vector3d>             node_norm;
    std::vector<Eigen::Quaterniond>          node_rots;
    std::vector<Eigen::Matrix3d>         node_rots_mat;
    std::vector<Eigen::Vector3d>            node_trans;
    std::vector<Eigen::Vector3d>              node_pos;
    std::vector<std::vector<int>>           node_neigh;
    std::vector<std::vector<double>> node_node_weights;
    std::vector<std::vector<double>> vert_node_weights;

    double max_dist_neigh;      // used for building neighbours
    int           k_neigh;
    
    DeformGraph() 
    {
        k_neigh = kDG_KNN;
    };

    ~DeformGraph(){};

    void BuildGraph(
        const TriMesh    &     mesh, 
        pcl::PolygonMesh & in_graph,
        int                    k_nn = kDG_KNN);

    inline void SetMGraphStructure(pcl::PolygonMesh & mesh)
    {   
        this->graph_structure_ = mesh;
    }

    inline void SetSampleRadius(const double & radius)
    {
        this->node_sample_radius_ = radius;
    }

    void BuildNeighbours();

    void CalculateNodeNodeWeights();

    std::vector<std::vector<double>> CalculateVertNodeWeights(
        const std::vector<pcl::PointXYZRGBNormal> &      verts,
        const std::vector<std::vector<int>>       & vert_neigh);
    
    std::vector<std::vector<double>> CalculateVertNodeWeights(
        const pcl::PolygonMesh              &    in_mesh,
        const std::vector<std::vector<int>> & vert_neigh);

	void UpdateGraph(const Eigen::VectorXd & X);

    // Per-vertex get the knn nodes of influence for input mesh
    // in_graph is used as it will provide better correspondence
    std::vector<std::vector<int>> GetVertNodeNeighbours(
        pcl::PolygonMesh & in_mesh);

    // Get vert-node matches using graph connectivity and sampling radius 
    // constraints. Requires cloud so that we don't have to 
    // calculate norms erry time
    std::vector<int> GetVertNodeMatches(
        const pcl::PointXYZRGBNormal &       source, 
        const pcl::PolygonMesh       &        query,
        const CloudPtr               & query_cloud);

    std::vector<std::vector<int>> GetConNodeNeighbours(
        const std::vector<pcl::PointXYZRGBNormal> & verts,
        int                                   search_size = 12);

    // Smooth out the deform space to prevent high frequency noise
    void NormalizeTforms();

	void Deform(
        const std::vector<std::vector<int>>    & node_influence_list,
        const std::vector<std::vector<double>> &   vert_node_weights,
        const pcl::PolygonMesh                 &         input_graph,
        TriMesh                                &              d_mesh);

    // Initialize the vector of stacked affine transforms
    // 12 * num_nodes, [R|t] : 9 R, 3 t
    Eigen::VectorXd InitAffineVector(size_t nodeSize) const;

    //constraints_index is in respect to graph nodes' index
    //TODO: we could probably make a pair<int,double> for index-weight access
    double OptimizeOnce(
        DeformGraph                      &            graph,
        GNParams                         &           params,
        std::vector<std::vector<int>>    &   con_node_neigh,
        std::vector<std::vector<double>> & con_node_weights,
        const NodeCorrs                  &             cons); 

    void DebugRandomVertNodeDeformInfluence(
        std::vector<std::vector<int>>  node_influence_list,
        std::vector<std::vector<double>> vert_node_weights,
        CloudPtr                              moving_cloud,
        pcl::PolygonMesh                        graph_mesh,
        std::string                               out_path);

};

#endif
////////////////////////////////////////////////////////////////////////////////
