#ifndef _MESH_PROCESSING_H
#define _MESH_PROCESSING_H

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <cstdio>
#include <map>

#include <pcl/io/ply_io.h>
#include <pcl/io/obj_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/point_types.h>
#include <pcl/surface/mls.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/PolygonMesh.h>
#include <pcl/geometry/polygon_mesh.h>
#include <pcl/surface/poisson.h>
#include <pcl/surface/vtk_smoothing/vtk_mesh_smoothing_laplacian.h>
#include <pcl/surface/vtk_smoothing/vtk_mesh_quadric_decimation.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/IO/Ostream_iterator.h>

// #include "cloud_processing.h"
# include "constants.h"
#include "utils.h"
#include "log.h"

typedef pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr                   CloudPtr;
typedef pcl::geometry::PolygonMesh<
    pcl::geometry::DefaultMeshTraits<pcl::PointXYZRGBNormal>>              Mesh;

typedef CGAL::Exact_predicates_inexact_constructions_kernel                   K;
typedef K::Point_3                                                        Point;
typedef K::Vector_3                                                      Vector;
typedef CGAL::Surface_mesh<Point>                                  Surface_mesh;
typedef CGAL::Surface_mesh<K::Point_3>                           Surface_Mesh_3;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor  vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor      face_descriptor;
typedef boost::graph_traits<
    Surface_Mesh_3>::vertex_descriptor                    sm3_vertex_descriptor;
typedef boost::graph_traits<
    Surface_Mesh_3>::face_descriptor                        sm3_face_descriptor;
typedef boost::graph_traits<
    Surface_Mesh_3>::halfedge_descriptor                sm3_halfedge_descriptor;
typedef boost::graph_traits<
    Surface_Mesh_3>::edge_descriptor                        sm3_edge_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

struct halfedge2edge
{
  const Surface_Mesh_3             &  m_mesh;
  std::vector<sm3_edge_descriptor> & m_edges;

  halfedge2edge(
    const Surface_Mesh_3             &     m, 
    std::vector<sm3_edge_descriptor> & edges): 
    m_mesh(m), 
    m_edges(edges)
  {}

  void operator()(
    const sm3_halfedge_descriptor & h) const
  {
    m_edges.push_back(edge(h, m_mesh));
  }
};


namespace mesh_processing
{

struct BBox
{

    double min_x     =  DBL_MAX;
    double min_y     =  DBL_MAX;
    double min_z     =  DBL_MAX;
    double max_x     = -DBL_MAX;
    double max_y     = -DBL_MAX;
    double max_z     = -DBL_MAX;
    double cx        =        0;
    double cy        =        0;
    double cz        =        0;
    double l_diag_sq =        0;

    BBox(const pcl::PolygonMesh & mesh){ this->get(mesh);};
    BBox(const CloudPtr         & cloud){ this->get(cloud);};

    void get(const pcl::PolygonMesh & mesh)
    {
        CloudPtr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
        pcl::fromPCLPointCloud2(mesh.cloud,*cloud);

        for (size_t pt = 0; pt < cloud->points.size(); pt++)
        {
            pcl::PointXYZRGBNormal p = cloud->points[pt];

            if (p.x < this->min_x) this->min_x = p.x; 
            if (p.x > this->max_x) this->max_x = p.x; 
            if (p.y < this->min_y) this->min_y = p.y; 
            if (p.y > this->max_y) this->max_y = p.y; 
            if (p.z < this->min_z) this->min_z = p.z; 
            if (p.z > this->max_z) this->max_z = p.z; 
        }

        this->cx = this->min_x + (this->max_x - this->min_x)*0.5;
        this->cy = this->min_y + (this->max_y - this->min_y)*0.5;
        this->cz = this->min_z + (this->max_z - this->min_z)*0.5;

        this->l_diag_sq = (max_x - min_x) + (max_y - min_y);
    };

    void get(const CloudPtr & cloud)
    {
        for (size_t pt = 0; pt < cloud->points.size(); pt++)
        {
            pcl::PointXYZRGBNormal p = cloud->points[pt];

            if (p.x < this->min_x) this->min_x = p.x; 
            if (p.x > this->max_x) this->max_x = p.x; 
            if (p.y < this->min_y) this->min_y = p.y; 
            if (p.y > this->max_y) this->max_y = p.y; 
            if (p.z < this->min_z) this->min_z = p.z; 
            if (p.z > this->max_z) this->max_z = p.z; 
        }

        this->cx = this->min_x + (this->max_x - this->min_x)*0.5;
        this->cy = this->min_y + (this->max_y - this->min_y)*0.5;
        this->cz = this->min_z + (this->max_z - this->min_z)*0.5;

        this->l_diag_sq = (max_x - min_x) + (max_y - min_y);
    };

    void print()
    {
        std::cout<<"Bbox:"<<std::endl;
        std::cout<<"x: ["<<this->min_x<<", "<<this->max_x<<"] l:"<<this->max_x - this->min_x<<std::endl;
        std::cout<<"y: ["<<this->min_y<<", "<<this->max_y<<"] l:"<<this->max_y - this->min_y<<std::endl;
        std::cout<<"z: ["<<this->min_z<<", "<<this->max_z<<"] l:"<<this->max_z - this->min_z<<std::endl;
        std::cout<<"c: ["<<this->cx<<", "<<this->cy<<", "<<this->cz<<"]"<<std::endl;
    };
};

}

class MeshProcessing
{
private:
    // specify decimation method for abstraction layer creation,
    // options are Quadratic Edge Collapse and Voxel Grid sampling
    enum DecimMethod{QEC, VG};

public:
    MeshProcessing(){};
    ~MeshProcessing(){};

    // Returns an highly-smoothed and simplified version of the input mesh 
    // to be used for abstraction-level processes. Small components are removed
    // and surface manifoldness is ensured. Normals are also recalculated.
    pcl::PolygonMesh CreateSmoothedAbstractionLayer(
        const pcl::PolygonMesh &         in_mesh,
        const int              & target_vertices = 4500,
        const int                     decim_type = QEC);

    // Calculates abstraction layer using poisson reconstruction,
    // Quadradtic Edge Collapse decimation and Laplacian smoothing
    pcl::PolygonMesh PoissonQECAbstraction(
        const pcl::PolygonMesh &        in_mesh,
        const int              & target_vertices
    );

    // Calculates abstraction layer using voxel grid resampling
    // and Laplacian smoothing
    pcl::PolygonMesh VGAbstraction(
        const pcl::PolygonMesh &        in_mesh,
        const int              & target_vertices
    );

    pcl::PolygonMesh SmoothMesh(
        const pcl::PolygonMesh & in_mesh,
        const int              &   iters =    100,
        const bool edge_smoothing        =   true,
        const float feat_angle           = M_PI/2);

    pcl::PolygonMesh DecimateMesh(
        const pcl::PolygonMesh &      in_mesh,
        const int              & target_verts);

    pcl::PolygonMesh PoissonRemesh(
        const pcl::PolygonMesh &   in_mesh,
        const int              & oct_depth);

    // Perform CGAL Isotropic remesh. If not provided, target edge lenth
    // will be calculated from the average edge length of the input mesh.
    pcl::PolygonMesh IsotropicRemesh(
        const pcl::PolygonMesh & in_mesh,
        double        target_edge_length = 0);

    // Calls InstantMeshes to remesh the input. Non manifold faces and small,
    // non-connected components are removed
    pcl::PolygonMesh RemeshInstantMeshes(
        const pcl::PolygonMesh &         in_mesh,
        const size_t           & target_vertices
    );

    // Apply poisson meshing to a point cloud
    // pcl::PolygonMesh PoissonMeshing(
    //     const CloudPtr& in_cloud,
    //     const int& oct_depth);

    // Removes any unconnected components of size less than the input 
    // percent with respect to the global mesh
    pcl::PolygonMesh RemoveSmallComponents(
        const pcl::PolygonMesh & in_mesh,
        const double           & percent = 0.1);

    // Return the average edge length across all faces in input mesh
    double CalculateAverageEdgeLength(
        const pcl::PolygonMesh & mesh);

    // Proportionally scales the bounding box to input height
    // with option to move bottom of mesh to world origin
    // Returns xyz tform
    // TODO: Specify hieight axis, for now assumes longest axis
    // should also have an option to normalize
    double GetAffineScaleFactor(
        const pcl::PolygonMesh &          in_mesh, 
        const double           &    target_height,
        const bool             & translate2origin = false);

    // scale input mesh by given factor 
    pcl::PolygonMesh ApplyScaleFactor(
    const pcl::PolygonMesh     & in_mesh, 
    const double               &      sf,
    const std::array<double,3> &  offset);

    // Fix any non-manifold faces in polygon mesh
    pcl::PolygonMesh FixManifoldness(
        const pcl::PolygonMesh & mesh);

    // Calculates the centroid of all 1-ring neighbours for a given point
    pcl::PointXYZRGBNormal GetOneRingCentroid(
        const pcl::PolygonMesh &    mesh,
        const int              & init_pt);

    // Returns list of indices for 1-ring connected vertices for a 
    // given vertex index
    std::vector<int> GetOneRingNeighbors(
        const pcl::PolygonMesh &    mesh,
        const int              & init_pt);

    // Returns list of indices for 1-ring connected vertices for a 
    // given vertex index
    std::vector<int> GetTwoRingNeighbors(
        const pcl::PolygonMesh &    mesh,
        const int              & init_pt);

    // Returns list of all k-ring connected vertices
    // WARNING: Scales terribly
    std::vector<int> GetKRingNeighbors(
        const pcl::PolygonMesh & mesh,
        const int              &  vid,
        const int              &    k);

    // Checks the first point to see if xyz normals are non-zero
    bool NormalsExist(
        const pcl::PolygonMesh & mesh);

    // Takes a polygon mesh object and calculates vertex normals, 
    // option to return normals separately
    pcl::PolygonMesh CalculateVertexNormals(
        const pcl::PolygonMesh & mesh);

    // Takes a PolygonMesh and returns a point cloud with normals
    CloudPtr GetCloudFromPolygonMesh(
        pcl::PolygonMesh & mesh);

    // Check for NaN points in vertex data
    bool HasNanData(
        const pcl::PolygonMesh & mesh);

    // Performs LERP between source and target mesh. Meshes must have the same
    // number of vertices and topology. Step is expressed as between [0,1]
    pcl::PolygonMesh LERPMesh(
        const pcl::PolygonMesh source,
        const pcl::PolygonMesh target,
        const float              step = 0.5);

    // Fix any nan points by using connectivity info to replace nan data
    // with the centroid of connecting points
    pcl::PolygonMesh FixNanPoints(
        const pcl::PolygonMesh & mesh_in);

    // // Given a mesh, a source and a vector of points, returns the input
    // // vector sorted by descrete geodesic distance wrt input point
    // std::vector<int> SortByGeodesicDistance(
    //     const pcl::PolygonMesh& mesh,
    //     const int& src_pt_idx,
    //     const std::vector<int>& pts_indxs);

    // // Get descrete geodesic distance between src and dst vertex along mesh
    // int GetGeodesicDistance(
    //     const pcl::PolygonMesh& mesh,
    //     const int& src_pt,
    //     const int& dst_pt);
};

namespace MeshProcDefs
{
    const int PQEC = 0;
    const int VG   = 1;
}




#endif