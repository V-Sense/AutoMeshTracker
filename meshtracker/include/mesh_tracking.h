#ifndef _MESH_TRACKING_
#define _MESH_TRACKING_

#include <pcl/io/ply_io.h>
#include <pcl/io/obj_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/common/centroid.h>
#include <pcl/PolygonMesh.h>
#include <pcl/features/normal_3d.h>
#include <pcl/geometry/mesh_circulators.h>
#include <pcl/geometry/polygon_mesh.h>
#include <pcl/geometry/mesh_conversion.h>
#include <pcl/surface/vtk_smoothing/vtk_utils.h>
#include <pcl/registration/icp.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>

#include "cpd_nonrigid.h"

#include "mesh_segmentation.h"
#include "mesh_processing.h"
#include "cloud_processing.h"
#include "tri_mesh.h"
#include "deform_graph.h"
#include "cpd_icp.h"
#include "nonrigid_icp.h"
#include "prog_opts.h"
#include "constants.h"
#include "matching.h"
#include "log.h"

typedef pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr                   CloudPtr;
typedef std::pair<int,std::vector<int>>                               GroupIdxs;

typedef pcl::geometry::PolygonMesh<
    pcl::geometry::DefaultMeshTraits<pcl::PointXYZRGBNormal>>              Mesh;

typedef CGAL::Exact_predicates_inexact_constructions_kernel            EPKernel;
typedef CGAL::Surface_mesh<EPKernel::Point_3>                                SM;

// Vertex_id , Segment_id mapping
typedef std::map<long unsigned int,long unsigned int>                 VidSidMap;

// Segment_id, centroid mapping
typedef std::map<std::size_t, pcl::PointXYZRGBNormal>                 CtrSidMap;

// corrspondence pairs
typedef std::vector<
    std::pair<pcl::PointXYZRGBNormal,pcl::PointXYZRGBNormal>>         NodeCorrs;


class MeshTracking{

public:

    MeshTracking(){};
    ~MeshTracking(){};

    // Avoid unless you know what you're doing!
    inline void SetCurrSegs(const int& num_segs)
    {
        this->curr_num_segs_ = num_segs;
    }

    inline void SetInitNumSegs( const int& num_segs)
    {
        this->init_num_segs_ = num_segs;
        this->curr_num_segs_ = num_segs;
    }

    inline void setMeshNames(
        const std::vector<std::string>& _meshNames)
    {
        meshNames_ = _meshNames;
    }
    inline const std::vector<std::string>& getMeshNames() const 
    {
        return meshNames_;
    }

    inline void setClouds(
        const std::vector<CloudPtr>& _clouds)
    {
        clouds_ = _clouds;
    }
    
    inline const std::vector<CloudPtr>& getClouds() const 
    {
        return clouds_;
    }
    
    inline void setKeyFrameIndices(
        const std::vector<unsigned int>& _keyframeIndices)
    {
        keyframeIndices_ = _keyframeIndices;
    }

    inline const std::vector<unsigned int>& getKeyFrameIndices() const 
    {
        return keyframeIndices_;
    }

    inline void IncrementIter()
    { 
        this->curr_iter_++; 
    }

    inline void Reset()
    { 
        this->curr_iter_    =   0; 
        this->curr_energy_  = 0.0;
        this->prev_energy_  = 0.0;
    }

    inline int GetCurrentIter(){ return this->curr_iter_; }

    // Read mesh list
    int ReadMeshList(
        const std::string & _fileName);

    // Read keyframes and regions from file
    int ReadRegions(
        const std::string & file_name);

    // Return a mesh from input filename. Throws error if unable
    pcl::PolygonMesh LoadMeshFromFile(
        const std::string & filename);

    // Generates abs layers from input meshlist and saves them to a subfolder 
    // in the out_path. Also generates a new meshlist in the out_path folder
    void GenAbslayersFromList(
        const std::string & out_path);

    // This approach tracks a keyframe back and forward through a region
    // Boolean flag switches on adaptive tracking usage
    void TrackMeshSequence(
        icp::IcpOptions icp_ops,
        std::string out_path);

    // Performs global tracking using the VTK tracker.
    // Call the optimal_nonrigid VTKMeshTracker for global tracking.
    pcl::PolygonMesh DoSurfaceAlignment(
        const pcl::PolygonMesh  & moving,
        const pcl::PolygonMesh  & fixed);

    // Function used to call segmentation-based tracking functions, perform
    // iterative fusion and seam blending.
    pcl::PolygonMesh DoAdaptiveTracking(
        const pcl::PolygonMesh  & _moving_mesh,
        const pcl::PolygonMesh  &  _fixed_mesh);

    // Returns an aligned deformation graph corresponding to the constraints
    // which aligned moving to fixed
    void GetDenseCorrespondences(
        pcl::PolygonMesh &      moving_mesh,
        pcl::PolygonMesh &       fixed_mesh,
        pcl::PolygonMesh &   aligned_graph_,
        pcl::PolygonMesh &        in_graph_,
        const int        & align_layer_size = 1750); //default 2500

    // Uses a list of dense vertex correspondences to match segmentation
    // information between meshes
    VidSidMap SegmentationFromMatches(
        const pcl::PolygonMesh   &   fixed,
        const VidSidMap       input_vs_map,
        const std::vector<int>   & matches) const;

    // Calculates a map of segment centroids and corresponding ids
    CtrSidMap GetSegmentCentroids(
        const VidSidMap         & vs_map,
        const pcl::PolygonMesh  &   mesh) const;

    // Returns the index of the first segment pair which shows 
    // a disparity equal or greater than member threshold. 
    // Returns -1 if all segments are coherent
    int SegmentProjectionCoherenceTest(
        VidSidMap & moving_vs_map,
        VidSidMap &  fixed_vs_map);

    // Exports segments for a given vs_map. Used for vizualisation and debug.
    // Currently only exports points
    void ExportSegments(
        const VidSidMap        &   vs_map,
        const pcl::PolygonMesh &     mesh,
        const std::string      & out_path) const;

    // Given seg_maps for two meshes, deforms the 
    CloudPtr AlignSegment(
        CloudPtr        &        moving,
        CloudPtr        &         fixed,
        const VidSidMap & moving_vs_map,
        const VidSidMap &  fixed_vs_map,
        const int       &        seg_id);

    // Rigid alignment of two input meshes, sometimes useful for initialization
    pcl::PolygonMesh AlignRigid(
        pcl::PolygonMesh & moving,
        pcl::PolygonMesh &  fixed);

    // Get correspondence constraints for solver
    NodeCorrs GetCorresConstraints(
        const DeformGraph &               graph,
        const CloudPtr    & aligned_graph_cloud,
        const CloudPtr    &   input_graph_cloud,
        const CloudPtr    &         fixed_cloud);

    // Given two polygon meshes, performs surface-based non-rigid alignment
    pcl::PolygonMesh TrackMesh(
        const pcl::PolygonMesh & _moving_mesh,
        const pcl::PolygonMesh &   _fixed_mesh,
        const float            &     _alpha_max = 15.0, 
        const float            &     _alpha_min =  5.0,  
        const int              &      _iter_num =    5,    
        const float            &         _gamma =  1.0) const;

    // Perform actual deformation
    pcl::PolygonMesh DeformMesh(
        const pcl::PolygonMesh              &             in_mesh,
        const pcl::PolygonMesh              &             d_graph,
        const std::vector<std::vector<int>> & node_influence_list);

    // Returns tform between two meshes
    std::vector<std::array<double,6>> GetTform(
        const pcl::PolygonMesh  & source,
        const pcl::PolygonMesh  & target);

    // Apply input tform to input mesh and return result
    pcl::PolygonMesh ApplyTform(
        const pcl::PolygonMesh                  & in_mesh,
        const std::vector<std::array<double,6>> &   tform,
        const bool                              &  invert = false);

    // save the per-vertex transformation between soure and target
    bool StoreTform(
        const pcl::PolygonMesh  & source,
        const pcl::PolygonMesh  & target,
        const std::string       &   path) const;

    // Perform dense wrapping of source to target to simulate details.
    // Source and target should be aligned as closely as possible. 
    pcl::PolygonMesh GetDetailLayer(
        const pcl::PolygonMesh  &         source,
        const pcl::PolygonMesh  &         target,
        const bool              & strict_matches = true);

    // save debug files, default is false
    void SaveDebugOutput(
        bool        & state, 
        std::string &  path)
    {
        this->save_debug_ = state;
        this->debug_path_ = path;
    }
    
private:

    std::vector<std::string> meshNames_;
    std::vector<pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr> clouds_;
    std::vector<unsigned int> keyframeIndices_;
    std::vector<std::pair<unsigned int, unsigned int> > regions_;
    
    // Only counts segment which fail matching coherence i.e, 
    // does not track segs which were fused due to small size
    std::vector<int> fused_segments_;
    
    int curr_num_segs_  =   0;
    int init_num_segs_  =   0;
    int curr_iter_      =   0;
    int max_iters_      =   0;
    double curr_energy_ = 0.0;
    double prev_energy_ = 0.0;

    // TODO: rename to max_num_vert_node_neighbors_ 
    // This is not a target, but more of an upper limit. 
    // The real influence is determined by the the scaler and sampling radius
    // in deform graph members
    int num_vert_node_neighbors_            =    25; //10
    const double default_CPD_tolerance_     =  1e-8;//1e-8;
    double m_segment_size_match_threshold   =  1.00; // % disparity in nummber of points //0.6
    double tracking_error_thresh_           =  1e-5; //1e-5; 
    double GN_solver_relax_thresh_          = 0.65; // 0.05
    double GN_solver_rigidity_w_min         = 100.0; //100.0

    bool save_debug_        = false;
    std::string debug_path_ = "/MeshTracker/resource/tests/oneshot_test/debug/";
};

#endif