#ifndef MESH_SEGMENTATION_H
#define MESH_SEGMENTATION_H

// vertex_id , segment_id
typedef std::map<long unsigned int, long unsigned int>                VidSidMap;
typedef CGAL::Exact_predicates_inexact_constructions_kernel              Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3>                                  SM;

class MeshSegmentation
{

private:
  // default paramters seem good for 25K meshes of typical human body poses
  // generally aim to create as many segments as possible and then fuse the small ones

  int m_num_clusters       =    20; // number of segment levels about large sdf values, recc = 20
  double m_lambda          =   0.1; // smoothness param (0.0 - 1.0), recc = 0.1
  size_t m_num_segments    =     0; //init member variable
  int m_seg_size_threshold =   100; // minimum allowed number of vertices for a segment
  bool m_seg_size_passed   = false;

public:
  MeshSegmentation(){};
  ~MeshSegmentation(){};

  // Takes a CGAL::Surface_mesh, performs segmentation based
  // on Shape diameter function. Outputs a vector of vert_id->seg_id pairs
  void Segment(
      const pcl::PolygonMesh &   mesh,
      VidSidMap              & output);

  // Looks for segments smaller than m_seg_size_threshold and fuse with nearest
  // connected neighbour
  VidSidMap FuseSmallSegments(
      const VidSidMap        &             vs_map,
      const pcl::PolygonMesh &               mesh,
      const int              & number_of_segments);

  // Finds connected segments and merges indices with the smallest neighbor. 
  // Returns the index of the fused segment
  int FuseWithSmallestConnectedSegment(
      VidSidMap              & vs_map_,
      const pcl::PolygonMesh &    mesh,
      const size_t           & ctr_sid);

  // Circulates each vertex in mesh to find a connected vertex with different
  // seg_id. Result is first unique segment. Returns -1 if none found
  int FindNeighbourSegment(
      const VidSidMap        & vs_map,
      const int              & seg_id,
      const pcl::PolygonMesh &  _mesh) const;

  inline void SetNumClusters(const int & nc)
  {
    m_num_clusters = nc;
  }

  inline void SetLambda(const double & l)
  {
    m_lambda = l;
  }

  inline size_t GetNumSegments() const
  {
    return m_num_segments;
  }

  inline void SetNumSegments(const size_t & ns)
  {
    m_num_segments = ns;
  }

  inline void SetSizeCheckPass(const bool & state)
  {
    m_seg_size_passed = state;
  }

  // Returns a list of all indices for given segment id
  std::vector<int> GetSegIndices(
      const VidSidMap &    map,
      const int       & seg_id);

  // Returns a vector of connected seg ids
  std::vector<int> GetSegmentConnectivity(
      const int              & seg_id,
      const VidSidMap        & vs_map,
      const pcl::PolygonMesh &  _mesh) const;

  // Reorders the seg_ids such that 0 is the most-connected descending
  void ReOrderSegMapByConnectivity(
      VidSidMap              &   vs_map,
      const pcl::PolygonMesh &     mesh,
      const int              & num_segs);

  // returns vertex indices of self and k-connected segments
  std::vector<int> GetSelfAndConnectedSegIndices(
      const VidSidMap        & vs_map,
      const pcl::PolygonMesh &   mesh,
      const int              & seg_id,
      const int              &      k = 2);

  // returns the seg_ids for input and connected segments
  std::vector<int> GetSelfAndConnectedSegIds(
      const VidSidMap        & vs_map,
      const pcl::PolygonMesh &   mesh,
      const int              & seg_id);

  // Given a VidSidMap vertex seg_id map, returns a polygonmesh for the given
  // seg_id consisting of all the faces around the input vertices. No conisderation
  // is given for overlapping faces.
  // NOTE: Ideally we would just add the required faces by index but for some
  // reason we get a seg fault from PCL when attempting to manually triangulate
  // input vertices. Hence, the implementation below takes the less efficient route
  // of adding every face and then slectively removing faces not associated with
  // the input vertex set.
  pcl::PolygonMesh MeshSegFromVertexIds(
      const pcl::PolygonMesh & global_mesh,
      const VidSidMap        &      vs_map,
      const int              &      seg_id);

  // Prints the distribution of points among discovered segments, can be returned
  std::map<int, int> GetSegSizes(
      const VidSidMap &   map,
      const bool      & print);

  // Save segments to file for debug
  void ExportSegments(
      const VidSidMap        &   vs_map,
      const pcl::PolygonMesh &     mesh,
      const std::string      & out_path) const;
};

#endif
