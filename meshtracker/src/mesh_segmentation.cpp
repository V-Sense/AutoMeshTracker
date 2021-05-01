#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/IO/PLY_reader.h>
#include <CGAL/IO/PLY_writer.h>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/PolygonMesh.h>
#include <pcl/geometry/polygon_mesh.h>
#include <pcl/geometry/mesh_conversion.h>

#include <boost/foreach.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "cloud_processing.h"
#include "mesh_segmentation.h"
#include "mesh_tracking.h"
#include "matching.h"
#include "utils.h"
#include "log.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SM;

typedef boost::graph_traits<SM>::face_descriptor face_descriptor;
typedef boost::graph_traits<SM>::vertex_descriptor vertex_descriptor;

typedef CGAL::Face_filtered_graph<SM> Filtered_graph;
typedef SM::Property_map<face_descriptor, double> Facet_double_map;
typedef SM::Property_map<face_descriptor, std::size_t> Facet_int_map;
typedef SM::Property_map<SM::Vertex_index, std::size_t> Vertex_id_map;

typedef pcl::geometry::PolygonMesh<pcl::geometry::DefaultMeshTraits<pcl::PointXYZRGBNormal>> Mesh;

////////////////////////////////////////////////////////////////////////////////

// Requires OFF file input
void MeshSegmentation::Segment(const pcl::PolygonMesh& pcl_mesh, VidSidMap& output)
{

  CloudProcessing cp;
  // CGAL segmentation requires OFF file format, so we convert ply to OFF
  // and let cgal read it from a stringstream
  std::stringstream ss_in = utils::PclMesh2CgalStream(pcl_mesh);
  SM mesh;
  ss_in >> mesh; 

  // Create vertex_id Property_map
  SM::Property_map<vertex_descriptor, std::size_t> v_id;
  bool created;
  boost::tie(v_id, created) = mesh.add_property_map<vertex_descriptor, std::size_t>("v:vid", 0);
  assert(created);

  // Assign each vertex an id
  std::size_t idx = 0;
  boost::graph_traits<SM>::vertex_iterator m_vit, m_vend;
  for (boost::tie(m_vit, m_vend) = vertices(mesh); m_vit != m_vend; ++m_vit)
  {
    v_id[*m_vit] = idx++;
  }

  // create sdf property map
  Facet_double_map sdf_property_map;
  sdf_property_map =
      mesh.add_property_map<face_descriptor, double>("f:sdf").first;
  CGAL::sdf_values(mesh, sdf_property_map);

  // create a property-map for segment-ids
  Facet_int_map segment_property_map =
      mesh.add_property_map<face_descriptor, std::size_t>("f:sid").first;

  // segment the mesh using default parameters for number of levels,
  // and smoothing lambda. Any other scalar values can be used instead of using
  // SDF values computed using the CGAL function
  std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh,
    sdf_property_map, segment_property_map, m_num_clusters, m_lambda);

  Filtered_graph segment_mesh(mesh, 0, segment_property_map);

  // BOOST_LOG_TRIVIAL(debug) << "Mesh divided into " << number_of_segments << " segments";

  output = {{0, 0}};
  SetNumSegments(number_of_segments);

  for (unsigned int id = 0; id < number_of_segments; ++id)
  {
    // select all faces for given seg_id
    segment_mesh.set_selected_faces(id, segment_property_map);

    // iterate about all vertices in given selection
    boost::graph_traits<Filtered_graph>::vertex_iterator vit, vend;
    for (boost::tie(vit, vend) = vertices(segment_mesh); vit != vend; ++vit)
    {

      // get the global index of the current vertex in the list of vertices
      long unsigned int vidx = (int)v_id[*vit];

      // don't overwrite existing values.
      if (output.count(vidx) == 0)
      {
        output[vidx] = id;
      }
    }
  }

  // outlier check
  // Outliers are determined on segment level but indices must be stored
  // on global scale in order to remap new seg_ids
  std::vector<int> outlier_global_indices;

  // Create global cloud object
  CloudPtr mesh_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  pcl::fromPCLPointCloud2(pcl_mesh.cloud, *mesh_cloud);

  for (int seg_id = 0; seg_id < number_of_segments; seg_id++)
  {

    // Get segment indices
    std::vector<size_t> seg_indices;
    for (VidSidMap::const_iterator it = output.begin(); it != output.end(); ++it)
    {
      if (it->second == seg_id)
      {
        seg_indices.push_back(it->first);
      }
    }

    // populate cloud object
    CloudPtr seg_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
    for (size_t idx = 0; idx < seg_indices.size(); idx++)
    {
      seg_cloud->push_back(mesh_cloud->points[seg_indices[idx]]);
    }

    // Get local outlier list
    std::vector<int> seg_outliers = cp.FindOutliers(seg_cloud);

    // Perform kdtree search to find the local outlier within global cloud and use
    // the index of its neighbour to remap the seg_id
    pcl::KdTreeFLANN<pcl::PointXYZRGBNormal> kdtree;
    kdtree.setInputCloud(mesh_cloud);

    for (size_t outlier = 0; outlier < seg_outliers.size(); outlier++)
    {

      // search length should include the point and a couple of neighbours
      int search_length = 3;
      std::vector<int> pointIdxNKNSearch(search_length);
      std::vector<float> pointNKNSquaredDistance(search_length);

      // define search point
      pcl::PointXYZRGBNormal search_pt = seg_cloud->points[seg_outliers[outlier]];

      // perform search
      int nres = kdtree.nearestKSearch(search_pt, search_length,
                                       pointIdxNKNSearch, pointNKNSquaredDistance);

      // assign id of it's nearest neighbour
      output[pointIdxNKNSearch[0]] = output.find(pointIdxNKNSearch[1])->second;

      // add first result to global list
      outlier_global_indices.push_back(pointIdxNKNSearch[0]);
    }

  } // end foreach seg

  // BOOST_LOG_TRIVIAL(info) << outlier_global_indices.size() << " outlier(s) found...";

  // apply a segment hierarchy from most-connected to least
  ReOrderSegMapByConnectivity(output,pcl_mesh,m_num_segments); 

  output = FuseSmallSegments(output, pcl_mesh, m_num_segments);

  // // Finally, apply the connectivity heirarchy
  // BOOST_LOG_TRIVIAL(debug) << "Applying segment heirarchy...";
  // ReOrderSegMapByConnectivity(output,pcl_mesh,m_num_segments);
  // BOOST_LOG_TRIVIAL(debug) << "Segment heirarchy applied.";

  //-Debug----------------------------------------------------------------------
  // std::string ofs = "/Meshtracker/resource/tests/oneshot_test/debug/segs/";
  // ExportSegments(output,pcl_mesh,ofs);
  // std::cerr<<"exported segs"<<'\n';
  // std::cin.get();
  // ---------------------------------------------------------------------------


  return;
}


////////////////////////////////////////////////////////////////////////////////

// Fuse small segments based on threshold set in header
VidSidMap MeshSegmentation::FuseSmallSegments(
  const VidSidMap& vs_map, 
  const pcl::PolygonMesh& mesh,
  const int& number_of_segments)
{
  time_t toc, tic; 

  // copy input
  VidSidMap output = vs_map;
  std::queue<int> segs_to_fuse;

  bool all_segs_pass = false;

  // keep fusing until all segments are above threshold
  while(!all_segs_pass)
  {
    // break this for loop everytime a corrective action is needed
    for (int seg_id = 0; seg_id < this->GetNumSegments(); seg_id++)
    {
      int count = 0;

      // count num of vertices in segment
      for (VidSidMap::const_iterator it = output.begin(); it != output.end(); ++it)
      {
        if (it->second == seg_id)
        {
          count++;
        }
      }

      // remove empty segs
      if (count == 0)
      {
        for (VidSidMap::iterator it = output.begin(); it != output.end(); ++it)
        {
          if (it->second > seg_id)
          {
            it->second--;
          }
        }
        this->SetNumSegments(this->GetNumSegments()-1);
        break;
      }

      // perform fuse if sufficiently small
      if (count < m_seg_size_threshold) 
      {
        FuseWithSmallestConnectedSegment(output, mesh, seg_id);
        this->SetNumSegments(this->GetNumSegments()-1);
        break;
      }

      // if we reach the end then all segs are within acceptable size
      if (seg_id == this->GetNumSegments()-1)
      {
        all_segs_pass = true;
      }
    }
  }

  return output;
}

////////////////////////////////////////////////////////////////////////////////

int MeshSegmentation::FuseWithSmallestConnectedSegment(
  VidSidMap &vs_map_,
  const pcl::PolygonMesh &mesh,
  const size_t &fuse_sid)
{

  // find all connected segments
  std::vector<int> con_list = GetSegmentConnectivity(fuse_sid,vs_map_,mesh);

  if (con_list.empty())
  {
      // TODO: add a better solution! Maybe hold it static and interpolate 
      // between current and last frame?
      BOOST_LOG_TRIVIAL(warning) 
        << "In FuseWithSmallestConnectedSegment: couldn't fuse seg_"<<fuse_sid;
      return 0;
  }
  

  // initialize
  int target_sid = con_list[0];

  // if more than 1 connected neighbor then find the smallest
  if (con_list.size() > 1)
  {
    // find smallest
    size_t min_seg_size = DBL_MAX;
    for (size_t seg_id = 0; seg_id < con_list.size(); seg_id++)
    {
      size_t count = 0;
      
      // count vertices
      for (VidSidMap::const_iterator it = vs_map_.begin(); it != vs_map_.end(); ++it)
      {
        if (it->second == con_list[seg_id])
        {
          count++;
        }
      }

      // set min and target if smaller
      if (count < min_seg_size)
      {
        min_seg_size = count;
        target_sid = con_list[seg_id];
      }
    }
  }
  
  // BOOST_LOG_TRIVIAL(debug) <<
  //   "Fusing seg_"<<std::to_string(fuse_sid)<<" into seg_"<<std::to_string(target_sid);

  if (target_sid == -1)
  {
    std::string e = "\tIn MeshSegmentation::FuseWithNearestSegment: Unable to match ID for target segment";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }

  // TODO: check if tmp var is necessary
  VidSidMap tmp_map = vs_map_;
  for (VidSidMap::iterator it = tmp_map.begin(); it != tmp_map.end(); ++it)
  {
    if (it->second == fuse_sid)
    {
      it->second = target_sid;
    }

    if (it->second > fuse_sid)
    {
      it->second--;
    }
  }
  vs_map_ = tmp_map;

  return target_sid;
}

////////////////////////////////////////////////////////////////////////////////

int MeshSegmentation::FindNeighbourSegment(
  const VidSidMap &vs_map,
  const int &seg_id, 
  const pcl::PolygonMesh &_mesh) const
{

  // convert mesh format
  Mesh mesh;
  pcl::geometry::toHalfEdgeMesh(_mesh, mesh);

  typedef Mesh::VertexIndices VertexIndices;
  VertexIndices vert_idxs;

  int border_idx = -1;

  // find relevant vertices
  for (VidSidMap::const_iterator it = vs_map.begin(); it != vs_map.end(); ++it)
  {

    if (it->second == seg_id)
    {

      Mesh::VertexIndex vi = Mesh::VertexIndex(it->first);
      Mesh::VertexAroundVertexCirculator vav =
          mesh.getVertexAroundVertexCirculator(vi);

      // validity test, expect errors otherwise
      if (!vav.isValid())
        continue;
      const Mesh::VertexAroundVertexCirculator vav_end = vav;

      do
      {
        int vav_v = vav.getTargetIndex().get();
        int test_id = vs_map.find(vav_v)->second;
        if (test_id != seg_id)
        {
          return test_id;
        }
      } while (++vav != vav_end);
    }
  }

  // return -1 if nothing found
  return border_idx;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> MeshSegmentation::GetSegIndices(
  const VidSidMap& map, 
  const int& seg_id
)
{
  std::vector<int> indices;
  for (VidSidMap::const_iterator it = map.begin(); it != map.end(); ++it)
  {
    if (it->second == seg_id)
    {
      indices.push_back(it->first);
    }
  }

  if (indices.empty())
  {
    BOOST_LOG_TRIVIAL(warning) << "In GetSegIndices: None found for seg_"<<seg_id;
  }
  
  return indices;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> MeshSegmentation::GetSegmentConnectivity(
    const int& seg_id,
    const VidSidMap& vs_map,
    const pcl::PolygonMesh& _mesh
) const
{

  // convert mesh to halfedge format
  Mesh mesh;
  pcl::geometry::toHalfEdgeMesh(_mesh, mesh);

  typedef Mesh::VertexIndices VertexIndices;
  VertexIndices vert_idxs;

  int connectivity = 0;
  std::vector<int> connected_segs;

  // find relevant vertices
  for (VidSidMap::const_iterator it = vs_map.begin(); it != vs_map.end(); ++it)
  {
    // for each vertex in seg
    if (it->second == seg_id)
    {

      // create a vertex circulator
      Mesh::VertexIndex vi = Mesh::VertexIndex(it->first);
      Mesh::VertexAroundVertexCirculator vav =
          mesh.getVertexAroundVertexCirculator(vi);

      // validity test, expect errors otherwise
      if (!vav.isValid())
        continue;
      const Mesh::VertexAroundVertexCirculator vav_end = vav;

      // if vert neighbor seg_id is different, store it
      do
      {
        int vav_v = vav.getTargetIndex().get();
        int test_id = vs_map.find(vav_v)->second;
        if (test_id != seg_id)
        {
          connected_segs.push_back(test_id);
        }
      } while (++vav != vav_end);
    }
  }

  if (connected_segs.empty())
  {
    BOOST_LOG_TRIVIAL(warning) 
      << "In GetSegmentConnectivity: no connected components";
    return connected_segs;
  }
  

  // Remove duplicates
  connected_segs.erase(
    utils::remove_duplicates(
      connected_segs.begin(),connected_segs.end()),
    connected_segs.end());

  return connected_segs;
}

////////////////////////////////////////////////////////////////////////////////

/// Re-orders the segmentation id's such that the 0th segment is the most 
/// connected segment w.r.t other segments and the nth is least i.e. torso 
/// is an ideal 0th segment while hands and feet would be nth. 
void MeshSegmentation::ReOrderSegMapByConnectivity(
  VidSidMap& vs_map, 
  const pcl::PolygonMesh& mesh,
  const int& num_segs
)
{
  // find and store connectivity of each segment <conn,seg>
  std::vector<std::pair<int,int>> connectivity;

  for (size_t seg = 0; seg < num_segs; seg++)
  {
    // get connectivity
    std::vector<int> conn = GetSegmentConnectivity(seg,vs_map,mesh);

    connectivity.push_back(std::make_pair(conn.size(),seg));
  }

  // sort based on connectivity (1st elem) in descending order
   std::sort(connectivity.begin(), connectivity.end());
   std::reverse(connectivity.begin(),connectivity.end());

  // Re-map vert ids
  for(VidSidMap::iterator it = vs_map.begin(); it != vs_map.end(); ++it)
  {
  
    // find index of of current seg within sorted connectivity vector
    for (size_t seg_id = 0; seg_id < connectivity.size(); seg_id++)
    {
      if (it->second == connectivity[seg_id].second)
      {
        it->second = seg_id;
        break;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> MeshSegmentation::GetSelfAndConnectedSegIndices(
  const VidSidMap& vs_map,
  const pcl::PolygonMesh& in_mesh,
  const int& seg_id,
  const int& k
)
{

  // convert mesh to halfedge format
  Mesh mesh;
  pcl::geometry::toHalfEdgeMesh(in_mesh, mesh);

  typedef Mesh::VertexIndices VertexIndices;
  VertexIndices vert_idxs;

  std::vector<int> connected_segs;
  std::vector<int> output;

  connected_segs = GetSelfAndConnectedSegIds(vs_map,in_mesh,seg_id);

  int iters = connected_segs.size();
  for (size_t sid = 1; sid < iters; sid++)
  {
    std::vector<int> kring = 
      GetSelfAndConnectedSegIds(vs_map,in_mesh,connected_segs[sid]);

    for (size_t k = 0; k < kring.size(); k++)
    {
      if(std::count (connected_segs.begin(), connected_segs.end(), kring[k]) == 0)
      {
        connected_segs.push_back(kring[k]);
      }
    }
  }
  
  // iterate again to add the connected segments
  for (VidSidMap::const_iterator it = vs_map.begin(); it != vs_map.end(); ++it)
  {    
    // if vertex seg_id matches a connected seg, add it to output indices
    if (std::count (connected_segs.begin(), connected_segs.end(), it->second) > 0)
    {
      output.push_back(it->first);
    }
  }

  return output;
}
////////////////////////////////////////////////////////////////////////////////

std::vector<int> MeshSegmentation::GetSelfAndConnectedSegIds(
  const VidSidMap& vs_map,
  const pcl::PolygonMesh& in_mesh,
  const int& seg_id
)
{
  // convert mesh to halfedge format
  Mesh mesh;
  pcl::geometry::toHalfEdgeMesh(in_mesh, mesh);

  typedef Mesh::VertexIndices VertexIndices;
  VertexIndices vert_idxs;

  int connectivity = 0;
  std::vector<int> connected_segs;
  if (std::count (connected_segs.begin(), connected_segs.end(), seg_id)==0)
  {
    connected_segs.push_back(seg_id);
  }

  // foreach vertex
  for (VidSidMap::const_iterator it = vs_map.begin(); it != vs_map.end(); ++it)
  {
    if (it->second == seg_id)
    {
      Mesh::VertexIndex vi = Mesh::VertexIndex(it->first);
      Mesh::VertexAroundVertexCirculator vav =
          mesh.getVertexAroundVertexCirculator(vi);

      // validity test, expect errors otherwise
      if (!vav.isValid())
        continue;
      const Mesh::VertexAroundVertexCirculator vav_end = vav;

      // circulate current vertex and add seg_id to output for each 
      // encoutered seg_id != current
      do
      {
        int vav_v = vav.getTargetIndex().get();
        int test_id = vs_map.find(vav_v)->second;
        if (test_id != seg_id)
        {
          if (std::count (connected_segs.begin(), connected_segs.end(), test_id)==0)
          {
            connected_segs.push_back(test_id);
          }
          connectivity++;
        }
      } while (++vav != vav_end);
    }
  }
  
  return connected_segs;
}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshSegmentation::MeshSegFromVertexIds(
  const pcl::PolygonMesh& global_mesh,
  const VidSidMap& vs_map,
  const int& seg_id)
{

  // find all vertices in vs_map for given seg_id.
  //  foreach vertex, knn match with vertices in cloud.
  CloudPtr vs_map_pts(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  CloudPtr mesh_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  pcl::fromPCLPointCloud2(global_mesh.cloud, *mesh_cloud);

  for (VidSidMap::const_iterator it = vs_map.begin(); it != vs_map.end(); ++it)
  {
    if (it->second == seg_id)
    {
      vs_map_pts->push_back(mesh_cloud->points[it->first]);
    }
  }

  // BOOST_LOG_TRIVIAL(debug) << "Creating seg mesh, size: "<<vs_map_pts->points.size();
  Matcher mt;
  // std::vector<match_conf> vert_ids_conf;
  // vert_ids_conf = mt.GetMatchesAndConfidence(vs_map_pts, mesh_cloud);
  std::vector<int> vert_ids = mt.GetKNNMatches(vs_map_pts,mesh_cloud);
  std::vector<VertexIndices> pcl_vert_idxs; 

  pcl::PolygonMesh out_mesh;

  // convert he_mesh format
  Mesh he_mesh, he_out_mesh;
  pcl::geometry::toHalfEdgeMesh(global_mesh, he_mesh);
  he_out_mesh = he_mesh; 

  //  use index to point to vertex for fav circulator.
  //  add fav faces to selected faces.
  std::map<Mesh::FaceIndex,bool> map_face_bool_keep;
  for (size_t face = 0; face < he_out_mesh.sizeFaces(); face++)
  {
     map_face_bool_keep.insert(std::make_pair(Mesh::FaceIndex(face),false)); 
  }

  for(size_t vert = 0; vert < vert_ids.size(); vert++)
  {

    // Mesh::VertexIndex vi = Mesh::VertexIndex(vert_ids[vert].first);
    Mesh::VertexIndex vi = Mesh::VertexIndex(vert_ids[vert]);
    
    if (!vi.isValid()){
      continue;
    }

    
    FAVC circ_fav = he_mesh.getFaceAroundVertexCirculator(vi);
    const FAVC circ_fav_end = circ_fav;
    if (int(circ_fav_end.getTargetIndex().get())==0)
    {
      continue;;
    }
    

    do
    {
      // Very important: Some half_edges are on the boundary
      //  -> have an invalid face index
      if (!he_mesh.isBoundary (circ_fav.getCurrentHalfEdgeIndex ()))
      {
        Mesh::FaceIndex face_index = circ_fav.getTargetIndex();
        map_face_bool_keep[face_index] = true;
      }
      else
      {
        continue;
      }
    } while (++circ_fav!=circ_fav_end);
  }

  size_t remaining_faces_count_expected = 0;
  for (size_t face = 0; face < map_face_bool_keep.size(); face++)
  {
    if(map_face_bool_keep[Mesh::FaceIndex(face)]){
      remaining_faces_count_expected++;
    }
  }
  // Remove all non-flagged faces by index
  // Warning: deleteFace() seems only to flag faces as deleted,
  // Call cleanUp() to actually commit changes :)
  for (size_t face = 0; face < he_out_mesh.sizeFaces(); face++)
  {

    if (!map_face_bool_keep[Mesh::FaceIndex(face)]){
      he_out_mesh.deleteFace(Mesh::FaceIndex(face));
      if(!he_out_mesh.isDeleted(Mesh::FaceIndex(face))){
        BOOST_LOG_TRIVIAL(warning) << "In MeshSegFromVertexIds: Face not properly flagged for removal. ";
      }
    }
  }

  he_out_mesh.cleanUp();
  
  // Confirm expected behaviour from cleanUp()
  if (he_out_mesh.sizeFaces() != remaining_faces_count_expected){
    std::string e = "In MeshSegFromVertexIds: Failed to return correct size for input segment.";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }

  pcl::geometry::toFaceVertexMesh(he_out_mesh,out_mesh);

  return out_mesh;
}

////////////////////////////////////////////////////////////////////////////////

std::map<int,int> MeshSegmentation::GetSegSizes(
  const VidSidMap& map, 
  const bool& print = false)
{
  std::map<int,int> counts; 
  for (VidSidMap::const_iterator it = map.begin(); it != map.end(); ++it)
  {
    if(counts.find(it->second)==counts.end())
    {
      counts.insert(std::make_pair(it->second,1));
    }
    else
    {
      counts[it->second]++;
    }
  }
  
  if(print)
  {
    std::cerr<<"seg_id   num_verts"<<'\n';
    for(std::map<int,int>::const_iterator it = counts.begin(); it != counts.end(); ++it)
    {
      std::cerr<<it->first<<" "<<it->second<<'\n';
    }
  }
  return counts;
}

////////////////////////////////////////////////////////////////////////////////

void MeshSegmentation::ExportSegments(
  const VidSidMap &vs_map,
  const pcl::PolygonMesh &mesh,
  const std::string &out_path
) const
{

  for (size_t seg_id = 0; seg_id < this->m_num_segments; seg_id++)
  {

    MeshSegmentation ms;
    CloudProcessing cp;
    pcl::PolygonMesh out_mesh = ms.MeshSegFromVertexIds(mesh,vs_map,seg_id);

    std::array<int,3> rgb = {
    utils::RandInt(0,255),
    utils::RandInt(0,255),
    utils::RandInt(0,255)};

    if (seg_id == -1)
    {
      rgb = {255,255,255};
    }

    CloudPtr tmp_c(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
    pcl::fromPCLPointCloud2(out_mesh.cloud,*tmp_c);
    cp.ColorizeCloud(tmp_c,rgb);
    pcl::toPCLPointCloud2(*tmp_c,out_mesh.cloud);

    std::string filename = out_path + "seg_" + std::to_string(seg_id) + ".ply";
    pcl::io::savePLYFile(filename, out_mesh);
  }
}

////////////////////////////////////////////////////////////////////////////////