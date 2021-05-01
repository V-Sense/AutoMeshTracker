#ifndef _UTILS_H
#define _UTILS_H

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <cstdio>
#include <random>
#include <chrono>

#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/surface/mls.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/PolygonMesh.h>
#include <pcl/search/kdtree.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/filters/bilateral.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/geometry/polygon_mesh.h>
#include <pcl/geometry/mesh_conversion.h>
#include <pcl/geometry/mesh_circulators.h>

#include "log.h"

typedef std::chrono::high_resolution_clock c_clock;

typedef pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr CloudPtr;
typedef pcl::geometry::PolygonMesh<pcl::geometry::DefaultMeshTraits<pcl::PointXYZRGBNormal>> Mesh;

typedef Mesh::VertexIndex VertexIndex;
typedef Mesh::VertexIndices VertexIndices;
typedef Mesh::FaceAroundVertexCirculator FAVC;
typedef Mesh::VertexAroundFaceCirculator VAFC;

// Vertex_id , Segment_id mapping
typedef std::map<long unsigned int,long unsigned int>                 VidSidMap;


namespace utils
{

static auto square = [](const float argu) { return argu * argu; };
static auto cube = [](const float argu) { return argu * argu * argu; };
static auto max = [](const float lhs, const float rhs) { return lhs > rhs ? lhs : rhs; };

static bool PairSortPredicate(const std::pair<int,int>& lhs, const std::pair<int,int>& rhs)
{
  return (lhs.first == rhs.first && lhs.second == rhs.second);
}

////////////////////////////////////////////////////////////////////////////////

// Return time difference in s using std::chrono::high_resolution_clock points
static double time_elapsed(c_clock::time_point& t1, c_clock::time_point& t2)
{
  return std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
}

////////////////////////////////////////////////////////////////////////////////

// Given PolygonMesh, write to stringstream for CGAL. More versatile to provide
// stream rather than a specific type as cgal can automatically assign type.
static std::stringstream PclMesh2CgalStream(const pcl::PolygonMesh &mesh)
{

  const unsigned int num_points = mesh.cloud.height * mesh.cloud.width;
  const unsigned int num_polygons = mesh.polygons.size();

  if (num_points == 0 || num_polygons == 0)
  {
    std::string e = "In PclMesh2CgalStream: Invalid input mesh";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }
  
  // convert to more manageable format
  CloudPtr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  pcl::fromPCLPointCloud2(mesh.cloud, *cloud);

  std::stringstream ss;

  // write header
  ss << "OFF\n";
  ss << num_points << " " << num_polygons << " 0"
     << "\n";

  // write points
  for (size_t pt = 0; pt < num_points; pt++)
  {
    double x = cloud->points[pt].x;
    double y = cloud->points[pt].y;
    double z = cloud->points[pt].z;
    ss << x << " " << y << " " << z << "\n";
  }

  // write polygons
  for (size_t poly = 0; poly < num_polygons; poly++)
  {
    int idx1 = mesh.polygons[poly].vertices[0];
    int idx2 = mesh.polygons[poly].vertices[1];
    int idx3 = mesh.polygons[poly].vertices[2];
    ss << "3 " << idx1 << " " << idx2 << " " << idx3 << "\n";
  }

  return ss;
}

////////////////////////////////////////////////////////////////////////////////

// Currently only accept triangle polygon objects
// TODO: Update for robustness, headers can vary
static pcl::PolygonMesh CgalStream2PclMesh(std::stringstream &ss)
{

  // create mesh and cloud objects
  pcl::PolygonMesh mesh;
  CloudPtr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  std::string line;
  unsigned int line_num = 0;
  unsigned int num_verts, num_faces;

  while (std::getline(ss, line))
  {

    if (line.empty())
    {
      ++line_num;
      continue;
    }

    // fancy way of reading whitespace delimited data
    std::istringstream iss(line);
    std::vector<std::string> tokens;

    copy(std::istream_iterator<std::string>(iss),
         std::istream_iterator<std::string>(),
         std::back_inserter(tokens));

    // extract header info
    if (line_num == 1)
    {
      num_verts = std::stoi(tokens[0]);
      num_faces = std::stoi(tokens[1]);
      ++line_num;
      continue;
    }

    if (tokens[0].compare("OFF") == 0)
    {
      ++line_num;
      continue;
    }

    // get polygon tri
    if (tokens[0].compare("3") == 0)
    {
      unsigned int v1, v2, v3;
      v1 = std::stoi(tokens[1]);
      v2 = std::stoi(tokens[2]);
      v3 = std::stoi(tokens[3]);

      pcl::Vertices verts;
      verts.vertices.push_back(v1);
      verts.vertices.push_back(v2);
      verts.vertices.push_back(v3);
      mesh.polygons.push_back(verts);
    }
    else
    {

      //get vertex
      pcl::PointXYZRGBNormal pt;
      pt.x = std::stof(tokens[0]);
      pt.y = std::stof(tokens[1]);
      pt.z = std::stof(tokens[2]);

      cloud->push_back(pt);
    }
    ++line_num;
  } // end lines

  // attach cloud to ouput mesh
  pcl::toPCLPointCloud2(*cloud, mesh.cloud);

  return mesh;
}

////////////////////////////////////////////////////////////////////////////////

// Checks if vectors are roughly aligned, 
static bool areAligned(
  Eigen::Vector3d a, 
  Eigen::Vector3d b, 
  double tolerance)
{

  // normalize
  a.normalize();
  b.normalize();
  
  bool ans;
  a.dot(b) > tolerance ? ans = true : ans = false;

  return ans;
}

////////////////////////////////////////////////////////////////////////////////

// Simple Linear interpolation, step normalized [0,1]
static double LERP(
  double min, 
  double max, 
  double step)
{
	return min + (max - min) * step;
}

////////////////////////////////////////////////////////////////////////////////

// double* Mat2Quat(const Eigen::MatrixXd& mat)
// {

//   // TODO: check that matrix is orothogonal and special orthogonal
//   float tr = mat(0,0) + mat(1,1) + mat(2,2);
//   double* q[4];

//   if (tr > 0) { 
//     float S = sqrt(tr+1.0) * 2; // S=4*qw 
//     *q[0] = 0.25 * S;
//     *q[1] = (mat(2,1) - mat(1,2)) / S;
//     *q[2] = (mat(0,2) - mat(2,0)) / S; 
//     *q[3] = (mat(1,0) - mat(0,1)) / S; 
//   } else if ((mat(0,0) > mat(1,1))&(mat(0,0) > mat(2,2))) { 
//     float S = sqrt(1.0 + mat(0,0) - mat(1,1) - mat(2,2)) * 2; // S=4*qx 
//     *q[0] = (mat(2,1) - mat(1,2)) / S;
//     *q[1] = 0.25 * S;
//     *q[2] = (mat(0,1) + mat(1,0)) / S; 
//     *q[3] = (mat(0,2) + mat(2,0)) / S; 
//   } else if (mat(1,1) > mat(2,2)) { 
//     float S = sqrt(1.0 + mat(1,1) - mat(0,0) - mat(2,2)) * 2; // S=4*qy
//     *q[0] = (mat(0,2) - mat(2,0)) / S;
//     *q[1] = (mat(0,1) + mat(1,0)) / S; 
//     *q[2] = 0.25 * S;
//     *q[3] = (mat(1,2) + mat(2,1)) / S; 
//   } else { 
//     float S = sqrt(1.0 + mat(2,2) - mat(0,0) - mat(1,1)) * 2; // S=4*qz
//     *q[0] = (mat(1,0) - mat(0,1)) / S;
//     *q[1] = (mat(0,2) + mat(2,0)) / S;
//     *q[2] = (mat(1,2) + mat(2,1)) / S;
//     *q[3] = 0.25 * S;
//   }
  
//   return *q;
// }

////////////////////////////////////////////////////////////////////////////////

// Eigen::MatrixXd Quat2Mat(const double* q)
// {
//   double sqw = q[0]*q[0];
//   double sqx = q[1]*q[1];
//   double sqy = q[2]*q[2];
//   double sqz = q[3]*q[3];

//   Eigen::MatrixXd mat;
//   // invs (inverse square length) 
//   // is only required if quaternion is not already normalised
//   double invs = 1 / (sqx + sqy + sqz + sqw);
//   mat(0,0) = ( sqx - sqy - sqz + sqw)*invs ; // since sqw + sqx + sqy + sqz =1/invs*invs
//   mat(1,1) = (-sqx + sqy - sqz + sqw)*invs ;
//   mat(2,2) = (-sqx - sqy + sqz + sqw)*invs ;
  
//   double tmp1 = q[1]*q[2];
//   double tmp2 = q[3]*q[0];
//   mat(1,0) = 2.0 * (tmp1 + tmp2)*invs ;
//   mat(0,1) = 2.0 * (tmp1 - tmp2)*invs ;
  
//   tmp1 = q[1]*q[3];
//   tmp2 = q[2]*q[0];
//   mat(2,0) = 2.0 * (tmp1 - tmp2)*invs ;
//   mat(0,2) = 2.0 * (tmp1 + tmp2)*invs ;
//   tmp1 = q[2]*q[3];
//   tmp2 = q[1]*q[0];
//   mat(2,1) = 2.0 * (tmp1 + tmp2)*invs ;
//   mat(1,2) = 2.0 * (tmp1 - tmp2)*invs ;  
  
//   return mat;
// }

////////////////////////////////////////////////////////////////////////////////

// used to sort nested vectors based on first index
// TODO: replace with a lambda
static bool colsort( const std::vector<int>& v1, 
               const std::vector<int>& v2 ) { 
 return v1[0] < v2[0]; 
}

////////////////////////////////////////////////////////////////////////////////

static int RandInt(const int& min = 0, const int& max = 100)
{
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_int_distribution<int>  distr(min, max);
  return distr(generator);
}

////////////////////////////////////////////////////////////////////////////////

// Remove duplicates from unsorted list object
template <typename ForwardIterator>
static ForwardIterator remove_duplicates( 
  ForwardIterator first, 
  ForwardIterator last )
{
    auto new_last = first;

    for ( auto current = first; current != last; ++current )
    {
        if ( std::find( first, new_last, *current ) == new_last )
        {
            if ( new_last != current ) *new_last = *current;
            ++new_last;
        }
    }

    return new_last;
}

////////////////////////////////////////////////////////////////////////////////

// from Evan Teran on SO: 
// http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring 

// trim from start
static inline std::string &ltrim(std::string &s) 
{
        s.erase(
          s.begin(), 
          std::find_if(s.begin(), 
          s.end(), 
          std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) 
{
        s.erase(
          std::find_if(s.rbegin(), 
          s.rend(), 
          std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) 
{
        return ltrim(rtrim(s));
}


////////////////////////////////////////////////////////////////////////////////

} // namespace utils
#endif

