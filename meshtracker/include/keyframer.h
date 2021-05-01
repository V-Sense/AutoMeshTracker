#ifndef _KEYFRAMER_H
#define _KEYFRAMER_H

#include <iostream>
#include <vector>

#include <pcl/io/ply_io.h>
#include <pcl/io/obj_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/PolygonMesh.h>
#include <pcl/geometry/polygon_mesh.h>
#include <pcl/geometry/mesh_conversion.h>

#include "mesh_processing.h"

typedef pcl::geometry::PolygonMesh<pcl::geometry::DefaultMeshTraits<pcl::PointXYZRGBNormal>> Mesh;


class KeyFramer{


public:
  KeyFramer(){}
  ~KeyFramer(){}

  // Takes list of meshes as input, calculates keyframe locations and returns
  // kf indices for input list.
  std::vector<int> GetKeyframeIndices(const std::string& _filename);

  // Read mesh list
  std::vector<std::string> ReadMeshList(const std::string& _filename);

  // Computes the area of a surface mesh
  double ComputeMeshArea(const pcl::PolygonMesh& mesh);

  unsigned int ComputeMeshGenus(const Mesh& mesh);

  unsigned int ComputeMeshGenus(const pcl::PolygonMesh& mesh);

  void GenerateSPHDescriptors(
    const std::string & meshlist,
    const std::string & exe);

  // This function only exists to implement feasibility score of Collet et. al
  void GenerateFeasScore(
    std::vector<std::string> meshlist,
    std::string outpath);

  void GenKeyframesAndRegions(
    const std::vector<std::string> & meshlist);

};
#endif
