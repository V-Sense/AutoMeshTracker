#ifndef _TRI_MESH_H_
#define _TRI_MESH_H_

#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/PolygonMesh.h>
#include <pcl/geometry/polygon_mesh.h>
#include <pcl/geometry/mesh_conversion.h>

typedef pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr CloudPtr;
typedef pcl::geometry::PolygonMesh<pcl::geometry::DefaultMeshTraits<pcl::PointXYZRGBNormal>> Mesh;

class TriMesh
{
public:
	struct PolyIndex
	{
		unsigned int vert_index[3];
		unsigned int norm_index[3];
	};
	std::vector<Eigen::Vector3d> vertex_coord;//
	std::vector<Eigen::Vector3d> norm_coord;
	std::vector<PolyIndex> polyIndex;
	std::vector<Eigen::Vector3d> face_norm;
	Eigen::Vector3d BoundingBox_Min, BoundingBox_Max;	//min max
	int vert_num;
	int poly_num;

	TriMesh() : vert_num(0), poly_num(0) {}
	TriMesh(const char* filename);
	TriMesh(pcl::PolygonMesh& pcl_mesh);
	CloudPtr ToPclPointCloud();
	void updateNorm();
	void prepareFaceNeighbours(std::vector<std::vector<int>>& neighbours);
	void saveOBJ(const char* filename);
	void savePLY(const char* filename);
	void getBoundingBox(Eigen::Vector3d& Min, Eigen::Vector3d& Max);
private:
	void readPLY(const char* filename);
	void readOBJ(const char* filename);
};

#endif/*_TRI_MESH_H_*/



