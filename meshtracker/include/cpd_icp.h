#ifndef CPD_ICP_H
#define CPD_ICP_H

using CloudPtr = pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr;
using MatCloudPtr = pcl::PointCloud<Eigen::MatrixXf>::Ptr;

struct CPD_Result{

  CloudPtr        aligned_cloud;
  pcl::PolygonMesh aligned_mesh;
  CloudPtr    subsampled_moving;
  CloudPtr     subsampled_fixed;
  CloudPtr       lo_res_aligned;

};

class CPDICP {

public:
    CPDICP(){};
    ~CPDICP(){};

    // conversion functions
    cpd::Matrix PCLToCpdMat(const CloudPtr & ptcloud);
    CloudPtr CpdMatToPCL(const cpd::Matrix & mat);

    cpd::NonrigidResult NonRigidRegistrationImpl(
      const cpd::Matrix &      fixed,
      const cpd::Matrix &     moving,
      const double      & tolerance);

    // encapsulated function call to cpd::getAffinity
    void getAffinity(
      const cpd::Matrix &  fixed, 
      const cpd::Matrix & moving,
      cpd::Matrix       &   _aff)
    {
        cpd::Nonrigid nr;
        nr.getAffinity(fixed,moving,_aff);
    }

    // Use cpd to align high res meshses by downsampling and applying low-res
    // alignment to a high res mesh.
    pcl::PolygonMesh AlignHighResMeshes(
      const pcl::PolygonMesh &  _mesh_fixed,
      const pcl::PolygonMesh & _mesh_moving,
      const double           &    tolerance);

    // Overrride for CloudPtr input
    CPD_Result AlignHighResClouds(
      const CloudPtr &  _cloud_fixed,
      const CloudPtr & _cloud_moving,
      const double   &    tolerance);

    inline void SetTrackedMeshNum(const int & mesh_num){
      tracked_mesh_num=mesh_num;
    }

    inline int GetTrackedMeshNum() const{
      return tracked_mesh_num;
    }



private:
    int       tracked_mesh_num   =   0;
    const int kSamplesUpperLimit = 800;
    const int kSamplesLowerLimit =  50;

};

#endif
