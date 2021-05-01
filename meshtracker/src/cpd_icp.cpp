#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/io/ply_io.h>
#include <pcl/search/kdtree.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/random_sample.h>

#include <Eigen/Geometry>

#include "cpd_nonrigid.h"

#include "utils.h"
#include "cloud_processing.h"
#include "cpd_icp.h"

////////////////////////////////////////////////////////////////////////////////

// Utility function for converting PCL point cloud ptrs to eigen matrices for
// use with cpd-ICP library
cpd::Matrix CPDICP::PCLToCpdMat(const CloudPtr &ptcloud){

    cpd::Matrix mat(ptcloud->size(),3);

    for (size_t pt = 0; pt < ptcloud->size(); pt++) {
        mat(pt,0)=ptcloud->points[pt].x;
        mat(pt,1)=ptcloud->points[pt].y;
        mat(pt,2)=ptcloud->points[pt].z;
    }

    return mat;
}

////////////////////////////////////////////////////////////////////////////////

// Utility function for converting eigen matrices back to pcl point clouds
CloudPtr CPDICP::CpdMatToPCL(const cpd::Matrix &mat){

    CloudPtr ptcloud (new pcl::PointCloud<pcl::PointXYZRGBNormal>);

    for (size_t pt = 0; pt < mat.rows(); pt++) {

        pcl::PointXYZRGBNormal pcl_pt;
        pcl_pt.x = (float) mat(pt,0);
        pcl_pt.y = (float) mat(pt,1);
        pcl_pt.z = (float) mat(pt,2);
        ptcloud->push_back(pcl_pt);
    }

    return ptcloud;
}

////////////////////////////////////////////////////////////////////////////////

// Nonrigid ICP registration of two eigen matrix point sets. Returns a results
// object which contains the output point set among other paramters
cpd::NonrigidResult CPDICP::NonRigidRegistrationImpl(
  const cpd::Matrix &fixed,
  const cpd::Matrix &moving,
  const double& tolerance){

    // Matrices must be equal size
    if (moving.rows() != fixed.rows()) {
      std::string e = "Size error In CPDICP::NonRigidRegistrationImpl, input matrices rows must match size\n";
      e += "\tfixed: "+std::to_string(fixed.rows())+","+std::to_string(fixed.cols())+"\n";
      e += "\tmoving: "+std::to_string(moving.rows())+","+std::to_string(moving.cols())+"\n";
      BOOST_LOG_TRIVIAL(error) << e;
      throw std::runtime_error(e);
    }

    cpd::Nonrigid nonrigid;
    nonrigid.tolerance(tolerance); //Essential!

    return nonrigid.run(fixed,moving);
}

////////////////////////////////////////////////////////////////////////////////
// Similar to AlignHigherResClouds except it uses affinity to estimate the 
// relationship between the lo-res tform and hi-res.

pcl::PolygonMesh CPDICP::AlignHighResMeshes(
  const pcl::PolygonMesh& _mesh_fixed,
  const pcl::PolygonMesh& _mesh_moving,
  const double& tolerance){

  // Create PCL Pointers
  CloudPtr cloud_moving (new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  CloudPtr cloud_fixed (new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  // extract cloud data from meshes
  pcl::fromPCLPointCloud2(_mesh_moving.cloud, *cloud_moving);
  pcl::fromPCLPointCloud2(_mesh_fixed.cloud, *cloud_fixed);

  // Copy our source before downsampling
  CloudPtr cloud_moving_hi (new pcl::PointCloud<pcl::PointXYZRGBNormal>(*cloud_moving));
  CloudPtr cloud_moving_lo (new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  // Don't bother sampling if it's a small enough segment
  if (cloud_moving->size()<this->kSamplesUpperLimit) {

    cloud_moving_lo = cloud_moving;

  } else {

    CloudProcessing cp;

    // Approximate uniform sampling
    std::optional<int> target_samples; // uses default if not init
    cloud_moving_lo = cp.UniformsamplePointCloud(
      cloud_moving,target_samples);

    // if uniform sample fails, default to random sample
    if (cloud_moving_lo->size() < this->kSamplesLowerLimit) {
      cloud_moving_lo = cp.RandomsamplePointCloud(cloud_moving,(int)this->kSamplesUpperLimit);
    }
  }
  
  // Create and populate downsampled fixed cloud
  CloudPtr cloud_fixed_lo (new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  if (cloud_fixed->size()<this->kSamplesUpperLimit) {

    cloud_fixed_lo = cloud_fixed;

  } else {

    CloudProcessing cp;

    // Approximate uniform sampling
    std::optional<int> target_samples; 
    cloud_fixed_lo = cp.UniformsamplePointCloud(cloud_fixed,target_samples);

    // if uniform sample fails, default to random sample
    if (cloud_fixed_lo->size() < this->kSamplesLowerLimit) {
      cloud_fixed_lo = cp.RandomsamplePointCloud(cloud_fixed,(int)this->kSamplesUpperLimit);
    }
  }

  // Convert from ptcloud to Matrix
  cpd::Matrix mat_moving_hi = PCLToCpdMat(cloud_moving_hi);
  cpd::Matrix mat_moving_lo = PCLToCpdMat(cloud_moving_lo);
  cpd::Matrix mat_fixed_lo = PCLToCpdMat(cloud_fixed_lo);

  // Matrices must be equal size
  if (mat_moving_lo.rows() > mat_fixed_lo.rows()) {
      mat_moving_lo.resize(mat_moving_lo.rows(),mat_fixed_lo.cols());
  } else {
      mat_fixed_lo.resize(mat_fixed_lo.rows(),mat_moving_lo.cols());
  }

  // perform lo-res alignment
  cpd::NonrigidResult nr_result = NonRigidRegistrationImpl(
    mat_fixed_lo,
    mat_moving_lo,
    tolerance);

  Eigen::MatrixXd weights;
  weights = nr_result.tform;

  Eigen::MatrixXd m_g;

  getAffinity(mat_moving_hi,mat_moving_lo,m_g);

  // apply new transformation
  mat_moving_hi = mat_moving_hi + m_g*weights;

  // Convert new points to cloud
  CloudPtr cloud_out (new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  cloud_out = CpdMatToPCL(mat_moving_hi);

  // apply new points to input mesh
  pcl::PolygonMesh mesh_out;
  pcl::PCLPointCloud2 cloud_aligned;
  pcl::toPCLPointCloud2(*cloud_out, cloud_aligned);
  mesh_out = _mesh_moving;
  mesh_out.cloud = cloud_aligned;

  return mesh_out;
}

////////////////////////////////////////////////////////////////////////////////

CPD_Result CPDICP::AlignHighResClouds(
  const CloudPtr& in_cloud_fixed,
  const CloudPtr& in_cloud_moving,
  const double& tolerance){

  // Create PCL Pointers, hard copy instances
  CloudPtr cloud_moving (new pcl::PointCloud<pcl::PointXYZRGBNormal>(*in_cloud_moving));
  CloudPtr cloud_fixed (new pcl::PointCloud<pcl::PointXYZRGBNormal>(*in_cloud_fixed));

  // Copy our source before downsampling
  CloudPtr cloud_moving_hi (new pcl::PointCloud<pcl::PointXYZRGBNormal>(*cloud_moving));
  CloudPtr cloud_moving_lo (new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  bool was_downsampled = false;
  // Don't bother sampling if it's a small enough segment
  if (cloud_moving->size() <= this->kSamplesUpperLimit) {
    cloud_moving_lo = cloud_moving;

  } else {

    CloudProcessing cp;
    // Approximate uniform sampling
    std::optional<int> target_samples;
    // cloud_moving_lo = cp.UniformsamplePointCloud(
    //   cloud_moving,target_samples);
    cloud_moving_lo = cp.RandomsamplePointCloud(
      cloud_moving,this->kSamplesUpperLimit-10);
    was_downsampled = true;
  }

  // Create and populate downsampled fixed cloud
  CloudPtr cloud_fixed_lo (new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  if (cloud_fixed->size() <= this->kSamplesUpperLimit) {
    cloud_fixed_lo = cloud_fixed;

  } else {
    CloudProcessing cp;
    // Approximate uniform sampling
    std::optional<int> target_samples;
    // cloud_fixed_lo = cp.UniformsamplePointCloud(
    //   cloud_fixed, target_samples);
    cloud_fixed_lo = cp.RandomsamplePointCloud(
      cloud_fixed, this->kSamplesUpperLimit-10);
    was_downsampled = true;
  }

  if (cloud_moving_lo->size() == 0)
  {
    std::string e = "In CPDICP::AlignHighResClouds: failed to resample input moving";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }

  if (cloud_fixed_lo->size() == 0)
  {
    std::string e = "In CPDICP::AlignHighResClouds: failed to resample input fixed";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }
  CPD_Result result;
  result.subsampled_moving = cloud_moving_lo;
  result.subsampled_fixed = cloud_fixed_lo;

  // Convert from ptcloud to Matrix
  cpd::Matrix mat_moving_hi = PCLToCpdMat(cloud_moving);
  cpd::Matrix mat_moving = PCLToCpdMat(cloud_moving_lo);
  cpd::Matrix mat_fixed = PCLToCpdMat(cloud_fixed_lo);

  cpd::Matrix test(mat_moving);

  // Matrices must be equal size, subsampling has a margin of error so we resize to smallest
  if (mat_moving.rows() != mat_fixed.rows()) {
    // std::string w = "\tIn CPDICP::AlignHighResClouds: correcting downsample mismatch";
    // BOOST_LOG_TRIVIAL(debug) << w;
    if (mat_moving.rows() > mat_fixed.rows()) {
      Eigen::MatrixXd nu_moving = mat_moving.topRows(mat_fixed.rows());
      mat_moving = nu_moving;
    } else {
      Eigen::MatrixXd nu_fixed = mat_fixed.topRows(mat_moving.rows());
      mat_fixed = nu_fixed;
    }
  }

  // perform lo-res alignment
  cpd::NonrigidResult nr_result = NonRigidRegistrationImpl(
    mat_fixed,
    mat_moving,
    tolerance);

  // convert back to pcl format
  result.lo_res_aligned = CpdMatToPCL(nr_result.points);

  // Only apply this if the input was downsampled
  if (was_downsampled)
  {
    Eigen::MatrixXd weights;
    weights = nr_result.tform;

    Eigen::MatrixXd m_g;
    getAffinity(mat_moving_hi,test,m_g);

    // apply new transformation
    mat_moving_hi = mat_moving_hi + m_g*weights;

    // Convert new points to cloud
    CloudPtr cloud_out (new pcl::PointCloud<pcl::PointXYZRGBNormal>);
    cloud_out = CpdMatToPCL(mat_moving_hi);
    result.aligned_cloud = cloud_out;
  } else 
  {
    result.aligned_cloud = result.lo_res_aligned;
    result.subsampled_moving = cloud_moving;
  }
  return result;
}

////////////////////////////////////////////////////////////////////////////////
