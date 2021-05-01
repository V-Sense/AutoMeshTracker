#include "cloud_processing.h"

////////////////////////////////////////////////////////////////////////////////

CloudPtr CloudProcessing::SmoothCloudMLS(
    const CloudPtr &cloud)
{

  CloudPtr out_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>(*cloud));
  pcl::PointCloud<pcl::PointXYZ>::Ptr xyz_cloud(new pcl::PointCloud<pcl::PointXYZ>);

  pcl::copyPointCloud(*out_cloud, *xyz_cloud);

  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);

  // Output has the PointNormal type in order to store the normals calculated by MLS
  pcl::PointCloud<pcl::PointNormal> mls_points;

  // Init object (second point type is for the normals, even if unused)
  pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;

  mls.setComputeNormals(true);

  // Set parameters
  mls.setInputCloud(xyz_cloud);
  mls.setPolynomialFit(true);
  mls.setSearchMethod(tree);
  mls.setSearchRadius(10); //TODO Define radius

  // Reconstruct
  mls.process(mls_points);

  pcl::copyPointCloud(*xyz_cloud, *out_cloud);

  return out_cloud;
}

////////////////////////////////////////////////////////////////////////////////

double CloudProcessing::ComputeCloudResolution(
    const CloudPtr &cloud)
{

  double res = 0.0;
  int n_points = 0;
  int nres;
  std::vector<int> indices(2);
  std::vector<float> sqr_distances(2);
  pcl::search::KdTree<pcl::PointXYZRGBNormal> tree;

  tree.setInputCloud(cloud);

  for (size_t i = 1; i < cloud->size(); ++i)
  {
    if (!pcl_isfinite((*cloud)[i].x))
    {
      continue;
    }
    //Considering the second neighbor since the first is the point itself.
    nres = tree.nearestKSearch(i, 2, indices, sqr_distances);
    if (nres == 2)
    {
      res += sqrt(sqr_distances[1]);
      ++n_points;
    }
  }

  if (n_points != 0)
  {
    res /= n_points;
  }

  return res;
}

////////////////////////////////////////////////////////////////////////////////

bool CloudProcessing::NormalsExist(
    const CloudPtr &cloud)
{

  if (cloud->points[0].normal_x == 0 &&
      cloud->points[0].normal_y == 0 &&
      cloud->points[0].normal_z == 0)
  {
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////

void CloudProcessing::CalculatePointCloudNormals(
    CloudPtr& cloud)
{
  pcl::NormalEstimationOMP<pcl::PointXYZRGBNormal,pcl::PointXYZRGBNormal> ne;
  ne.setInputCloud (cloud);

  pcl::search::KdTree<pcl::PointXYZRGBNormal>::Ptr tree 
    (new pcl::search::KdTree<pcl::PointXYZRGBNormal> ());
  ne.setSearchMethod (tree);

  double res = ComputeCloudResolution(cloud);
  ne.setRadiusSearch (res * 3);

  ne.compute (*cloud);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> CloudProcessing::FindOutliers(
    const CloudPtr &cloud)
{

  // Radius Outlier removal removes all points that have less than the specified
  // number of neighbours within a given radius
  pcl::RadiusOutlierRemoval<pcl::PointXYZRGBNormal> outrem(true);

  // build the filter
  outrem.setInputCloud(cloud);

  // use 2x density to specify radius
  double radius = ComputeCloudResolution(cloud);
  outrem.setRadiusSearch(radius * 4);
  outrem.setMinNeighborsInRadius(1);
  outrem.setNegative(true); // set true so that output returns the removed points

  // apply filter
  std::vector<int> filtered_indices;
  outrem.filter(filtered_indices);

  return filtered_indices;
}

////////////////////////////////////////////////////////////////////////////////

double CloudProcessing::ComputeHausdorffDistance(
    const CloudPtr& cloud_a,
     const CloudPtr& cloud_b)
{

  // compare A to B
  pcl::search::KdTree<pcl::PointXYZRGBNormal> tree_b;
  tree_b.setInputCloud(cloud_b);
  double max_dist_a = -std::numeric_limits<double>::max ();

  for (const auto &point : cloud_a->points)
  {
    std::vector<int> indices (1);
    std::vector<float> sqr_distances (1);

    tree_b.nearestKSearch (point, 1, indices, sqr_distances);
    if (sqr_distances[0] > max_dist_a)
      max_dist_a = sqr_distances[0];
  }

  // compare B to A
  pcl::search::KdTree<pcl::PointXYZRGBNormal> tree_a;
  tree_a.setInputCloud(cloud_a);
  double max_dist_b = -std::numeric_limits<double>::max ();

  for (const auto &point : cloud_b->points)
  {
    std::vector<int> indices (1);
    std::vector<float> sqr_distances (1);

    tree_a.nearestKSearch (point, 1, indices, sqr_distances);
    if (sqr_distances[0] > max_dist_b)
      max_dist_b = sqr_distances[0];
  }

  double max_dist = max_dist_a > max_dist_b ? max_dist_a : max_dist_b;
  return max_dist;
}

////////////////////////////////////////////////////////////////////////////////

double CloudProcessing::ComputeHausdorffDistance(
  const pcl::PointXYZRGBNormal& pt,
  const CloudPtr& cloud)
{

  pcl::search::KdTree<pcl::PointXYZRGBNormal> tree;
  tree.setInputCloud(cloud);

  std::vector<int> indices (1);
  std::vector<float> sqr_distances (1);

  tree.nearestKSearch(pt, 1, indices, sqr_distances);

  return std::sqrt(sqr_distances[0]);
}

////////////////////////////////////////////////////////////////////////////////
// Uses a voxel grid to approximate a subsampled surface. Can give errors for
// thin objects when the voxel size is larger than the object as the approximation
// will be a point lying in between

CloudPtr CloudProcessing::ResamplePointCloud(CloudPtr& cloud){

  CloudProcessing cp;

  // Create the filtering object
  CloudPtr out_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>(*cloud));
  pcl::VoxelGrid<pcl::PointXYZRGBNormal> sor;
  sor.setInputCloud (out_cloud);

  // set voxel grid size as small as possible. 1.5 has a high chance of being
  // large enough without badly approximating thin surfaces
  double leaf_size = cp.ComputeCloudResolution(cloud);
  sor.setLeafSize (leaf_size, leaf_size, leaf_size);
  sor.filter (*out_cloud);

  return out_cloud;
}


////////////////////////////////////////////////////////////////////////////////

CloudPtr CloudProcessing::UniformsamplePointCloud(
  const CloudPtr& cloud,
  std::optional<int> max_size)
{
  if (!max_size.has_value())
  {
    max_size = this->kSamplesUpperLimit;
  }

  CloudPtr res_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>(*cloud));

  // Create the filtering object
  pcl::VoxelGrid<pcl::PointXYZRGBNormal> sor;
  sor.setInputCloud(res_cloud);

  int cloud_size = cloud->size();
  double reduction_factor = 0.005;

  CloudPtr out_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  // iteratively approach the ideal cloud size
  while (cloud_size > max_size.value()) {

      // get downsampled cloud 
      // (this is quite uniform but it creates new points instead of actually sampling)
      double leaf_size = (reduction_factor);
      sor.setLeafSize (leaf_size, leaf_size, leaf_size);
      sor.filter (*res_cloud);

      // perform knn search for each resampled point toward input cloud
      pcl::KdTreeFLANN<pcl::PointXYZRGBNormal> kdtree;
      kdtree.setInputCloud(cloud);

      const int search_length=2; // just has to be > 1
      std::vector<int> visited_indices;// store visited points to avoid returning duplicates
      for (size_t point = 0; point < res_cloud->size(); point++) {

        std::vector<int> pointIdxNKNSearch(search_length);
        std::vector<float> pointNKNSquaredDistance(search_length);

        pcl::PointXYZRGBNormal search_pt = res_cloud->points[point];

        int nres = kdtree.nearestKSearch(search_pt, search_length,
           pointIdxNKNSearch, pointNKNSquaredDistance);

        if (std::count(visited_indices.begin(),visited_indices.end(),pointIdxNKNSearch[0])==0) {
          out_cloud->push_back(cloud->points[pointIdxNKNSearch[0]]);
          visited_indices.push_back(pointIdxNKNSearch[0]);
        }
      }

      cloud_size = out_cloud->size();
      if (cloud_size > this->kSamplesUpperLimit) {
        out_cloud->clear();
        reduction_factor+=0.005;
        // BOOST_LOG_TRIVIAL(debug) 
        // << "In CloudProcessing::UniformsamplePointCloud: increasing reduction factor to "
        // <<reduction_factor;
      }
  }

  return out_cloud;
}

////////////////////////////////////////////////////////////////////////////////

CloudPtr CloudProcessing::RandomsamplePointCloud(CloudPtr& cloud,
  const unsigned int& samples){

  // Create the filtering object
  CloudPtr out_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>(*cloud));
  pcl::RandomSample<pcl::PointXYZRGBNormal> sor;
  sor.setInputCloud(out_cloud);
  sor.setSample(samples);
  sor.setSeed(rand());
  sor.filter (*out_cloud);

  return out_cloud;
}

////////////////////////////////////////////////////////////////////////////////

CloudPtr CloudProcessing::AlignCentres(
  const CloudPtr& moving,
  const CloudPtr& fixed)
{

    CloudPtr out_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>(*moving));

    cloud_processing::BBox moving_bb(moving);
    cloud_processing::BBox fixed_bb(fixed);

    Eigen::Vector3f tform(
      fixed_bb.cx - moving_bb.cx,
      fixed_bb.cy - moving_bb.cy,
      fixed_bb.cz - moving_bb.cz);

    for (size_t pt = 0; pt < out_cloud->size(); pt++)
    {
      out_cloud->points[pt].x += tform[0];
      out_cloud->points[pt].y += tform[1];
      out_cloud->points[pt].z += tform[2];
    }
        
    return out_cloud;
}

////////////////////////////////////////////////////////////////////////////////

CloudPtr CloudProcessing::TransformToOrigin(
  const CloudPtr& cloud)
{
  CloudPtr out(new pcl::PointCloud<pcl::PointXYZRGBNormal>(*cloud));

  cloud_processing::BBox bb(out);
  double tform[3] = {-bb.cx,-bb.cy,-bb.cz};

  for (size_t pt = 0; pt < out->size(); pt++)
  {
    out->points[pt].x += tform[0];
    out->points[pt].y += tform[1];
    out->points[pt].z += tform[2];
  }

  return out;
}

////////////////////////////////////////////////////////////////////////////////

CloudPtr CloudProcessing::RigidICP(
  const CloudPtr& moving,
  const CloudPtr& fixed,
  const double& max_pair_dist,
  const int& max_iters,
  const double& tform_epsilon,
  const double& fitness_epsilon)
{
  pcl::IterativeClosestPoint<pcl::PointXYZRGBNormal, pcl::PointXYZRGBNormal> icp;
  
  // Set the input source and target
  icp.setInputSource(moving);
  icp.setInputTarget(fixed);

  // Set the max correspondence distance (e.g., correspondences with higher distances will be ignored)
  icp.setMaxCorrespondenceDistance(max_pair_dist); 

  // Set the maximum number of iterations
  icp.setMaximumIterations(max_iters);

  // Set the transformation epsilon
  icp.setTransformationEpsilon(tform_epsilon);
  
  // Set the euclidean distance difference epsilon
  icp.setEuclideanFitnessEpsilon(fitness_epsilon);

  // Perform the alignment
  CloudPtr aligned_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  icp.align(*aligned_cloud);

  return aligned_cloud;
  
}

////////////////////////////////////////////////////////////////////////////////

void CloudProcessing::ColorizeCloud(
  CloudPtr& cloud, 
  std::array<int,3>& rgb)
{

  if (rgb.size() != 3 || rgb[0]<0 || rgb[1]<0 || rgb[2]<0)
  {
    BOOST_LOG_TRIVIAL(warning) 
      << "In ColorizeCloud: Please supply a valid color";
    return;
  }
  
  for (size_t pt = 0; pt < cloud->size(); pt++)
  {
    cloud->points[pt].r = uint8_t(rgb[0]);
    cloud->points[pt].g = uint8_t(rgb[1]);
    cloud->points[pt].b = uint8_t(rgb[2]);
  }
}

////////////////////////////////////////////////////////////////////////////////