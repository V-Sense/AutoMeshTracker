#include "matching.h"

////////////////////////////////////////////////////////////////////////////////

std::vector<int> Matcher::GetKNNMatches(
    const CloudPtr &source, 
    const CloudPtr &query,
    const int& search_length)
{

  pcl::KdTreeFLANN<pcl::PointXYZRGBNormal> kdtree;

  std::vector<int> matches;

  kdtree.setInputCloud(query);

  for (size_t point = 0; point < source->size(); point++)
  {

    std::vector<int> pointIdxNKNSearch(search_length);
    std::vector<float> pointNKNSquaredDistance(search_length);

    pcl::PointXYZRGBNormal search_pt = source->points[point];

    int nres = kdtree.nearestKSearch(search_pt, search_length,
                                     pointIdxNKNSearch, pointNKNSquaredDistance);

    matches.push_back(pointIdxNKNSearch[0]);
  }

  if (matches.size() != source->size())
  {
    std::string e = "\tIn GetVertexMatchesExact: output doesn't match input size.\n";
    throw std::runtime_error(e);
  }

  return matches;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> Matcher::GetKNNMatches(
    const pcl::PolygonMesh& source, 
    const pcl::PolygonMesh& query,
    const int& search_length)
{


  CloudPtr c_source(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  CloudPtr c_query(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  pcl::fromPCLPointCloud2(source.cloud,*c_source);
  pcl::fromPCLPointCloud2(query.cloud,*c_query);

  pcl::KdTreeFLANN<pcl::PointXYZRGBNormal> kdtree;

  std::vector<int> matches;

  kdtree.setInputCloud(c_query);

  for (size_t point = 0; point < c_source->size(); point++)
  {

    std::vector<int> pointIdxNKNSearch(search_length);
    std::vector<float> pointNKNSquaredDistance(search_length);

    pcl::PointXYZRGBNormal search_pt = c_source->points[point];

    int nres = kdtree.nearestKSearch(search_pt, search_length,
                                     pointIdxNKNSearch, pointNKNSquaredDistance);

    matches.push_back(pointIdxNKNSearch[0]);
  }

  if (matches.size() != c_source->size())
  {
    std::string e = "\tIn GetVertexMatchesExact: output doesn't match input size.\n";
    throw std::runtime_error(e);
  }

  return matches;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> Matcher::GetKNNMatches(
    const pcl::PointXYZRGBNormal& source, 
    const CloudPtr& query,
    const int& _nres,
    int search_length)
{

  if (_nres > search_length)
  {
    search_length = _nres;
  }
  

  pcl::KdTreeFLANN<pcl::PointXYZRGBNormal> kdtree;

  std::vector<int> matches;

  kdtree.setInputCloud(query);

  std::vector<int> pointIdxNKNSearch(search_length);
  std::vector<float> pointNKNSquaredDistance(search_length);

  pcl::PointXYZRGBNormal search_pt = source;

  int nres = kdtree.nearestKSearch(search_pt, search_length,
                                    pointIdxNKNSearch, pointNKNSquaredDistance);

  for (size_t res = 0; (res < nres && res < _nres); res++)
  {
    matches.push_back(pointIdxNKNSearch[res]);
  }
  return matches;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> Matcher::GetConstrainedMatches(
  const pcl::PolygonMesh& source,
  const pcl::PolygonMesh& query,
  double align_tolerance)
{

  // non-destructive cloud copy w/ normals
  MeshProcessing mp;
  pcl::PolygonMesh source_(source);
  pcl::PolygonMesh query_(query);
  CloudPtr cloud_s = mp.GetCloudFromPolygonMesh(source_);
  CloudPtr cloud_q = mp.GetCloudFromPolygonMesh(query_); 

  std::vector<int> matches  = GetConstrainedMatches(
    cloud_s,
    cloud_q, 
    align_tolerance); 
  return matches;

}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> Matcher::GetConstrainedMatches(
  const CloudPtr& source,
  const CloudPtr& query,
  double align_tolerance,
  float radius)
{

  // init for size of query and -1 to flag false matches
  std::vector<int> matches(source->size(),-1); 

  // auto-calc flag
  if(radius == 0)
  {
    mesh_processing::BBox bb(query);
    radius = float(0.3*bb.l_diag_sq);
  }

  typedef std::pair<double,int> npair;

  // foreach point in input cloud
  #pragma omp parallel for
  for (size_t s_pt = 0; s_pt < source->size(); s_pt++)
  {

    // sanity check
    if (
      source->points[s_pt].normal_x == 0 && 
      source->points[s_pt].normal_y == 0 &&
      source->points[s_pt].normal_z == 0)
    {
      std::string e = 
        "In Matcher::GetConstrainedMatches, invalid normal for input point set";
      BOOST_LOG_TRIVIAL(error) << e;
      throw std::runtime_error(e);
    }
    

    // init vector pair: dist,idx
    std::vector<npair> dist_vid;

    // brute force, EVERY point in destination cloud
    for (size_t q_pt = 0; q_pt < query->size(); q_pt++)
    {
      double sq_dist = 0.0;
      sq_dist += (source->points[s_pt].x - query->points[q_pt].x)*(source->points[s_pt].x - query->points[q_pt].x);
      sq_dist += (source->points[s_pt].y - query->points[q_pt].y)*(source->points[s_pt].y - query->points[q_pt].y);
      sq_dist += (source->points[s_pt].z - query->points[q_pt].z)*(source->points[s_pt].z - query->points[q_pt].z);
      double dist = sqrt(sq_dist);

      if ( dist <= radius)
      {
        npair dist_pt(dist,q_pt);
        dist_vid.push_back(dist_pt);
      }
    }
    
    // sort by distance
    std::sort(dist_vid.begin(), dist_vid.end(),
    [](const npair& l, const npair& r) {
      if (l.first != r.first){
        return l.first < r.first;}
      return l.second < r.second; //unlikely to need
    });    

    Eigen::Vector3d src_n{
      source->points[s_pt].normal_x,
      source->points[s_pt].normal_y,
      source->points[s_pt].normal_z
    };

    // insert first one that meets match criteria
    for (size_t entry = 0; entry < dist_vid.size(); entry++)
    {
      Eigen::Vector3d query_n{
        query->points[dist_vid[entry].second].normal_x,
        query->points[dist_vid[entry].second].normal_y,
        query->points[dist_vid[entry].second].normal_z
      };

      if (utils::areAligned(src_n,query_n,align_tolerance))
      {
        matches[s_pt] = dist_vid[entry].second;
        break;
      } 
    }  
  }
  
  return matches;

}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> Matcher::GetConstrainedMatches(
  const pcl::PointXYZRGBNormal& source, 
  const CloudPtr& query,
  const double& tolerance,
  int search_length)
{
  std::vector<int> matches;
  mesh_processing::BBox bb(query);
 
  // TODO: add option to use average edge length if we think that the 
  // input with be uniformly distributed
  float radius = float(0.1*bb.l_diag_sq);

  // sanity check 
  Eigen::Vector3d src_n{
    source.normal_x,
    source.normal_y,
    source.normal_z
  };

  if (src_n[0] == 0.0 && src_n[1] == 0.0 && src_n[2] == 0.0 )
  {
    BOOST_LOG_TRIVIAL(error) 
      << "In GetConstrainedMatches(point): No input normal!"
      << "\n\t Returning empty result.";
    return matches;
  }
  

  // init vector pair: dist,idx
  typedef std::pair<double,int> npair;
  std::vector<npair> dist_vid;

  // brute force, EVERY point in destination cloud
  for (size_t q_pt = 0; q_pt < query->size(); q_pt++)
  {

    double sq_dist = 0.0;
    sq_dist += (source.x - query->points[q_pt].x)*(source.x - query->points[q_pt].x);
    sq_dist += (source.y - query->points[q_pt].y)*(source.y - query->points[q_pt].y);
    sq_dist += (source.z - query->points[q_pt].z)*(source.z - query->points[q_pt].z); 
    double dist = sqrt(sq_dist);

    if ( dist <= radius)
    {
      npair dist_pt(dist,q_pt);
      dist_vid.push_back(dist_pt);
    }
  }

  if(dist_vid.size() == 0)
  {
    BOOST_LOG_TRIVIAL(error) 
      << "In GetConstrainedMatches(point):\tNo matches, returning empty";
    return matches;
  }
  
  // sort by distance
  std::sort(dist_vid.begin(), dist_vid.end(),
  [](const npair& l, const npair& r) {
    if (l.first != r.first){
      return l.first < r.first;}
    return l.second < r.second; //unlikely to need
  });        

  // insert any that match alignment criteria up to requested
  for (size_t entry = 0; entry < dist_vid.size(); entry++)
  {
    Eigen::Vector3d query_n{
      query->points[dist_vid[entry].second].normal_x,
      query->points[dist_vid[entry].second].normal_y,
      query->points[dist_vid[entry].second].normal_z
    };

    // if (utils::areAligned(src_n,query_n,tolerance))
    // {
    matches.push_back(dist_vid[entry].second);
    // } 

    if (matches.size() == search_length)
    {
      break;
    }
  }

  if (matches.size() < 1)
  {
    std::stringstream warn;
    warn << "In GetConstrainedMatches(point):\tNo matches found" 
      << "\n\t size: " + std::to_string(matches.size());
    BOOST_LOG_TRIVIAL(warning) << warn.str();
  }

  return matches;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<match_conf> Matcher::GetMatchesAndConfidence(
  CloudPtr& source, 
  CloudPtr& query,
  const float normal_tol,
  const float max_ray_dist) 
{

  CloudProcessing cp;
  // Calculate normals if needed
  if (!cp.NormalsExist(source))
  {
    BOOST_LOG_TRIVIAL(debug) << "Calculating normals for input source.";
    cp.CalculatePointCloudNormals(source);
  }
  if (!cp.NormalsExist(query))
  {
    BOOST_LOG_TRIVIAL(debug) << "Calculating normals for input query.";
    cp.CalculatePointCloudNormals(query);
  }

  const int search_length=query->size();

  pcl::KdTreeFLANN<pcl::PointXYZRGBNormal> kdtree;

  std::vector<match_conf> matches(source->size());

  kdtree.setInputCloud(query);
  int bad_match_count=0;

  // For debug: mark bad points
  CloudPtr bad_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());

  #pragma omp parallel for
  for (size_t point = 0; point < source->size(); point++)
  {

    std::vector<int> pointIdxNKNSearch(search_length);
    std::vector<float> pointNKNSquaredDistance(search_length);

    pcl::PointXYZRGBNormal search_pt = source->points[point];

    int nres = kdtree.nearestKSearch(search_pt, search_length,
       pointIdxNKNSearch, pointNKNSquaredDistance);

    // specify unit length normal vector for search pt
    float mov_x = search_pt.normal_x;
    float mov_y = search_pt.normal_y;
    float mov_z = search_pt.normal_x;

    Eigen::Vector3f mov(mov_x, mov_y, mov_z);
    mov.normalize();

    // Check for normals
    if (search_pt.normal_x == 0 && search_pt.normal_y == 0 && search_pt.normal_z == 0)
    {
      BOOST_LOG_TRIVIAL(warning) << "Warning In GetICPCorrespondences: Input normals (src) should be non-zero!";
    }
    int skipped = 0;

    // Check for nearest point that obeys normal constraint AND isn't already matched
    if (nres > 1)
    {
      for (int res = 0; res < nres; res++)
      {
        
        pcl::PointXYZRGBNormal pt = query->points[pointIdxNKNSearch[res]];

        float test_x = pt.normal_x;
        float test_y = pt.normal_y;
        float test_z = pt.normal_z;

        // Check for normals
        if (test_x == 0 && test_y == 0 && test_z == 0)
        {
          BOOST_LOG_TRIVIAL(warning) << "Warning In GetICPCorrespondences: Input normals (target) must be non-zero!";
        }

        Eigen::Vector3f test(test_x, test_y, test_z);
        test.normalize();

        match_conf m;
        float dot_prod = mov.dot(test);

        // positive dot product should indicate relatively aligned normals
        if (dot_prod > normal_tol)
        {
          m.first = pointIdxNKNSearch[res];
          m.second = dot_prod;
          matches[point] = m;
          break;
        // else perform raycast search
        } 
        else if(res == nres-1)
        {
          ++bad_match_count;
          bad_cloud->push_back(source->points[point]);

          // Can constrain the ray by adding max_ray_dist to the source_pt normal
          Eigen::Vector3f ray(
            mov_x+max_ray_dist,
            mov_y+max_ray_dist,
            mov_z+max_ray_dist);

          // Get ray eqn using source point and normal 

          // For each point in query cloud, get distance to ray

          m.first = pointIdxNKNSearch[0];
          m.second = dot_prod;
          matches[point]=(m);
          break;
        }
      } // end for:results
    } // end if:results
    else 
    {
        match_conf m;
        m.first = pointIdxNKNSearch[0];
        m.second = 0.0f;
        matches[point] = m;
        ++bad_match_count;
        BOOST_LOG_TRIVIAL(warning) << "In GetVertexMatches: No match found for point, no search performed...";
    } // end match w/normal constraint
  }// end for:points

  BOOST_LOG_TRIVIAL(debug) << bad_match_count<<"/"<<source->size()<<" matches without normal constraint";
  return matches;
}

////////////////////////////////////////////////////////////////////////////////

// TODO: Does this need to be static?
static std::vector<int> GetSegMapCorrespondences(CloudPtr& moving, CloudPtr& fixed) 
{

  CloudProcessing cp;
  // Calculate normals if needed
  if (!cp.NormalsExist(moving))
  {
    BOOST_LOG_TRIVIAL(debug) << "Calculating normals for input moving.";
    cp.CalculatePointCloudNormals(moving);
  }
  if (!cp.NormalsExist(fixed))
  {
    BOOST_LOG_TRIVIAL(debug) << "Calculating normals for input fixed.";
    cp.CalculatePointCloudNormals(fixed);
  }

  const int search_length=10;

  pcl::KdTreeFLANN<pcl::PointXYZRGBNormal> kdtree;

  std::vector<int> matches;

  kdtree.setInputCloud(fixed);
  int bad_match_count=0;

  for (size_t point = 0; point < moving->size(); point++) {

    std::vector<int> pointIdxNKNSearch(search_length);
    std::vector<float> pointNKNSquaredDistance(search_length);

    pcl::PointXYZRGBNormal search_pt = moving->points[point];

    int nres = kdtree.nearestKSearch(search_pt, search_length,
       pointIdxNKNSearch, pointNKNSquaredDistance);
 
    // Check for nearest point that obeys normal constraint
    if (nres > 1){

      // specify unit length normal vector
      float mov_x = search_pt.normal_x;
      float mov_y = search_pt.normal_y;
      float mov_z = search_pt.normal_x;

      Eigen::Vector3f mov(mov_x, mov_y, mov_z);
      mov.normalize();

      // Check for normals
      if (search_pt.normal_x == 0 && search_pt.normal_y == 0 && search_pt.normal_z == 0) {
        BOOST_LOG_TRIVIAL(warning) << "Warning In GetICPCorrespondences: Input normals (src) should be non-zero!";
       }

      for (int res = 0; res < nres; res++) {
        pcl::PointXYZRGBNormal pt = fixed->points[pointIdxNKNSearch[res]];

        float test_x = pt.normal_x;
        float test_y = pt.normal_y;
        float test_z = pt.normal_z;

        // Check for normals
        if (test_x == 0 && test_y == 0 && test_z == 0) {
          BOOST_LOG_TRIVIAL(warning) << "Warning In GetICPCorrespondences: Input normals (target) must be non-zero!";
        }

        Eigen::Vector3f test(test_x, test_y, test_z);
        test.normalize();

        // positive dot product should indicate relatively aligned normals
        if (mov.dot(test)>0){
          matches.push_back(pointIdxNKNSearch[res]);
          break;

        // else just push back a flag
        } else if(res == nres-1){
          ++bad_match_count;
          matches.push_back(-1);
          break;
        }
      } // end for

    } else {
      matches.push_back(pointIdxNKNSearch[0]);
      ++bad_match_count;
      BOOST_LOG_TRIVIAL(warning) << "In GetSegMapCorrespondences: No match found for point, no search performed...";
    } // end match w/normal constraint
  }

  BOOST_LOG_TRIVIAL(info) << bad_match_count<<"/"<<moving->size()<<" matches without normal constraint";
  return matches;
}

// ////////////////////////////////////////////////////////////////////////////////

// void Matcher::VisualizeMatches(
//   const std::vector<match_conf>& matches, 
//   const CloudPtr& source,
//   const CloudPtr& query,
//   CloudPtr& _colored_source,
//   CloudPtr& _colored_query)
// {

//   // Map RGB cube to bounding box
//   mesh_processing::BBox bb(source);
//   double min = 0;
//   double max = 255;

//   // Initialize inputs
//   *_colored_source = pcl::PointCloud<pcl::PointXYZRGBNormal>(*source);
//   *_colored_query = pcl::PointCloud<pcl::PointXYZRGBNormal>(*query);

//   for (size_t pt = 0; pt < source->size(); pt++)
//   {
//     double step_x = (source->points[pt].x - bb.min_x)/(bb.max_x - bb.min_x);
//     double step_y = (source->points[pt].y - bb.min_y)/(bb.max_y - bb.min_y);
//     double step_z = (source->points[pt].z - bb.min_z)/(bb.max_z - bb.min_z);

//     double r = utils::LERP(min,max,step_x);
//     double g = utils::LERP(min,max,step_y);
//     double b = utils::LERP(min,max,step_z);

//     // set colour in source cloud
//     _colored_source->points[pt].r = r;
//     _colored_source->points[pt].g = g;
//     _colored_source->points[pt].b = b;

//     // set colour according to found match
//     _colored_query->points[matches[pt].first].r = r;
//     _colored_query->points[matches[pt].first].g = g;
//     _colored_query->points[matches[pt].first].b = b;
//   }
// }

////////////////////////////////////////////////////////////////////////////////

void Matcher::VisualizeMatches(
  const std::vector<match_conf>& matches, 
  const pcl::PolygonMesh& source,
  const pcl::PolygonMesh& query,
  pcl::PolygonMesh& _colored_source,
  pcl::PolygonMesh& _colored_query)
{

  // Initialize output if not already
  _colored_source = source;
  _colored_query = query;

  // pull cloud data
  CloudPtr c_source(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  CloudPtr c_query(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  CloudPtr c_col_source(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  CloudPtr c_col_query(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  pcl::fromPCLPointCloud2(source.cloud,*c_source);
  pcl::fromPCLPointCloud2(query.cloud, *c_query);
  pcl::fromPCLPointCloud2(_colored_source.cloud,*c_col_source);
  pcl::fromPCLPointCloud2(_colored_query.cloud,*c_col_query);
  
  // colorize
  VisualizeMatches(matches,
    c_source,
    c_query,
    c_col_source,
    c_col_query);

  // put back onto mesh
  pcl::toPCLPointCloud2(*c_col_source,_colored_source.cloud);
  pcl::toPCLPointCloud2(*c_col_query,_colored_query.cloud);
}

////////////////////////////////////////////////////////////////////////////////

void Matcher::VisualizeMatches(
  const std::vector<match_conf>& matches, 
  const CloudPtr& source,
  const CloudPtr& query,
  CloudPtr& _colored_source,
  CloudPtr& _colored_query)
{

  // Map RGB cube to bounding box
  mesh_processing::BBox bb(query);
  double min = 0;
  double max = 255;

  // Initialize outputs
  *_colored_source = pcl::PointCloud<pcl::PointXYZRGBNormal>(*source);
  *_colored_query = pcl::PointCloud<pcl::PointXYZRGBNormal>(*query);
  for (size_t pt = 0; pt < query->size(); pt++)
  {
      double step_x = (query->points[pt].x - bb.min_x)/(bb.max_x - bb.min_x);
      double step_y = (query->points[pt].y - bb.min_y)/(bb.max_y - bb.min_y);
      double step_z = (query->points[pt].z - bb.min_z)/(bb.max_z - bb.min_z);

      double r = utils::LERP(min,max,step_x);
      double g = utils::LERP(min,max,step_y);
      double b = utils::LERP(min,max,step_z);

      // set colour in query cloud
      _colored_query->points[pt].r = r;
      _colored_query->points[pt].g = g;
      _colored_query->points[pt].b = b;
  }
  for (size_t pt = 0; pt < source->size(); pt++)
  {
    _colored_source->points[pt].r = 255.0;
    _colored_source->points[pt].g = 255.0;
    _colored_source->points[pt].b = 255.0;
  }

  for (size_t pt = 0; pt < source->size(); pt++)
  {
    if (matches[pt].first == -1)
    {
      // no match, set to white
    _colored_source->points[pt].r = 255.0;
    _colored_source->points[pt].g = 255.0;
    _colored_source->points[pt].b = 255.0;
    std::cerr<<"no matches"<<"\n";
      continue;
    } 
    else
    {
      _colored_source->points[pt].r = _colored_query->points[matches[pt].first].r;
      _colored_source->points[pt].g = _colored_query->points[matches[pt].first].g;
      _colored_source->points[pt].b = _colored_query->points[matches[pt].first].b;

    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Matcher::VisualizeMatches(
  const std::vector<int>& matches, 
  const CloudPtr& source,
  const CloudPtr& query,
  CloudPtr& _colored_source,
  CloudPtr& _colored_query)
{

  // Map RGB cube to bounding box
  mesh_processing::BBox bb(query);
  double min = 0;
  double max = 255;

  // Initialize outputs
  *_colored_source = pcl::PointCloud<pcl::PointXYZRGBNormal>(*source);
  *_colored_query = pcl::PointCloud<pcl::PointXYZRGBNormal>(*query);
  for (size_t pt = 0; pt < query->size(); pt++)
  {
      double step_x = (query->points[pt].x - bb.min_x)/(bb.max_x - bb.min_x);
      double step_y = (query->points[pt].y - bb.min_y)/(bb.max_y - bb.min_y);
      double step_z = (query->points[pt].z - bb.min_z)/(bb.max_z - bb.min_z);

      double r = utils::LERP(min,max,step_x);
      double g = utils::LERP(min,max,step_y);
      double b = utils::LERP(min,max,step_z);

      // set colour in query cloud
      _colored_query->points[pt].r = r;
      _colored_query->points[pt].g = g;
      _colored_query->points[pt].b = b;
  }
  for (size_t pt = 0; pt < source->size(); pt++)
  {
    _colored_source->points[pt].r = 255.0;
    _colored_source->points[pt].g = 255.0;
    _colored_source->points[pt].b = 255.0;
  }

  for (size_t pt = 0; pt < source->size(); pt++)
  {
    if (matches[pt] == -1)
    {
      // no match, set to white
    _colored_source->points[pt].r = 255.0;
    _colored_source->points[pt].g = 255.0;
    _colored_source->points[pt].b = 255.0;
      continue;
    } 
    else
    {
      _colored_source->points[pt].r = _colored_query->points[matches[pt]].r;
      _colored_source->points[pt].g = _colored_query->points[matches[pt]].g;
      _colored_source->points[pt].b = _colored_query->points[matches[pt]].b;

    }
  }
}

////////////////////////////////////////////////////////////////////////////////