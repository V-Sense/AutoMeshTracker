
#include "keyframer.h"
#include "log.h"

////////////////////////////////////////////////////////////////////////////////
// Determine keyframe locations.

std::vector<int> KeyFramer::GetKeyframeIndices(const std::string& _filename){

  // Read mesh list
  std::vector<std::string> mesh_names = ReadMeshList(_filename);

  std::vector<int> kf_indices;

  return kf_indices;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::string> KeyFramer::ReadMeshList(const std::string& _filename){

  std::vector<std::string> mesh_names_;

  std::ifstream inputfile(_filename.c_str());

  if (inputfile.is_open()){
    std::string line;
    do {
      std::getline(inputfile,line);
      if (line.empty()) continue;
      // std::cerr << line << std::endl;
      mesh_names_.push_back(line);

    } while (!inputfile.eof());

  } else {
    std::string e = "\t In KeyFramer::ReadMeshList: couldn't open input file";
    throw std::runtime_error(e);
  }

  BOOST_LOG_TRIVIAL(info) << "Number of meshes in list: " << mesh_names_.size();
  return mesh_names_;
}

////////////////////////////////////////////////////////////////////////////////

double KeyFramer::ComputeMeshArea(const pcl::PolygonMesh& mesh) {

  double totalarea = 0.0;

  pcl::PointCloud<pcl::PointXYZRGBNormal> mesh_cloud;
  pcl::fromPCLPointCloud2 (mesh.cloud, mesh_cloud);

  // #pragma omp parallel for
  for (int i = 0; i < mesh.polygons.size(); i++){
    const pcl::PointXYZRGBNormal pclv0 = mesh_cloud.points[mesh.polygons[i].vertices[0]];
    const pcl::PointXYZRGBNormal pclv1 = mesh_cloud.points[mesh.polygons[i].vertices[1]];
    const pcl::PointXYZRGBNormal pclv2 = mesh_cloud.points[mesh.polygons[i].vertices[2]];
    const Eigen::Vector3f v0 = pclv0.getVector3fMap();
    const Eigen::Vector3f v1 = pclv1.getVector3fMap();
    const Eigen::Vector3f v2 = pclv2.getVector3fMap();
    Eigen::Vector3f A, B;
    A = v1 - v0;
    B = v2 - v0;
    Eigen::Vector3f eignormal = A.cross(B);
    totalarea += 0.5 * eignormal.norm();
  }

  return totalarea;
}

////////////////////////////////////////////////////////////////////////////////

unsigned int KeyFramer::ComputeMeshGenus(const Mesh& mesh) {

    const int v = mesh.sizeVertices();
    const int f = mesh.sizeFaces();
    const int e = mesh.sizeEdges();

    // return 1 - 0.5*(v-e+f);
    return 0.5*(e-v-f+2);
    // return v + f - e + 2;
}

////////////////////////////////////////////////////////////////////////////////

unsigned int KeyFramer::ComputeMeshGenus(const pcl::PolygonMesh& mesh) {

    Mesh he_mesh;
    pcl::geometry::toHalfEdgeMesh(mesh, he_mesh);
    return ComputeMeshGenus(he_mesh);

}

////////////////////////////////////////////////////////////////////////////////

  void KeyFramer::GenerateSPHDescriptors
  (
    const std::string & meshlist,
    const std::string & exe
  )
  {

    BOOST_LOG_TRIVIAL(info) << "Generating SPH descriptors";
    std::vector<std::string> meshnames = ReadMeshList(meshlist);
    if (meshnames.empty())
    {
      std::string e 
        = "In KeyFramer::GenerateSPHDescriptors: couldn't load meshlist";
      BOOST_LOG_TRIVIAL(error) << e;
      throw std::runtime_error(e);
    }

    // get save location
    size_t path_end = meshnames[0].find_last_of("/\\");
    std::string outpath = meshnames[0].substr(0, path_end);

    // remove any existing file
    std::string cmd = "rm "+outpath+"descriptors.csv";
    system(cmd.c_str());
  
    for (size_t mesh_idx = 0; mesh_idx < meshnames.size(); mesh_idx++)
    {
      // generate descriptor
      cmd = exe+" "+meshnames[mesh_idx]+" tmp_out.txt";
      system(cmd.c_str());

      // cleanup and append to csv file
      cmd = "echo $( sed -n '2p' tmp_out.txt) >> "+ outpath + "descriptors.csv";
      system(cmd.c_str());
    } 

    // remove temp file
    cmd = "rm tmp_out.txt";
    system(cmd.c_str());
  }

////////////////////////////////////////////////////////////////////////////////

void KeyFramer::GenerateFeasScore(
  std::vector<std::string> meshlist,
  std::string outpath)
{

  std::vector<int> genii;
  std::vector<double> areas;
  int max_genus = 0;
  double max_area = 0.0;

  std::ofstream ofs_debug;
  std::string debug_path = outpath + "/result/abs/debug_feas.txt";
  ofs_debug.open(debug_path);

  MeshProcessing mp;

  // get params
  for (size_t mesh_id = 0; mesh_id < meshlist.size(); mesh_id++)
  {
    pcl::PolygonMesh mesh;
    pcl::io::loadPolygonFilePLY(meshlist[mesh_id],mesh);

    // pre-clean
    mesh = mp.RemoveSmallComponents(mesh);
    mesh = mp.FixManifoldness(mesh);

    int genus = ComputeMeshGenus(mesh);
    double area = ComputeMeshArea(mesh);
    genii.push_back(genus);
    areas.push_back(area);
    ofs_debug<<"i: "<<mesh_id<<"\n";
    ofs_debug<<"genus: "<<genus<<"\n";
    ofs_debug<<"area: "<<area<<"\n";
    
    max_genus = genus > max_genus ? genus : max_genus;
    max_area = area > max_area ? area : max_area;

  }
  ofs_debug.close();

  // get score
  std::vector<double> scores;
  for (size_t i = 0; i < meshlist.size(); i++)
  {
    double score = 1 + (max_genus - genii[i]) + (areas[i]/(max_area + 1) );
    
    scores.push_back(score);

  }

  // normalize score and save
  std::string feas_path = outpath + "/result/abs/feas.txt";
  std::ofstream ofs;
  ofs.open(feas_path);
  double max_feas = *std::max_element(scores.begin(), scores.end());
  double min_feas = *std::min_element(scores.begin(), scores.end());
  for (size_t i = 0; i < meshlist.size(); i++)
  {
    double score = (scores[i] - min_feas)/(max_feas - min_feas);
    ofs<<score<<"\n";
  }
  ofs.close();
}

////////////////////////////////////////////////////////////////////////////////

void KeyFramer::GenKeyframesAndRegions
(
  const std::vector<std::string> & meshlist
)
{
  size_t num_meshes = meshlist.size();

  size_t path_end = meshlist[0].find_last_of("/\\");
  const std::string meshpath = meshlist[0].substr(0,path_end);

  // for each frame in meshlist
  for (size_t mesh_idx = 0; mesh_idx < num_meshes; mesh_idx++)
  {
    // calculate shd


  }

  // execute matlab script

}