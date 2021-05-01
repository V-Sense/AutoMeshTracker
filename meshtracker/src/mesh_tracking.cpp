#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <typeinfo>

#include "mesh_tracking.h"

void _save_pcd2ply(
    std::vector<Eigen::Vector3d> node_pos,
    std::vector<Eigen::Vector3d> node_norm, 
    const char* filename) 
{
    assert (node_pos.size() == node_norm.size());
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "warning : unable to open %s when saving pcd to ply file.\n", filename);
        return;
    }
    fprintf(fp, "ply\n");
    fprintf(fp, "format ascii 1.0\n");
    fprintf(fp, "element vertex %d\n", node_pos.size());
    fprintf(fp, "property float x\n");
    fprintf(fp, "property float y\n");
    fprintf(fp, "property float z\n");
    fprintf(fp, "property float nx\n");
    fprintf(fp, "property float ny\n");
    fprintf(fp, "property float nz\n");
    fprintf(fp, "property uchar red\n");
    fprintf(fp, "property uchar green\n");
    fprintf(fp, "property uchar blue\n");
    fprintf(fp, "end_header\n");
    
    for (size_t i=0; i<node_pos.size(); ++i) {
        fprintf(fp, "%f %f %f %f %f %f %d %d %d\n", 
            node_pos[i][0], node_pos[i][1], node_pos[i][2],
            node_norm[i][0], node_norm[i][1], node_norm[i][2],
            0, 255, 255);
    }
    fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////

int MeshTracking::ReadMeshList(const std::string &_fileName)
{
  std::ifstream inputfile(_fileName.c_str());
  if (inputfile.is_open())
  {
    std::string line;
    do
    {
      std::getline(inputfile, line);
      if (line.empty())
        continue;
      meshNames_.push_back(line);

    } while (!inputfile.eof());
  }
  else
  {
    // BOOST_LOG_TRIVIAL(error) << "Mesh sequence reader: couldn't open file " << _fileName;
    return -1;
  }

  BOOST_LOG_TRIVIAL(info) << "Number of meshes in list: " << meshNames_.size();
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int MeshTracking::ReadRegions(const std::string &file_name)
{
  if (meshNames_.empty())
  {
    BOOST_LOG_TRIVIAL(error) << "Load mesh sequence before finding keyframes!";
    return -1;
  }

  // regions.txt file format: r_start,r_end,kf (per line)
  std::ifstream inputfile(file_name.c_str());

  if (inputfile.is_open())
  {

    std::string line;

    while (std::getline(inputfile, line))
    {
      if (line.empty())
        continue;

      unsigned int r_start, r_end;
      std::stringstream linestream(line);
      std::string value;

      //region start
      std::getline(linestream,value,',');
      r_start = std::stoi(value);

      //region end
      std::getline(linestream,value,',');
      r_end = std::stoi(value);

      regions_.push_back(std::make_pair(r_start,r_end));

      //keyframe
      std::getline(linestream,value,',');
      keyframeIndices_.push_back(std::stoi(value));
    }

    // Print regions
    BOOST_LOG_TRIVIAL(info) << "Regions: ";
    for (size_t region = 0; region < regions_.size(); region++)
    {
      BOOST_LOG_TRIVIAL(info) 
        <<regions_[region].first<< " "<<regions_[region].second<<", kf:"
        <<keyframeIndices_[region];
    }
    

    // Safety checks
    if(regions_.begin()->first != 0)
    {
      BOOST_LOG_TRIVIAL(warning) << "Warning: Check regions.txt.\n"
        << "\tFirst index must always be 0";
      regions_.begin()->first = 0;
    }

    if(regions_.back().second + 1 < meshNames_.size())
    {
      BOOST_LOG_TRIVIAL(warning) 
        << "Warning: End of specified regions is less than input meshes!"
        <<"["<<regions_.back().second + 1<<","<<meshNames_.size()<<"]\n"
        << "\tEnd of last region will be extended to end of sequence.";
        regions_.end()->second = meshNames_.size() - 1;
    }
  }
  else
  {

    BOOST_LOG_TRIVIAL(error) <<  "Keyframe reader: couldn't open file " << file_name;
    return -1;
  }


  return 1;
}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshTracking::LoadMeshFromFile(const std::string &filename)
{

  pcl::PolygonMesh mesh;
  std::string suf = filename.substr(filename.length()-3,filename.length());

  if (suf == "ply")
  {
    try
    {
      pcl::io::loadPolygonFilePLY(filename, mesh);
    }
    catch(const std::exception& e)
    {
      std::string ms = "Mesh loader: couldn't open file " + filename + "\n";
      BOOST_LOG_TRIVIAL(error) << ms;
      BOOST_LOG_TRIVIAL(error) << e.what();
      throw std::runtime_error(ms);
    }
  } 
  else if (suf == "obj")
  {
    try
    {
      pcl::io::loadPolygonFileOBJ(filename, mesh);
    }
    catch(const std::exception& e)
    {
      std::string ms = "Mesh loader: couldn't open file " + filename + "\n";
      BOOST_LOG_TRIVIAL(error) << ms;
      BOOST_LOG_TRIVIAL(error) << e.what();
      throw std::runtime_error(ms);
    }
  } 
  else
  {
    std::string ms = "Mesh loader: Unrecognised filetype, \'" + suf + "\'\n";
    BOOST_LOG_TRIVIAL(error) << ms;
    throw std::runtime_error(ms);
  }
  

  return mesh;
}

////////////////////////////////////////////////////////////////////////////////

void MeshTracking::GenAbslayersFromList
(
  const std::string & out_path
)
{

  // Check meshlist is populated
  if (meshNames_.empty())
  {
    std::string e = "In MeshTracking::GenAbslayersFromList: meshlist is empty";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }
  
  // Create a directory to store abs layers
  std::string abs_path = out_path + "abs/";
  if (mkdir(abs_path.c_str(), 0777) == -1)
  {
    BOOST_LOG_TRIVIAL(warning) 
      << "Couldn't create abs_layer storage directory at:\n" 
      << abs_path
      << "\n path may already exist...";
  }

  BOOST_LOG_TRIVIAL(info) 
    << "Generating abs layers";

  // io lock for atomic save
  omp_lock_t iolock;
  omp_init_lock(&iolock);

  // Pre-define out save paths
  std::vector<std::string> out_names;
  for (size_t i = 0; i < meshNames_.size(); i++)
  {
    size_t path_end = meshNames_[i].find_last_of("/\\");
    std::string f_name = 
      meshNames_[i].substr(path_end+1, meshNames_[i].size());
    f_name = f_name.substr(0, f_name.size()-4);
    std::string outname = abs_path + f_name + "_a.ply";
    out_names.push_back(outname);
  }

  MeshProcessing mp;

  #pragma omp parallel for
  for (size_t mesh_idx = 0; mesh_idx < meshNames_.size(); mesh_idx++)
  {  
    // omp_set_lock(&iolock);
    // pcl::PolygonMesh in_mesh = LoadMeshFromFile(meshNames_[mesh_idx]);
    // omp_unset_lock(&iolock);

    // pcl::PolygonMesh abs_mesh = mp.CreateSmoothedAbstractionLayer(in_mesh,5000,0);
    std::string remesh_cmd 
      = "/MeshTracker/resource/scripts/poisson_6.sh " + meshNames_[mesh_idx] + " " +
        out_names[mesh_idx];
    system(remesh_cmd.c_str());

    // // save
    // omp_set_lock(&iolock);
    // size_t ret = pcl::io::savePLYFile(out_names[mesh_idx], abs_mesh);
    // omp_unset_lock(&iolock);
    
    // if ( ret < 0)
    // {
    //   std::string e = "MeshTracking::GenAbslayersFromList: failed to save mesh";
    //   BOOST_LOG_TRIVIAL(error) << e;
    //   throw std::runtime_error(e);
    // }
  }

  // some redundancy but it'll still be an absolute path
  std::string cmd = "ls -d $PWD/"+abs_path+"/*_a.ply > "+abs_path+"meshlist.txt"; 
  system(cmd.c_str());

}

////////////////////////////////////////////////////////////////////////////////

void MeshTracking::TrackMeshSequence
(
  icp::IcpOptions icp_ops, 
  std::string    out_path
)
{

  // Create a directory to store transformations
  std::string t_form_path = out_path + "tforms";
  if (mkdir(t_form_path.c_str(), 0777) == -1)
  {
    BOOST_LOG_TRIVIAL(warning) 
      << "Couldn't create tform storage directory at:\n" 
      << t_form_path
      << "\n path may already exist...";
  }

  for (unsigned int r = 0; r < regions_.size(); r++)
  {

    BOOST_LOG_TRIVIAL(info) << "begin region: " << r;

    // store region tforms
    std::string r_t_form_path = t_form_path + "/r_"+std::to_string(r)+"/";
    if (mkdir(r_t_form_path.c_str(), 0777) == -1)
    {
      BOOST_LOG_TRIVIAL(warning) 
        << "Couldn't create tform storage directory at:\n" 
        << r_t_form_path
        << "\n path may already exist...";
    }

    // We do two passes:
    // KF -> Region start (bwd)
    // KF -> Region end (fwd)

    // sequence indices
    unsigned int init = keyframeIndices_[r]; //first kf_indx
    unsigned int end1 = regions_[r].first;   // start of region
    unsigned int end2 = regions_[r].second;  // end of region

    // region indices
    unsigned int end1_i = 0;
    unsigned int end2_i = end2 - end1; // length of region
    unsigned int init_i = init - end1; // kf index within region

    // copy mesh names for thread instance
    std::vector<std::string> meshNames = meshNames_;
    
    // Save keyframe with suffix
    pcl::PolygonMesh kf_mesh;

    std::size_t path_end = meshNames[init].find_last_of("/\\");
    std::string f_name = meshNames[init].substr(path_end+1, meshNames[init].size());
    f_name = f_name.substr(0, f_name.size()-4);
    std::string outnamei = out_path + f_name;
    std::string filename = meshNames[init];

    std::cerr<<"load: "<<filename<<"\n";
    kf_mesh = LoadMeshFromFile(filename);

    // Clean up the keyframe
    MeshProcessing mp;

    // Only use if input is absolute garbage, this will severely reduce detail
    #if false
    size_t target_vertices = 25000;
    // kf_mesh = mp.RemeshInstantMeshes(kf_mesh,25000);
    kf_mesh = mp.SmoothMesh(kf_mesh,10,true);
    #endif
    
    kf_mesh = mp.FixManifoldness(kf_mesh);
    kf_mesh = mp.CalculateVertexNormals(kf_mesh);

    if ( pcl::io::saveOBJFile(outnamei + "_t.obj", kf_mesh) < 0)
    {
      std::cerr<<"size: "<<kf_mesh.cloud.width<<"\n";
      std::cerr<<"file: "<<outnamei + "_t.obj"<<"\n";
      throw std::runtime_error("failed to save keyframe");
    }

    pcl::PolygonMesh fixed_mesh, moving_mesh, mesh_aligned;
    moving_mesh = kf_mesh;

    // First pass runs backwards from keyframe to start of region
    for (int i = init_i - 1; i >= (int)end1_i; i--)
    {
      // target is always initialised
      fixed_mesh = LoadMeshFromFile(meshNames_[end1 + i]);

      // Apply tracking to target
      if (icp_ops.use_adaptive)
      {
        try
        {
          mesh_aligned = DoAdaptiveTracking(moving_mesh, fixed_mesh);
        }
        catch(const std::runtime_error& e)
        {
          // Just do global nonrigid icp if something went wrong
          BOOST_LOG_TRIVIAL(warning) << "Default to global registration: "<<e.what();
          mesh_aligned = TrackMesh(moving_mesh, fixed_mesh);
        }
      }
      else
      {
        mesh_aligned = TrackMesh(moving_mesh, fixed_mesh);
      }

      // Export the mesh
      std::string f_name = meshNames[i + end1].substr(path_end+1, meshNames[init].size());
      f_name = f_name.substr(0, f_name.size()-4);
      std::string outname = out_path + f_name  + "_t.obj";
      std::cerr << outname << std::endl;
      BOOST_LOG_TRIVIAL(info) << outname;
      
      if (pcl::io::saveOBJFile(outname, mesh_aligned) < 0)
      {
        std::string e = "In MeshTracking::TrackMeshSequence: Couldn't save output";
        BOOST_LOG_TRIVIAL(error) << e;
        throw std::runtime_error(e);
      }

      // store the tform
      std::string tform_fname = 
        r_t_form_path + std::to_string(end1 + i + 1)+"_"+std::to_string(end1 + i)
        + ".txt";
      StoreTform(moving_mesh,mesh_aligned,tform_fname);

      // set moving to aligned for next iteration
      moving_mesh = mesh_aligned;

    } //end bwd tracking
    
    moving_mesh = kf_mesh;

    // Second pass forward from keyframe to end of region
    for (int i = init_i + 1; i <= end2_i; i++)
    {
      fixed_mesh = LoadMeshFromFile(meshNames_[end1 + i]);
      if (icp_ops.use_adaptive)
      {
        try
        { 
          mesh_aligned = DoAdaptiveTracking(moving_mesh, fixed_mesh);
        }
        catch(const std::runtime_error& e)
        {
          BOOST_LOG_TRIVIAL(warning) << "Default to global registration: "<<e.what();
          mesh_aligned = TrackMesh(moving_mesh, fixed_mesh);
        }
      }
      else
      {
        mesh_aligned = TrackMesh(moving_mesh, fixed_mesh);
      }

      // Export Mesh
      std::string f_name = meshNames[i + end1].substr(path_end+1, meshNames[init].size());
      f_name = f_name.substr(0, f_name.size()-4);
      std::string outname = out_path + f_name  + "_t.obj";
      std::cerr << outname << std::endl;
      BOOST_LOG_TRIVIAL(info) << outname;
      
      if (pcl::io::saveOBJFile(outname, mesh_aligned) < 0)
      {
        std::string e = "In MeshTracking::TrackMeshSequence: Couldn't save output";
        BOOST_LOG_TRIVIAL(error) << e;
        throw std::runtime_error(e);
      }

      // store the tform
      std::string tform_fname = 
        r_t_form_path + std::to_string(end1 + i - 1)+"_"+std::to_string(end1 + i)
        + ".txt";
      StoreTform(moving_mesh,mesh_aligned,tform_fname);

      // set moving to aligned for next iteration
      moving_mesh = mesh_aligned;

    } // end fwd tracking

    #if false
    // LERP and Detail Synthesis
    //--------------------------------------------------------------------------
    if (r > 0)
    {
      BOOST_LOG_TRIVIAL(info) << "applying boundary LERP...";

      unsigned int src_indx, tgt_indx;

      // first region
      //---------------
      src_indx = regions_[r-1].second;
      tgt_indx = regions_[r].first;

      // get the output files
      std::string s_name = meshNames[src_indx].substr(path_end+1, meshNames[init].size());
      s_name = s_name.substr(0, s_name.size()-4);
      const std::string s1_filename = out_path + s_name  + "_t.obj";

      std::string t_name = meshNames[tgt_indx].substr(path_end+1, meshNames[init].size());
      t_name = t_name.substr(0, t_name.size()-4);
      const std::string t1_filename = out_path + t_name  + "_t.obj";

      pcl::PolygonMesh source_mesh = LoadMeshFromFile(s1_filename);
      pcl::PolygonMesh target_mesh = LoadMeshFromFile(t1_filename);

      BOOST_LOG_TRIVIAL(debug) << "r0_s: "<<s1_filename;
      BOOST_LOG_TRIVIAL(debug) << "r0_t: "<<t1_filename;

      // Get LERP target
      pcl::PolygonMesh lerp_target(source_mesh);
      // this->max_iters_=3;
      lerp_target = DoAdaptiveTracking(lerp_target,target_mesh);
      
      // store tform and apply inverse to return original pose
      std::vector<std::array<double,6>> s2d_tform;
      s2d_tform = GetTform(source_mesh,lerp_target);

      pcl::PolygonMesh detail_syn_in = lerp_target;

      // Get detail layer from lerp_target (coarse first, then balls out)
      pcl::PolygonMesh r0_back_prop_detail = GetDetailLayer(detail_syn_in,target_mesh);
      r0_back_prop_detail = GetDetailLayer(r0_back_prop_detail,target_mesh,true);

      // initialise storage for LERP states 
      // (Keep distance to KF relatively short to avoid bad memory scaling!)
      std::vector<pcl::PolygonMesh> lerp_states;
      int num_states = regions_[r-1].second - keyframeIndices_[r-1];
      BOOST_LOG_TRIVIAL(debug) << "num_states: "<< num_states;

      // Generate LERP steps
      for (size_t step = num_states; step > 0; step--)
      {
        float l_step = float(step/float(num_states)); // floaty floaty float float
        BOOST_LOG_TRIVIAL(debug) << "l% = "<<l_step;
        pcl::PolygonMesh detail_lerp = 
          mp.LERPMesh(detail_syn_in,r0_back_prop_detail,l_step);
          
        lerp_states.push_back(detail_lerp);
      }

      // Propagate detail synth ( decreasing influence )
      for (size_t step = 0; step < lerp_states.size(); step++)
      {
        // load current mesh
        pcl::PolygonMesh curr_state_mesh = lerp_states[step];

        // Remove the synth state temporal offset
        bool invert_tform = true;
        curr_state_mesh = ApplyTform(curr_state_mesh,s2d_tform,invert_tform);
        
        // makea temp cloud to apply the tforms 
        CloudPtr inter_cloud = mp.GetCloudFromPolygonMesh(curr_state_mesh);

        // iteratively tform it until we hit the correct timestep (backwards from kf boundary)
        int kf_index = keyframeIndices_[r-1];
        int r_prop_init = regions_[r-1].second;
        int prop_limit = r_prop_init - step;
        BOOST_LOG_TRIVIAL(debug) << "step: "<<step;
        for (size_t tform_state = r_prop_init; tform_state > prop_limit; tform_state--)
        {

          BOOST_LOG_TRIVIAL(debug) << "state: "<<tform_state;

          std::string r0_tform_s_path = meshNames[tform_state-1].substr(path_end+1, meshNames[init].size());
          r0_tform_s_path = out_path + r0_tform_s_path.substr(0, r0_tform_s_path.size()-4) + "_t.obj";
          std::string r0_tform_t_path = meshNames[tform_state].substr(path_end+1, meshNames[init].size());
          r0_tform_t_path = out_path + r0_tform_t_path.substr(0, r0_tform_t_path.size()-4) + "_t.obj";

          pcl::PolygonMesh r0_tform_s = LoadMeshFromFile(r0_tform_s_path);
          pcl::PolygonMesh r0_tform_t = LoadMeshFromFile(r0_tform_t_path);
          std::vector<std::array<double,6>> r0_tforms = GetTform(r0_tform_s,r0_tform_t);

          // apply tform
          for (size_t pt = 0; pt < inter_cloud->size(); pt++)
          {

            inter_cloud->points[pt].x -= r0_tforms[pt][0];
            inter_cloud->points[pt].y -= r0_tforms[pt][1];
            inter_cloud->points[pt].z -= r0_tforms[pt][2];
            inter_cloud->points[pt].normal_x -= r0_tforms[pt][3];
            inter_cloud->points[pt].normal_y -= r0_tforms[pt][4];
            inter_cloud->points[pt].normal_z -= r0_tforms[pt][5];
          }
                  
          // tform_ifs.close();
        }

        // attach to output and return
        pcl::toPCLPointCloud2(*inter_cloud, curr_state_mesh.cloud);

        // save at correct timestep
        std::string f_name = meshNames[regions_[r-1].second-step].substr(path_end+1, meshNames[init].size());
        f_name = f_name.substr(0, f_name.size()-4);
        std::string outname = out_path + f_name  + "_t.obj";
        std::cerr << outname << std::endl;
        BOOST_LOG_TRIVIAL(info) << outname;
        
        if (pcl::io::saveOBJFile(outname, curr_state_mesh) < 0)
        {
          std::string e = "In MeshTracking::TrackMeshSequence: Couldn't save output";
          BOOST_LOG_TRIVIAL(error) << e;
          throw std::runtime_error(e);
        }
        
      }

      // second region
      //---------------
      src_indx = regions_[r].first;
      tgt_indx = regions_[r-1].second;

      // get the output files
      s_name = meshNames[src_indx].substr(path_end+1, meshNames[init].size());
      s_name = s_name.substr(0, s_name.size()-4);
      const std::string s2_filename = out_path + s_name  + "_t.obj";

      t_name = meshNames[tgt_indx].substr(path_end+1, meshNames[init].size());
      t_name = t_name.substr(0, t_name.size()-4);
      const std::string t2_filename = out_path + t_name  + "_t.obj";

      source_mesh = LoadMeshFromFile(s2_filename);
      target_mesh = LoadMeshFromFile(t2_filename);

      BOOST_LOG_TRIVIAL(debug) << "source: "<<s2_filename;
      BOOST_LOG_TRIVIAL(debug) << "target: "<<t2_filename;

      // Get LERP target
      lerp_target = source_mesh;
      lerp_target = DoAdaptiveTracking(lerp_target,target_mesh);
      // lerp_target = LoadMeshFromFile(debug_path_+"r1_lerp_target.ply");

      std::string asduasdu = debug_path_ + "r1_lerp_target.ply";
      pcl::io::savePLYFile(asduasdu,lerp_target); 

      // store tform and apply inverse to return original pose
      s2d_tform = GetTform(source_mesh,lerp_target);
      
      detail_syn_in = lerp_target;

      // Get detail layer from lerp_target
      pcl::PolygonMesh r1_back_prop_detail = GetDetailLayer(detail_syn_in,target_mesh);
      r1_back_prop_detail = GetDetailLayer(r1_back_prop_detail,target_mesh,true);
      std::string dsuasdud = debug_path_ + "r1_detail_layer.ply";
      pcl::io::savePLYFile(dsuasdud,r1_back_prop_detail); 

      // initialise storage for LERP states 
      // (Keep distance to KF relatively short to avoid bad memory scaling!)
      lerp_states.clear();
      num_states = keyframeIndices_[r] - regions_[r].first;
      BOOST_LOG_TRIVIAL(debug) << "num_states: "<< num_states;

      // Generate LERP steps
      for (size_t step = num_states; step > 0; step--)
      {
        float l_step = float(step/float(num_states)); // floaty floaty float float
        BOOST_LOG_TRIVIAL(debug) << "l% = "<<l_step;
        pcl::PolygonMesh detail_lerp = 
          mp.LERPMesh(detail_syn_in,r1_back_prop_detail,l_step);
          
        lerp_states.push_back(detail_lerp);
      }

      // Propagate detail synth ( decreasing influence )
      for (size_t step = 0; step < lerp_states.size(); step++)
      {
        // load current mesh
        pcl::PolygonMesh curr_state_mesh = lerp_states[step];

        // Remove the synth state temporal offset
        bool invert_tform = true;
        curr_state_mesh = ApplyTform(curr_state_mesh,s2d_tform,invert_tform);
        
        // makea temp cloud to apply the tforms 
        CloudPtr inter_cloud = mp.GetCloudFromPolygonMesh(curr_state_mesh);

        // iteratively tform it until we hit the correct timestep (backwards from kf boundary)
        int kf_index = keyframeIndices_[r];
        int f_prop_init = regions_[r].first;
        int f_prop_limit = f_prop_init + step;
        BOOST_LOG_TRIVIAL(debug) << "step: "<<step;
        for (size_t tform_state = f_prop_init; tform_state < f_prop_limit; tform_state++)
        {

          BOOST_LOG_TRIVIAL(debug) << "state: "<<tform_state;

          std::string r1_tform_s_path = meshNames[tform_state+1].substr(path_end+1, meshNames[init].size());
          r1_tform_s_path = out_path + r1_tform_s_path.substr(0, r1_tform_s_path.size()-4) + "_t.obj";
          std::string r1_tform_t_path = meshNames[tform_state].substr(path_end+1, meshNames[init].size());
          r1_tform_t_path = out_path + r1_tform_t_path.substr(0, r1_tform_t_path.size()-4) + "_t.obj";

          pcl::PolygonMesh r1_tform_s = LoadMeshFromFile(r1_tform_s_path);
          pcl::PolygonMesh r1_tform_t = LoadMeshFromFile(r1_tform_t_path);
          std::vector<std::array<double,6>> r1_tforms = GetTform(r1_tform_s,r1_tform_t);

           // apply tform
          for (size_t pt = 0; pt < inter_cloud->size(); pt++)
          {

            inter_cloud->points[pt].x -= r1_tforms[pt][0];
            inter_cloud->points[pt].y -= r1_tforms[pt][1];
            inter_cloud->points[pt].z -= r1_tforms[pt][2];
            inter_cloud->points[pt].normal_x -= r1_tforms[pt][3];
            inter_cloud->points[pt].normal_y -= r1_tforms[pt][4];
            inter_cloud->points[pt].normal_z -= r1_tforms[pt][5];
          }
        }

        // attach to output and return
        pcl::toPCLPointCloud2(*inter_cloud, curr_state_mesh.cloud);

        // save at correct timestep
        std::string f_name = meshNames[regions_[r].first+step].substr(path_end+1, meshNames[init].size());
        f_name = f_name.substr(0, f_name.size()-4);
        std::string outname = out_path + f_name  + "_t.obj";
        std::cerr << outname << std::endl;
        BOOST_LOG_TRIVIAL(info) << outname;
        
        if (pcl::io::saveOBJFile(outname, curr_state_mesh) < 0)
        {
          std::string e = "In MeshTracking::TrackMeshSequence: Couldn't save output";
          BOOST_LOG_TRIVIAL(error) << e;
          throw std::runtime_error(e);
        }
        
      }

      BOOST_LOG_TRIVIAL(info) << "end LERP for for boundary region " << r-1 <<"/"<< r;

    } // end reverse tracking 
    #endif
    //--------------------------------------------------------------------------
  } // end region 
}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshTracking::DoAdaptiveTracking
(
  const pcl::PolygonMesh &  _moving_mesh,
  const pcl::PolygonMesh &  _fixed_mesh
)
{

  time_t toc, tic; 

	MeshProcessing mp;
	CloudProcessing cp;
	MeshSegmentation ms;
	Matcher mt;

	// Copy input and clean up the target mesh
  pcl::PolygonMesh moving_mesh(_moving_mesh);
	pcl::PolygonMesh fixed_mesh = mp.FixManifoldness(_fixed_mesh);
	fixed_mesh = mp.RemoveSmallComponents(fixed_mesh,0.01); // as percent wrt input size
	fixed_mesh = mp.CalculateVertexNormals(fixed_mesh);

	BOOST_LOG_TRIVIAL(debug) << "Started Adaptive tracking...";
	int input_size = moving_mesh.cloud.width*moving_mesh.cloud.height;
	if (input_size > 50000)
	{
    BOOST_LOG_TRIVIAL(warning) << "Input mesh is too large ("
      <<moving_mesh.cloud.width<<"), decimating to 50k";
    moving_mesh = mp.DecimateMesh(moving_mesh,50000);
	}

	if (input_size < 2000)
	{
    BOOST_LOG_TRIVIAL(warning) 
      << "Input mesh too small. Defaulting to global registration";
    return TrackMesh(moving_mesh, fixed_mesh);
	}

  GNParams solver_params;

  //----------------------------------------------------------------------------
  // Main Loop
  //----------------------------------------------------------------------------

  while(true)
  {
    tic = time(0); 
    BOOST_LOG_TRIVIAL(debug) << "\n\n";
    BOOST_LOG_TRIVIAL(debug) << "Beginning iteration: "<<curr_iter_;

    TriMesh mesh_in(moving_mesh);
    TriMesh mesh_target(fixed_mesh);

    //--------------------------------------------------------------------------
    // Get Correspondences
    //--------------------------------------------------------------------------


    // pull cloud data, using MeshProcessing method because pcl::FromPointCloud2
    // discards normals
    CloudPtr fixed_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
    CloudPtr moving_cloud (new pcl::PointCloud<pcl::PointXYZRGBNormal>());
    fixed_cloud = mp.GetCloudFromPolygonMesh(fixed_mesh);
    moving_cloud = mp.GetCloudFromPolygonMesh(moving_mesh);

    // Get dense correspondences from abstraction layer. 
    pcl::PolygonMesh corr_layer_aligned, corr_layer_input;
    CloudPtr corr_aligned_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
    CloudPtr corr_input_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());

    GetDenseCorrespondences(
      moving_mesh,
      fixed_mesh,
      corr_layer_aligned,
      corr_layer_input);
    
    corr_aligned_cloud = mp.GetCloudFromPolygonMesh(corr_layer_aligned);
    corr_input_cloud = mp.GetCloudFromPolygonMesh(corr_layer_input);
    
    toc = time(0); 
    BOOST_LOG_TRIVIAL(debug) << "Dense Correspondences found, "<< difftime(toc, tic) <<"(s)";
    tic = time(0);

    //--------------------------------------------------------------------------
    // Build Graph
    //--------------------------------------------------------------------------

    pcl::PolygonMesh graph_mesh = corr_layer_input;
    CloudPtr graph_cloud = mp.GetCloudFromPolygonMesh(graph_mesh);

    // Partial initialization of graph object
    DeformGraph graph;
    graph.SetSampleRadius(mp.CalculateAverageEdgeLength(graph_mesh));
    graph.BuildGraph(
      mesh_in,
      graph_mesh, 
      kDG_KNN);

    // Graph Densification
    //--------------------------------------------------------------------------
    BOOST_LOG_TRIVIAL(debug) << 
      "Performing graph densification near distant vertices";

    typedef std::pair<double,int> npair;
    std::vector<npair> distances;
    double radius =  2.0 * mp.CalculateAverageEdgeLength(graph_mesh); // 1.5
    double qualify = 0.2 * radius; // e.g. radius/knn //0.2

    // for each vertex which lies outside the qualifying radius, 
    // add a point to the graph cloud which lies halfway between vertex and 
    // nearest existing node. 
    // Then we'll call meshlab to perform screened poisson on the denser cloud,
    // followed by a uniform remeshing from instant meshes
    for (size_t vid = 0; vid < moving_cloud->size(); vid++)
    {
      pcl::PointXYZRGBNormal source(moving_cloud->points[vid]);
      std::vector<double> dists_vid;

      for (size_t qvid = 0; qvid < graph_cloud->size(); qvid++)
      {
        pcl::PointXYZRGBNormal query(graph_cloud->points[qvid]);
        double sq_dist = 0.0;
        sq_dist += (source.x - query.x)*(source.x - query.x);
        sq_dist += (source.y - query.y)*(source.y - query.y);
        sq_dist += (source.z - query.z)*(source.z - query.z); 
        double dist = sqrt(sq_dist);
        
        if (dist <= radius)
        {
          dists_vid.push_back(dist);
        }
        
      }

      if (dists_vid.size() == 0)
      {
        throw std::runtime_error("Graph Densification: no points within radius");
      }
      
      // sort by distance
      std::sort(dists_vid.begin(), dists_vid.end());

      npair dist_vid(dists_vid[0],vid);
      distances.push_back(dist_vid);
    }
    // sort by distance
    std::sort(distances.begin(), distances.end(),
    [](const npair& l, const npair& r) {
      if (l.first != r.first){
        return l.first < r.first;}
      return l.second < r.second;
    }); 

    int max_dist_vid = distances.back().second;      
    
    int idx = distances.size()-1;
    while(true)
    {
      if (distances[idx].first < qualify)
      {
        break;
      }
              
      pcl::PointXYZRGBNormal pt_max(moving_cloud->points[distances[idx].second]);
      pcl::PointXYZRGBNormal tmp_pt(pt_max);
      std::vector<int> max2graph = mt.GetConstrainedMatches(tmp_pt,graph_cloud,0.0);
      if (max2graph.empty()) // should never be empty
      {
        // std::cerr<<"distance idx: "<<distances[idx].second<<"\n";
        // std::cerr<<"tmp_pt: "<<pt_max<<"\n";
        std::string e = "\tGraph Densification: max2graph search returned empty.\n";
        throw std::runtime_error(e);
      }
      
      tmp_pt.x = 0.5*(graph_cloud->points[max2graph[0]].x + pt_max.x);
      tmp_pt.y = 0.5*(graph_cloud->points[max2graph[0]].y + pt_max.y);
      tmp_pt.z = 0.5*(graph_cloud->points[max2graph[0]].z + pt_max.z);
      tmp_pt.g = 255.0;
      graph_cloud->push_back(tmp_pt);
      --idx;
    }
  
    // save
    std::string f_dense_graph_out = debug_path_+"densified_graph.ply";
    pcl::io::savePLYFile(f_dense_graph_out,*graph_cloud);
  
    // meshlab poisson and remesh 
    //(we use meshlab because PCL produces some...interesting results)
    std::string f_final_graph = debug_path_ + "uniform_dense_graph.ply";
    std::string remesh_cmd 
      = "/MeshTracker/resource/scripts/remesh_graph.sh " + f_dense_graph_out + " " +
        f_final_graph + " " + 
        std::to_string(0.75*mp.CalculateAverageEdgeLength(graph_mesh)); //0.95
    
    system(remesh_cmd.c_str());
    pcl::io::loadPolygonFilePLY(f_final_graph,graph_mesh);
    graph_mesh = mp.FixManifoldness(graph_mesh);
    graph_mesh = mp.RemoveSmallComponents(graph_mesh,0.01);
    pcl::io::savePLYFile(f_final_graph,graph_mesh);
    //--------------------------------------------------------------------------
    
    // Use graph mesh to build the deform graph object
    graph.BuildGraph(
      mesh_in,
      graph_mesh, 
      kDG_KNN);

    toc = time(0); 
    BOOST_LOG_TRIVIAL(debug) << "Graph Built, "<< difftime(toc, tic) <<"(s)";
    tic = time(0);
    
    int iter = 0;
    double pre_energy = std::numeric_limits<double>::max();

    // Build constraint data
    NodeCorrs cons = GetCorresConstraints(
      graph,
      corr_aligned_cloud,
      corr_input_cloud,
      fixed_cloud);

    toc = time(0); 
    BOOST_LOG_TRIVIAL(debug) << "constraints set, "<< difftime(toc, tic) <<"(s)";
    tic = time(0);

    // Get constraint-node neighbours and weights
    std::vector<pcl::PointXYZRGBNormal> cons_src;
    std::vector<std::vector<int>> con_node_neigh;

    for (size_t con_idx = 0; con_idx < cons.size(); con_idx++)
    {
      cons_src.push_back(cons[con_idx].first);
    }

    con_node_neigh = graph.GetConNodeNeighbours(
      cons_src,num_vert_node_neighbors_);

    toc = time(0); 
    BOOST_LOG_TRIVIAL(debug) << "con neighbours set, "<< difftime(toc, tic) <<"(s)";
    tic = time(0);

    std::vector<std::vector<double>> con_node_weights = 
      graph.CalculateVertNodeWeights(cons_src, con_node_neigh);

    toc = time(0); 
    BOOST_LOG_TRIVIAL(debug) << "con weights set, "<< difftime(toc, tic) <<"(s)";
    tic = time(0);

    // Build vertex-node influence list
    std::vector<std::vector<int>> node_influence_list = 
      graph.GetVertNodeNeighbours(moving_mesh);

    toc = time(0); 
    BOOST_LOG_TRIVIAL(debug) << "vert neighbours set, "<< difftime(toc, tic) <<"(s)";
    tic = time(0);

    std::vector<std::vector<double>> vert_node_weights = 
      graph.CalculateVertNodeWeights(moving_mesh, node_influence_list);

    toc = time(0); 
    BOOST_LOG_TRIVIAL(debug) << "vert weights set, "<< difftime(toc, tic) <<"(s)";
    tic = time(0);
    
    //--------------------------------------------------------------------------
    // Perfom Solve
    //--------------------------------------------------------------------------

    BOOST_LOG_TRIVIAL(debug) << "Nodes: "<<graph.node_pos.size();
    BOOST_LOG_TRIVIAL(debug) << "Constraints: "<<cons.size();

    curr_energy_ = graph.OptimizeOnce(
      graph,
      solver_params,
      con_node_neigh,
      con_node_weights,
      cons);

    toc = time(0); 
    BOOST_LOG_TRIVIAL(debug) << "Graph Solve Finished, "<< difftime(toc, tic) <<"(s)";
    tic = time(0);

    //--------------------------------------------------------------------------
    // Perfom Deform
    //--------------------------------------------------------------------------

    TriMesh deformed_mesh(moving_mesh);

    // Actual deform operation
    graph.Deform(
      node_influence_list,
      vert_node_weights,
      graph_mesh,
      deformed_mesh);

    toc = time(0); 
    BOOST_LOG_TRIVIAL(debug) << "Deformation Applied, "<< difftime(toc, tic) <<"(s)";
    tic = time(0);

    // Put deformed data points back onto the pcl mesh
    CloudPtr tmp_cloud = deformed_mesh.ToPclPointCloud();
    pcl::toPCLPointCloud2(*tmp_cloud,moving_mesh.cloud);

    // Final check
    if (mp.HasNanData(moving_mesh))
    {
      std::string e = "In Main Loop: output contains NaN data";
      BOOST_LOG_TRIVIAL(error) << e;
      throw std::runtime_error(e);
    }

    //--------------------------------------------------------------------------
    // Parameter Optimization For Next Iter
    //--------------------------------------------------------------------------

    double delta_energy = 0.0;
    if (curr_iter_ > 0)
    {
      delta_energy = std::abs(curr_energy_ - prev_energy_) / (prev_energy_);
    }
    
    //	Relax rigid term and smooth term
    if (curr_energy_ < GN_solver_relax_thresh_)
    {
      solver_params.RelaxParams();
    }

    BOOST_LOG_TRIVIAL(info) << "Registration Energy = " << curr_energy_ 
      << ", Energy Delta = " << delta_energy;

    //	Break the loop if rigid term is less than 100 or iter > max
    if ( solver_params.m_alpha_rigid < GN_solver_rigidity_w_min || 
      curr_iter_ == max_iters_)
    {
      Reset();
      break;
    }

    prev_energy_ = curr_energy_;
    IncrementIter();

  }

	return moving_mesh;
}

////////////////////////////////////////////////////////////////////////////////

void MeshTracking::GetDenseCorrespondences(
  pcl::PolygonMesh &moving_mesh,
  pcl::PolygonMesh &fixed_mesh,
  pcl::PolygonMesh &aligned_graph_,
  pcl::PolygonMesh &in_graph_,
  const int& align_layer_size)
{

  CloudProcessing cp;
  MeshProcessing mp;
  MeshSegmentation ms;
  Matcher mt;
  time_t toc, tic; 
  tic = time(0);

  // Generate abstraction layers
  pcl::PolygonMesh moving_abs = mp.CreateSmoothedAbstractionLayer(moving_mesh,align_layer_size);
  pcl::PolygonMesh fixed_abs = mp.CreateSmoothedAbstractionLayer(fixed_mesh,align_layer_size);
 
  toc = time(0);
  BOOST_LOG_TRIVIAL(debug) 
    << "\tGenerated correspondence abstraction layers, "<< difftime(toc, tic)<<"(s)";
  tic = time(0);

  // Extract ptClouds from meshes
  CloudPtr moving_cloud = mp.GetCloudFromPolygonMesh(moving_mesh);
  CloudPtr moving_abs_cloud = mp.GetCloudFromPolygonMesh(moving_abs);
  CloudPtr fixed_cloud = mp.GetCloudFromPolygonMesh(fixed_mesh);
  CloudPtr fixed_abs_cloud = mp.GetCloudFromPolygonMesh(fixed_abs);


  // get segmentation from moving_abs mesh
  VidSidMap moving_vs_map,fixed_vs_map;
  ms.Segment(moving_abs,moving_vs_map);

  toc = time(0);
  BOOST_LOG_TRIVIAL(debug) 
    << "\tSegmented abstraction layer, "<< difftime(toc, tic)<<"(s)";
  tic = time(0);

  // initialize seg count
  SetInitNumSegs(ms.GetNumSegments());

  VidSidMap fused_moving_vs_map = moving_vs_map;

  toc = time(0);
  BOOST_LOG_TRIVIAL(debug) 
    << "\tPre-CPD, "<< difftime(toc, tic)<<"(s)";
  tic = time(0);

  // Initialize alignment with global cpd
  CPDICP tmp_cpd;
  CPD_Result tmp_res = tmp_cpd.AlignHighResClouds(
      fixed_abs_cloud,
      moving_abs_cloud,
      this->default_CPD_tolerance_);
  CloudPtr init_moving_abs_cloud = tmp_res.aligned_cloud;

  toc = time(0);
  BOOST_LOG_TRIVIAL(debug) 
    << "\tInitialised using global semi-rigid ICP, "<< difftime(toc, tic)<<"(s)";
  tic = time(0);
  
  // FIXME: this is a dumb quick fix, but the alignment kills normals
  // to fix it, reattach to mesh, then calculate normals and pull cloud out again
  pcl::PolygonMesh tmp_moving_abs(moving_abs);
  pcl::toPCLPointCloud2(*init_moving_abs_cloud,tmp_moving_abs.cloud);
  init_moving_abs_cloud = mp.GetCloudFromPolygonMesh(tmp_moving_abs);

  double dot_prod_tolerance = 0.7; //0.5
  std::vector<int> sparse_matches = 
    mt.GetConstrainedMatches(
      fixed_abs_cloud, 
      init_moving_abs_cloud,dot_prod_tolerance);

  toc = time(0);
  BOOST_LOG_TRIVIAL(debug) 
    << "\tGot normal constrained, sparse matches, "<< difftime(toc, tic)<<"(s)";
  tic = time(0);

  // get segmentation info for fixed target
  fixed_vs_map = SegmentationFromMatches(
    fixed_abs, 
    fused_moving_vs_map, 
    sparse_matches);

  toc = time(0);
  BOOST_LOG_TRIVIAL(debug) 
    << "\tInitial segment map to target, "<< difftime(toc, tic)<<"(s)";
  tic = time(0);


  int incoherence_val = -1;
  try
  {
    incoherence_val =
      SegmentProjectionCoherenceTest(fused_moving_vs_map, fixed_vs_map);
  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    ExportSegments(moving_vs_map,moving_abs,debug_path_+"segs");
    throw std::runtime_error(e.what());
  }

  // flags segment projection as incoherent for any value about -1
  // Value returned will identify the first incoherent segment
  while (incoherence_val > -1 && this->curr_num_segs_ > 2)
  {

    int seg_id = incoherence_val;

    int fused_id = 
      ms.FuseWithSmallestConnectedSegment(fused_moving_vs_map, moving_abs, seg_id);

    BOOST_LOG_TRIVIAL(debug) 
      << "\tmissing geometry detected, performing fusion for seg: " << seg_id
      <<" -> "<<fused_id;

    // saved fused id if not already stored
    std::vector<int>::iterator it = std::find( 
      this->fused_segments_.begin(),  this->fused_segments_.end(), fused_id);
    if (it == this->fused_segments_.end())
    {
      this->fused_segments_.push_back(fused_id);
    }

    // remove the initial seg_id from fused list if it exists
    it = std::find( 
      this->fused_segments_.begin(),  this->fused_segments_.end(), seg_id);
    if (it != this->fused_segments_.end())
    {
      this->fused_segments_.erase(it);
    }

    // num_segs needs to be reduced after fusion
    this->curr_num_segs_--; 

    // Recalculate the target mapping
    fixed_vs_map =
        SegmentationFromMatches(fixed_abs, fused_moving_vs_map, sparse_matches);

    incoherence_val =
        SegmentProjectionCoherenceTest(fused_moving_vs_map, fixed_vs_map);
  }

  toc = time(0);
  BOOST_LOG_TRIVIAL(debug)
    << "\tFinal Seg map transfered, "<< difftime(toc, tic)<<"(s)";
  tic = time(0);

  // initialize output
  pcl::PolygonMesh aligned_abs(moving_abs);
  CloudPtr aligned_abs_cloud = mp.GetCloudFromPolygonMesh(aligned_abs);

  // initialize most connected segment
  aligned_abs_cloud = AlignSegment(
    aligned_abs_cloud,
    fixed_abs_cloud, 
    fused_moving_vs_map, 
    fixed_vs_map, 
    0);

  // perform segment-based alignment
  for (int seg_id = this->curr_num_segs_-1; seg_id > 0; seg_id--)
  {
    aligned_abs_cloud = AlignSegment(
      aligned_abs_cloud,
      fixed_abs_cloud, 
      fused_moving_vs_map, 
      fixed_vs_map, 
      seg_id);  
  } 

  toc = time(0);
  BOOST_LOG_TRIVIAL(debug) << "\tSegments aligned, "<< difftime(toc, tic)<<"(s)";

  // alignment correction step 
  //----------------------------------------------------------------------------
  CloudPtr ag_cloud = mp.GetCloudFromPolygonMesh(aligned_graph_);
  std::vector<int> al2tgt_matches = 
    mt.GetConstrainedMatches(ag_cloud,fixed_cloud,0.85);

  double res = cp.ComputeCloudResolution(fixed_abs_cloud);
  for (size_t m = 0; m < al2tgt_matches.size(); m++)
  {
    // Ignore if part of a fused segment
    int seg_id = fixed_vs_map[al2tgt_matches[m]];
    std::vector<int>::iterator it = std::find( 
      this->fused_segments_.begin(),  this->fused_segments_.end(), seg_id);

    if (it != this->fused_segments_.end())
    {
      continue;
    }

    // ignore if nearest match is too far
    Eigen::Vector3d d; d.setZero();
    d(0) = fixed_cloud->points[al2tgt_matches[m]].x - ag_cloud->points[m].x;
    d(1) = fixed_cloud->points[al2tgt_matches[m]].y - ag_cloud->points[m].y;
    d(2) = fixed_cloud->points[al2tgt_matches[m]].z - ag_cloud->points[m].z;
    double dist = std::abs(d(0))+std::abs(d(1))+std::abs(d(2));
    if (dist < 5*res)
    {
      // Align to tangent plane
      Eigen::Vector3d delta_vec; delta_vec.setZero();
      Eigen::Vector3d norm; norm.setZero();
      norm(0) = ag_cloud->points[m].normal_x;
      norm(1) = ag_cloud->points[m].normal_y;
      norm(2) = ag_cloud->points[m].normal_z;
      double arg = d.dot(norm);
      if (arg > M_PI)
      {
        arg -= M_PI;
      }
      delta_vec = d/std::cos(arg);
      ag_cloud->points[m].x += delta_vec(0);//fixed_cloud->points[al2tgt_matches[m]];
      ag_cloud->points[m].y += delta_vec(1);
      ag_cloud->points[m].z += delta_vec(2);
    }
  }
  pcl::toPCLPointCloud2(*ag_cloud,aligned_graph_.cloud);
  //----------------------------------------------------------------------------


  // Reapply points
  pcl::toPCLPointCloud2(*aligned_abs_cloud,aligned_abs.cloud);

  // smooth out some janky seams
  // aligned_graph_ = mp.SmoothMesh(aligned_abs,25);
  aligned_graph_ = aligned_abs;
  in_graph_ = moving_abs;

  return;
}

////////////////////////////////////////////////////////////////////////////////
VidSidMap MeshTracking::SegmentationFromMatches(
    const pcl::PolygonMesh &fixed,
    const VidSidMap input_vs_map,
    const std::vector<int> &matches) const
{

  CloudPtr fixed_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  pcl::fromPCLPointCloud2(fixed.cloud, *fixed_cloud);

  VidSidMap output;
  std::vector<int> unmatched;

  bool sparse_map = false;
  for (int point = 0; point < fixed_cloud->size(); point++)
  {
    // valid matches are positive indices
    if (matches[point] >= 0)
    {
      output[point] = input_vs_map.find(matches[point])->second;
    }
    else
    {
      unmatched.push_back(point);
      output[point] = -1;      
      sparse_map = true;
    }
  }

  if (sparse_map)
  {
    BOOST_LOG_TRIVIAL(warning) << "Sparse map used in seg transfer. "
      <<"Attepting to densify before alignment!";
  }

  // Get he_mesh for target
  Mesh he_mesh;
  pcl::geometry::toHalfEdgeMesh(fixed, he_mesh);
  MeshProcessing mp;

  // For all invalid points, assign segment of nearest valid point.
  int max_iter = 15;
  for (int pt = 0; pt < unmatched.size(); pt++)
  {
    // Geodesic search for valid points
    bool found = false;
    int k = 1;
    do
    {
      std::vector<int> nn = mp.GetKRingNeighbors(fixed,unmatched[pt],k);
      for (size_t n = 0; n < nn.size(); n++)
      {
        int sid = output.find(nn[n])->second;
        if(sid>0)
        {
          output[nn[n]] = sid;
          found = true;
        }
      }
      ++k;
    }
    while(!found && k < max_iter);

    if(k==max_iter && !found)
    {
      std::vector<int> nn = mp.GetKRingNeighbors(fixed,unmatched[pt],6);
      for (size_t ch = 0; ch < nn.size(); ch++)
      {
        fixed_cloud->points[nn[ch]].r=255.0;
        fixed_cloud->points[nn[ch]].g=0.0;
        fixed_cloud->points[nn[ch]].b=0.0;
      }
      fixed_cloud->points[unmatched[pt]].r=0.0;
      fixed_cloud->points[unmatched[pt]].g=0.0;
      fixed_cloud->points[unmatched[pt]].b=255.0;

    }

    // search neighbours via kd tree and assign best seg_id
    pcl::KdTreeFLANN<pcl::PointXYZRGBNormal> kdtree;
    kdtree.setInputCloud(fixed_cloud);

    // set params
    const float search_radius = 0.05;
    std::vector<int> pointIdxRadiusSearch;
    std::vector<float> pointNKNSquaredDistance;
    pcl::PointXYZRGBNormal search_pt = fixed_cloud->points[unmatched[pt]];
    int nres = kdtree.radiusSearch(
      search_pt, 
      search_radius,
      pointIdxRadiusSearch, 
      pointNKNSquaredDistance);

    // store normal
    Eigen::Vector3d norm_s = {
      search_pt.normal_x,
      search_pt.normal_y,
      search_pt.normal_z};

    
    for (size_t res = 0; res < nres; res++)
    {
      // first check normal constraint
      Eigen::Vector3d norm_r = {
        fixed_cloud->points[pointIdxRadiusSearch[res]].normal_x,
        fixed_cloud->points[pointIdxRadiusSearch[res]].normal_y,
        fixed_cloud->points[pointIdxRadiusSearch[res]].normal_z};

      if(norm_s.dot(norm_r) > 0 )
      {
        // then check for a valid id
        int test_id = output.find(pointIdxRadiusSearch[res])->second;
        if(test_id>=0)
        {
          output[unmatched[pt]] = test_id;
          found = true;
          break;
        }
      }
    }

    if (!found)
    {
      std::string e = "In SegmentationFromMatches: no valid seg_id for point in target";
      BOOST_LOG_TRIVIAL(error) << e;
      throw std::runtime_error(e);
    }
  }

  if (output.size() <= 0)
  {
    std::string e = "In SegmentationFromMatches: couldn't populate segmap";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }
  
  return output;
}

////////////////////////////////////////////////////////////////////////////////

int MeshTracking::SegmentProjectionCoherenceTest(
  VidSidMap &moving_vs_map,
  VidSidMap &fixed_vs_map)
{

  int ret_state = -1;

  for (size_t seg_id = 0; seg_id < this->curr_num_segs_; seg_id++)
  {
    int moving_count = 0;
    int fixed_count = 0;

    // TODO: use a std count function
    for (VidSidMap::iterator it = moving_vs_map.begin(); it != moving_vs_map.end(); ++it)
    {
      if (it->second == seg_id)
        moving_count++;
    }

    for (VidSidMap::iterator it = fixed_vs_map.begin(); it != fixed_vs_map.end(); ++it)
    {
      if (it->second == seg_id)
        fixed_count++;
    }

    if (moving_count == 0)
    {
      std::string e = 
      "In SegmentProjectionCoherenceTest: failed to count seg"+std::to_string(seg_id)+" size";
      BOOST_LOG_TRIVIAL(error) << e;
      throw std::runtime_error(e);
    }

    // In the event that an entire segment is completely missing in the fixed
    if (fixed_count == 0)
    {
      BOOST_LOG_TRIVIAL(warning) 
        << "In CoherenceTest: segment entirely missing in match, seg: "<<seg_id;
      return seg_id;
    }

    double diff = std::abs(moving_count - fixed_count);
    double total = moving_count + fixed_count;
    double overlap = diff / total;

    if (overlap > m_segment_size_match_threshold)
    {
      return seg_id;
    }
  }

  return ret_state;
}

////////////////////////////////////////////////////////////////////////////////

void MeshTracking::ExportSegments(
  const VidSidMap &vs_map,
  const pcl::PolygonMesh &mesh,
  const std::string &out_path) const
{

  CloudProcessing cp;
  MeshProcessing mp;
  pcl::PolygonMesh tmp(mesh);
  CloudPtr in_cloud = mp.GetCloudFromPolygonMesh(tmp);
  for (int seg_id = -1; seg_id < this->curr_num_segs_; seg_id++)
  {
    // MeshSegmentation ms;
    // pcl::PolygonMesh out_mesh = ms.MeshSegFromVertexIds(mesh,vs_map,seg_id);
    CloudPtr tmp(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
    for (VidSidMap::const_iterator it = vs_map.begin(); it != vs_map.end(); ++it)
    {
      if (it->second == seg_id)
      {
        tmp->push_back(in_cloud->points[it->first]);
      }
    }
    
    std::array<int,3> rgb = {
      utils::RandInt(0,255),
      utils::RandInt(0,255),
      utils::RandInt(0,255)};

    if (seg_id == -1)
    {
      rgb = {255,255,255};
    }

    // CloudPtr tmp = mp.GetCloudFromPolygonMesh(out_mesh);
    cp.ColorizeCloud(tmp,rgb);
    // pcl::toPCLPointCloud2(*tmp, out_mesh.cloud);

    std::string filename = out_path + "seg_" + std::to_string(seg_id) + ".ply";
    int res = pcl::io::savePLYFile(filename, *tmp);// out_mesh);
  }
}

////////////////////////////////////////////////////////////////////////////////

CloudPtr MeshTracking::AlignSegment(
  CloudPtr& moving_cloud,
  CloudPtr& fixed_cloud,
  const VidSidMap& moving_vs_map,
  const VidSidMap& fixed_vs_map,
  const int& seg_id)
{

  MeshProcessing mp;
  CloudProcessing cp;
  MeshSegmentation ms;
  Matcher mt;

  // if (!cp.NormalsExist(fixed_cloud) || !cp.NormalsExist(moving_cloud))
  // {
  //   BOOST_LOG_TRIVIAL(error) << "In MeshTracking::AlignSegment: normals could not be calculated for input clouds";
  // }

  // copy input
  CloudPtr out_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>(*moving_cloud));

  // populate moving segment group with vertex indices
  std::vector<int> moving_indices = ms.GetSegIndices(moving_vs_map,seg_id);

  // populate fixed segment group with vertex indices
  std::vector<int> fixed_indices = ms.GetSegIndices(fixed_vs_map,seg_id);

  if (moving_indices.size() == 0 || fixed_indices.size() == 0)
  {
    std::string e = 
      "In AlignSegment: no indices for seg_:"+std::to_string(seg_id)+"...";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }

  CloudPtr moving_seg(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  CloudPtr fixed_seg(new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  // populate new seg_clouds with vertex info from segment.indices
  for (size_t idx = 0; idx < moving_indices.size(); idx++)
  {
    moving_seg->push_back(moving_cloud->points[moving_indices[idx]]);
  }
  for (size_t idx = 0; idx < fixed_indices.size(); idx++)
  {
    fixed_seg->push_back(fixed_cloud->points[fixed_indices[idx]]);
  }

  // Perform registration
  CloudPtr aligned_points(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  CloudPtr prior_points(new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  // Use a rigid alignment if fused
  std::vector<int>::iterator it = std::find( 
  this->fused_segments_.begin(),  this->fused_segments_.end(), seg_id);

  if (it != this->fused_segments_.end())
  {
    BOOST_LOG_TRIVIAL(debug) << "using rigid ICP for seg: "<<seg_id;
    aligned_points = cp.RigidICP(
      moving_seg,
      fixed_seg,
      0.05, // max dist match
      50, // max_iters
      1e-8, // tform_epsilon
      5e-3); // fitness_epsilon
  }
  else
  { 
    CPDICP cpd;
    CPD_Result cpd_result = cpd.AlignHighResClouds(
      fixed_seg, 
      moving_seg, 
      this->default_CPD_tolerance_);

    aligned_points = cpd_result.aligned_cloud;
  }

  // store new positions in output cloud
  for (size_t idx = 0; idx < aligned_points->size(); idx++)
  {
    out_cloud->points[moving_indices[idx]] = aligned_points->points[idx];
  }

  return out_cloud; 
}

////////////////////////////////////////////////////////////////////////////////

CtrSidMap MeshTracking::GetSegmentCentroids(const VidSidMap &vs_map,
                                            const pcl::PolygonMesh &mesh) const
{

  // Strip ptcloud from mesh
  CloudPtr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  pcl::fromPCLPointCloud2(mesh.cloud, *cloud);

  CtrSidMap cs_map;

  for (size_t seg_id = 0; seg_id < this->curr_num_segs_; seg_id++)
  {

    pcl::CentroidPoint<pcl::PointXYZRGBNormal> centroid_calculator;

    // populate the pcl::centroid class with points belonging to current segment
    for (VidSidMap::const_iterator it = vs_map.begin(); it != vs_map.end(); ++it)
    {
      if (it->second == seg_id)
      {
        centroid_calculator.add(cloud->points[it->first]);
      }
    }

    pcl::PointXYZRGBNormal centroid;
    centroid_calculator.get(centroid);

    cs_map[seg_id] = centroid;
  }

  return cs_map;
}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshTracking::AlignRigid(
  pcl::PolygonMesh& moving,
  pcl::PolygonMesh& fixed)
{

  CloudPtr m_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  CloudPtr f_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  pcl::fromPCLPointCloud2(moving.cloud,*m_cloud);
  pcl::fromPCLPointCloud2(fixed.cloud,*f_cloud);

  CloudProcessing cp;
  CloudPtr a_cloud = cp.AlignCentres(m_cloud,f_cloud);
  // CloudPtr a_cloud = cp.RigidICP(m_cloud,f_cloud);
  pcl::PolygonMesh outmesh(moving);
  pcl::toPCLPointCloud2(*a_cloud,outmesh.cloud);
  return outmesh;
}

////////////////////////////////////////////////////////////////////////////////

// Modified implementation of VTKMeshTracker. Performs a surface-aware iterative
// non-rigid ICP
pcl::PolygonMesh MeshTracking::TrackMesh
(
  const pcl::PolygonMesh & _moving_pcl_mesh,
  const pcl::PolygonMesh &   _fixed_pcl_mesh,
  const float            &     _alpha_max, // default = 15
  const float            &     _alpha_min, // deafult = 5
  const int              &      _iter_num, // default = 5
  const float            &         _gamma  // default = 1.0
) const
{  
  if (_alpha_min <= 0 || _alpha_min > _alpha_max)
  {
    std::string e = "In TrackMesh: Provided alphas are unviable";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }
  
  CloudProcessing cp;

  vtkSmartPointer<vtkPolyData>    moving_vtk;
  vtkSmartPointer<vtkPolyData>    fixed_vtk;

  pcl::VTKUtils::convertToVTK(_moving_pcl_mesh, moving_vtk);
  pcl::VTKUtils::convertToVTK(  _fixed_pcl_mesh,   fixed_vtk);

  // Cloud for correspondences
  CloudPtr fixed_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  pcl::fromPCLPointCloud2(_fixed_pcl_mesh.cloud,*fixed_cloud);

  // Tmp cloud use to maintain normals and calc init error
  // TODO: recalculate normals from surface mesh
  CloudPtr moving_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  pcl::fromPCLPointCloud2(_moving_pcl_mesh.cloud,*moving_cloud);

  // initialises edges, vertices and KNN search (unused if matches supplied)
  NonRigidICP nri(moving_vtk, fixed_vtk);
  nri.init();

  // linear interpolation between max and min for _iter_num steps
  std::vector<float> linspaced_alphas;
  const float delta = (_alpha_min - _alpha_max) / (_iter_num - 1);
  for (int iter_idx = 0; iter_idx < _iter_num - 1; ++iter_idx) {
      linspaced_alphas.push_back(_alpha_max + delta * static_cast<float>(iter_idx));
  }
  
  // Ensure that _alpha_max and _alpha_min are exactly the same as the input
  linspaced_alphas.push_back(_alpha_min);

  pcl::PolygonMesh result_mesh;
  unsigned int iter = 0;
  double error_prev = cp.ComputeHausdorffDistance(moving_cloud,fixed_cloud);

  pcl::VTKUtils::convertToPCL(moving_vtk, result_mesh);
  CloudPtr res_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  pcl::fromPCLPointCloud2(result_mesh.cloud,*res_cloud);

  // Linear interpolate for given steps between alpha range
  for (const float & alpha : linspaced_alphas) {

    BOOST_LOG_TRIVIAL(debug) << "Icp iteration: "<<iter<<"/"<<linspaced_alphas.size();
    iter++;

    // Perform initial Rigid ICP alignment to compensate for bad alpha
    // in the NonRigidICP code
    // BOOST_LOG_TRIVIAL(debug) << "Rigid ICP";

    // Tried beginning with aligned centres to see if it would avoid
    // local minima.... it did not
    // res_cloud = cp.AlignCentres(moving_cloud,fixed_cloud);
    // res_cloud = cp.RigidICP(moving_cloud,fixed_cloud);

    // pcl/VTK conversion deletes normals. Seeing as it's a rigid transform we should be 
    // safe to just iterate and restore normals in a pointwise manner
    // TODO: assert that tmp has normals
    for (size_t pt = 0; pt < res_cloud->size(); pt++)
    {
      res_cloud->points[pt].normal_x = moving_cloud->points[pt].normal_x;
      res_cloud->points[pt].normal_y = moving_cloud->points[pt].normal_y;
      res_cloud->points[pt].normal_z = moving_cloud->points[pt].normal_z;
    }
    
    // perform conversion and correspondence match updates, fyi this function
    // will calculate normals for the input cloud if none exist.
    BOOST_LOG_TRIVIAL(debug) << "correspondance match"; 
    Matcher mt; 
    std::vector<match_conf> matches = mt.GetMatchesAndConfidence(
      res_cloud,
      fixed_cloud,
      0.25,
      5.0);

    std::vector<std::array<double,3>> mt_matches;
    for (size_t match = 0; match < matches.size(); match++)
    {
      std::array<double,3> point;
      int idx = matches[match].first;
      point[0] = fixed_cloud->points[idx].x;
      point[1] = fixed_cloud->points[idx].y;
      point[2] = fixed_cloud->points[idx].z;
      mt_matches.push_back(point);
    }

    BOOST_LOG_TRIVIAL(debug) << "NonRigid ICP";
    nri.initCompute(mt_matches);
    nri.compute(alpha, _gamma);

    // Check error difference for early convergence
    moving_vtk = nri.GetTemplate();
    pcl::VTKUtils::convertToPCL(moving_vtk, result_mesh);
    pcl::fromPCLPointCloud2(result_mesh.cloud,*res_cloud);

    double error = cp.ComputeHausdorffDistance(res_cloud,fixed_cloud);
    double error_diff = std::abs(error - error_prev);
    if (error_diff < tracking_error_thresh_*(1+error_prev))   
    {
      BOOST_LOG_TRIVIAL(debug) << "Tracking converged early!";
      break;
    }
    error_prev = error;

  }
  return result_mesh;
}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshTracking::DeformMesh(
    const pcl::PolygonMesh& in_mesh,
    const pcl::PolygonMesh& d_graph,
    const std::vector<std::vector<int>>& node_influence_list)
{

  CloudPtr in_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  CloudPtr d_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  pcl::fromPCLPointCloud2(in_mesh.cloud, *in_cloud);
  pcl::fromPCLPointCloud2(d_graph.cloud, *d_cloud);

  // output
  CloudPtr nu_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>(*in_cloud));

	for (size_t pt=0; pt<in_cloud->size(); ++pt) 
  {
    std::vector<double> weights;
    weights.resize(kDG_KNN);

    Eigen::Vector3d curr_pt(
      in_cloud->points[pt].x,
      in_cloud->points[pt].y,
      in_cloud->points[pt].z); 

    // get sqr distances to nn_nodes
    std::vector<double> dists_node;
    for (size_t node = 0; node < kDG_KNN; node++)
    {
        int nid = node_influence_list[pt][node];

        Eigen::Vector3d d_pt(
          d_cloud->points[nid].x,
          d_cloud->points[nid].y,
          d_cloud->points[nid].z);
        double x_sqr = (curr_pt[0]-d_pt[0])*(curr_pt[0]-d_pt[0]);
        double y_sqr = (curr_pt[1]-d_pt[1])*(curr_pt[1]-d_pt[1]);
        double z_sqr = (curr_pt[2]-d_pt[2])*(curr_pt[2]-d_pt[2]);
        dists_node.push_back(x_sqr + y_sqr + z_sqr);
    }
      
    // find furthest
    auto it = max_element(std::begin(dists_node), std::end(dists_node));
    double d_max = sqrt(dists_node[it - dists_node.begin()]); 
    double weight_sum = 0;

    for (int j=0; j<kDG_KNN; ++j) 
    {
        // weights[j] = 0.9 - j; 
        weights[j] = (1-sqrt(dists_node[j])/ d_max) ;
        weights[j] = weights[j]*weights[j];
        weight_sum += weights[j];
    }

    for (int j=0; j<kDG_KNN; ++j)
    {
        weights[j] /= weight_sum;
    }


    Eigen::Vector3d nu_pt(0.0, 0.0, 0.0);
    Eigen::Matrix3d R;
    R << 1,0,0,0,1,0,0,0,1;

    for (int n=0; n<kDG_KNN; ++n) 
    {
      int dn = node_influence_list[pt][n];

      Eigen::Vector3d d_pt(
        d_cloud->points[dn].x,
        d_cloud->points[dn].y,
        d_cloud->points[dn].z);

      nu_pt += weights[n]*(R*(curr_pt-d_pt)+d_pt+(d_pt-curr_pt));
    }

    nu_cloud->points[pt].x =nu_pt[0];
    nu_cloud->points[pt].y =nu_pt[1];
    nu_cloud->points[pt].z =nu_pt[2];
  }

  pcl::PolygonMesh out_mesh(in_mesh);
  pcl::toPCLPointCloud2(*nu_cloud,out_mesh.cloud);

  return out_mesh;

}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::array<double,6>> MeshTracking::GetTform
(
  const pcl::PolygonMesh & source,
  const pcl::PolygonMesh & target
)
{

  // check that inputs have same num vertices
  if (source.cloud.width != target.cloud.width)
  {
    std::string e = "In MeshTracking::GetTform: Inputs don't match";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }

  // extract clouds
  MeshProcessing mp;
  pcl::PolygonMesh tmp_s = source;
  pcl::PolygonMesh tmp_t = target;
  CloudPtr s_cloud = mp.GetCloudFromPolygonMesh(tmp_s);
  CloudPtr t_cloud = mp.GetCloudFromPolygonMesh(tmp_t);

  // Get tform mapping src to tgt
  std::vector<std::array<double, 6>> tforms;

  for (size_t pt = 0; pt < s_cloud->size(); pt++)
  {
    std::array<double,6> tform;
    tform[0] = (t_cloud->points[pt].x - s_cloud->points[pt].x);
    tform[1] = (t_cloud->points[pt].y - s_cloud->points[pt].y);
    tform[2] = (t_cloud->points[pt].z - s_cloud->points[pt].z);

    tform[3] = (t_cloud->points[pt].normal_x - s_cloud->points[pt].normal_x);
    tform[4] = (t_cloud->points[pt].normal_y - s_cloud->points[pt].normal_y);
    tform[5] = (t_cloud->points[pt].normal_z - s_cloud->points[pt].normal_z);
  
    tforms.push_back(tform);
  }
  return tforms;
}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshTracking::ApplyTform
(
  const pcl::PolygonMesh                  & in_mesh,
  const std::vector<std::array<double,6>> &   tform,
  const bool                              &  invert
)
{
  MeshProcessing mp;
  pcl::PolygonMesh tmp_mesh(in_mesh);
  CloudPtr cloud = mp.GetCloudFromPolygonMesh(tmp_mesh);

  // check that inputs have same num vertices
  if (cloud->size() != tform.size())
  {
    std::string e = "In MeshTracking::ApplyTform: Inputs don't match";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }

  std::vector<std::array<double,6>> tform_(tform);

  if(invert)
  {for (size_t pt = 0; pt < tform.size(); pt++)
  {
    for (size_t i = 0; i < 6; i++)
    {
      tform_[pt][i] *= -1;
    }
  }
  }

  // apply tform
  for (size_t pt = 0; pt < cloud->size(); pt++)
  {

    cloud->points[pt].x += tform_[pt][0];
    cloud->points[pt].y += tform_[pt][1];
    cloud->points[pt].z += tform_[pt][2];
    cloud->points[pt].normal_x += tform_[pt][3];
    cloud->points[pt].normal_y += tform_[pt][4];
    cloud->points[pt].normal_z += tform_[pt][5];
  }

  // attach to output and return
  pcl::toPCLPointCloud2(*cloud, tmp_mesh.cloud);
  return(tmp_mesh);   

}

////////////////////////////////////////////////////////////////////////////////

bool MeshTracking::StoreTform
(
  const pcl::PolygonMesh  & source,
  const pcl::PolygonMesh  & target,
  const std::string       &   path
) const 
{

  // check that inputs have same num vertices
  if (source.cloud.width != target.cloud.width)
  {
    std::string e = "In MeshTracking::StoreTform: Inputs don't match";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }

  // extract clouds
  MeshProcessing mp;
  pcl::PolygonMesh tmp_s = source;
  pcl::PolygonMesh tmp_t = target;
  CloudPtr s_cloud = mp.GetCloudFromPolygonMesh(tmp_s);
  CloudPtr t_cloud = mp.GetCloudFromPolygonMesh(tmp_t);

  // Get tform mapping src to tgt
  std::vector<std::array<double, 6>> tforms;

  // start output
  std::ofstream ofs(path);

  for (size_t pt = 0; pt < s_cloud->size(); pt++)
  {
    ofs << (t_cloud->points[pt].x - s_cloud->points[pt].x) << "\n";
    ofs << (t_cloud->points[pt].y - s_cloud->points[pt].y) << "\n";
    ofs << (t_cloud->points[pt].z - s_cloud->points[pt].z) << "\n";

    ofs << (t_cloud->points[pt].normal_x - s_cloud->points[pt].normal_x) << "\n";
    ofs << (t_cloud->points[pt].normal_y - s_cloud->points[pt].normal_y) << "\n";
    ofs << (t_cloud->points[pt].normal_z - s_cloud->points[pt].normal_z) << "\n";
  }
  
  ofs.close();

  return true;

}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshTracking::GetDetailLayer
(
  const pcl::PolygonMesh  &         source,
  const pcl::PolygonMesh  &         target,
  const bool              & strict_matches
)
{


  // Get point to surface correspondences 

  // Build NonLin solver

  // Solve (apply vanishing weight)

  // Return 

  // or......be lazy
  // super non-rigid surface alignment
  //----------------------------------------------------------------------------
  vtkSmartPointer<vtkPolyData>    moving_vtk;
  vtkSmartPointer<vtkPolyData>    fixed_vtk;
  pcl::VTKUtils::convertToVTK(source, moving_vtk);
  pcl::VTKUtils::convertToVTK(target, fixed_vtk);

  // Cloud for correspondences
  MeshProcessing mp;
  pcl::PolygonMesh tmp_target(target);
  CloudPtr fixed_cloud = mp.GetCloudFromPolygonMesh(tmp_target);

  // Tmp cloud use to maintain normals and calc init error
  // TODO: recalculate normals from surface mesh
  pcl::PolygonMesh tmp_source(source);
  CloudPtr moving_cloud = mp.GetCloudFromPolygonMesh(tmp_source);

  NonRigidICP nri(moving_vtk, fixed_vtk);
  nri.init();

  // TODO: Refine a bit, some weird results still possible from bad matches
  float _alpha_min = 5.00; // floss:5
  float _alpha_max = 15.00; // floss:15
  float _gamma = 0.5; // 0.5
  int _iter_num = 5; // floss: 3

  std::vector<float> linspaced_alphas;
  const float delta = (_alpha_min - _alpha_max) / (_iter_num - 1);
  for (int iter_idx = 0; iter_idx < _iter_num - 1; ++iter_idx) {
      linspaced_alphas.push_back(_alpha_max + delta * static_cast<float>(iter_idx));
  }

  pcl::PolygonMesh result_mesh;

  // Ensure that _alpha_max and _alpha_min are exactly the same as the input
  linspaced_alphas.push_back(_alpha_min);

  unsigned int iter = 0;

  pcl::VTKUtils::convertToPCL(moving_vtk, result_mesh);
  CloudPtr res_cloud = mp.GetCloudFromPolygonMesh(result_mesh);

  // Linear interpolate for given steps between alpha range
  for (const float & alpha : linspaced_alphas) 
  {
    BOOST_LOG_TRIVIAL(debug) << "Icp iteration: "<<iter<<"/"<<linspaced_alphas.size();
    iter++;

    // pcl/VTK conversion deletes normals. Seeing as it's a rigid transform we should be 
    // safe to just iterate and restore normals in a pointwise manner
    // TODO: assert that tmp has normals
    for (size_t pt = 0; pt < res_cloud->size(); pt++)
    {
      res_cloud->points[pt].normal_x = moving_cloud->points[pt].normal_x;
      res_cloud->points[pt].normal_y = moving_cloud->points[pt].normal_y;
      res_cloud->points[pt].normal_z = moving_cloud->points[pt].normal_z;
    }

    // perform conversion and correspondence match updates, fyi this function
    // will calculate normals for the input cloud if none exist.
    BOOST_LOG_TRIVIAL(debug) << "correspondance match"; 
    Matcher mt; 
    std::vector<int> matches = mt.GetConstrainedMatches(
      moving_cloud,fixed_cloud,0.25,0.25);

    std::vector<std::array<double,3>> mt_matches;
    for (size_t match = 0; match < matches.size(); match++)
    {
      std::array<double,3> point;
      int idx = matches[match];

      // if -1 don't move the point.....test scenario
      if(matches[match]==-1)
      {
        point[0] = moving_cloud->points[match].x;
        point[1] = moving_cloud->points[match].y;
        point[2] = moving_cloud->points[match].z;
      }else
      {
        point[0] = fixed_cloud->points[idx].x;
        point[1] = fixed_cloud->points[idx].y;
        point[2] = fixed_cloud->points[idx].z;
      }
      mt_matches.push_back(point);
    }

    BOOST_LOG_TRIVIAL(debug) << "NonRigid ICP, alpha: "<<alpha;
    strict_matches ? nri.initCompute(mt_matches) : nri.initCompute();
    nri.compute(alpha, _gamma);

    moving_vtk = nri.GetTemplate();
    pcl::VTKUtils::convertToPCL(moving_vtk, result_mesh);
  }

  return result_mesh;
  //----------------------------------------------------------------------------

}

////////////////////////////////////////////////////////////////////////////////

NodeCorrs MeshTracking::GetCorresConstraints
(
  const DeformGraph & graph,
  const CloudPtr    & corr_aligned_cloud,
  const CloudPtr    & corr_input_cloud,
  const CloudPtr    & fixed_cloud
)
{

  NodeCorrs corres_constraints;
  CloudProcessing cp;

  if (graph.node_pos.size()==0) 
  {
    std::string e = 
      "In GetCorresConstraints: Deformation graph not properly initialized";
      BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }

  for (size_t pt = 0; pt < corr_input_cloud->size(); pt++)  
  {

      pcl::PointXYZRGBNormal src, dst;
      src = corr_input_cloud->points[pt];
      dst = corr_aligned_cloud->points[pt];
      
      corr_aligned_cloud->points[pt].r = 255;
      corr_aligned_cloud->points[pt].g = 0;
      corr_aligned_cloud->points[pt].b = 0;
      corr_input_cloud->points[pt].r = 0;
      corr_input_cloud->points[pt].g = 255;
      corr_input_cloud->points[pt].b = 0;

      std::pair<pcl::PointXYZRGBNormal,pcl::PointXYZRGBNormal> corr = 
        std::make_pair(src, dst);
      
      corres_constraints.push_back(corr); // match pos [xyz]
    // }

  }

  return corres_constraints;
}
////////////////////////////////////////////////////////////////////////////////