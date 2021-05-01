
#include <iostream>
#include <string>
#include <ctime>
#include <sys/stat.h> 
#include <sys/types.h> 
#include <bits/stdc++.h> 

#include "utils.h"
#include "constants.h"
#include "keyframer.h"
#include "mesh_tracking.h"
#include "kalman_smoothing.h"
#include "mesh_segmentation.h"
#include "prog_opts.h"
#include "log.h"
#include "tests.h"

#include <optional>
#include <boost/program_options.hpp>

#include <pcl/pcl_config.h>

int main(int argc, char *argv[])
{
  std::cout << "Cxx: " << __cplusplus << std::endl;
  std::cout<<"Eigen version: "
  <<std::to_string(EIGEN_WORLD_VERSION)<<"."
  <<std::to_string(EIGEN_MAJOR_VERSION)<<"."
  <<std::to_string(EIGEN_MINOR_VERSION)<<"."<<std::endl;
  std::cout<<"PCL version: "<< PCL_VERSION << std::endl;

  // Disable all warnings. Otherwise pcl displays a list of missing fields
  // when loading each ply file.
  pcl::console::setVerbosityLevel(pcl::console::L_ALWAYS);

  bool use_adapt_tracking, parallel,testing;
  uint verbosity_level;

  try
  {

    namespace prog_opts = boost::program_options;

    prog_opts::options_description desc("Options");
    desc.add_options()
      ("help,-h", "Print help message")
      ("testing,-t", "Perform tests in resource folder.")
      ("run,-r", "Run tracking from '/input' folder" )
      ("verbosity,-v",prog_opts::value<uint>(&verbosity_level)->default_value(2),"set verbosity level: (5) fatal, (4) error, (3) warning, (2) <DEFAULT> info, (1) debug, (0) trace");

    prog_opts::variables_map var_map;
    prog_opts::store(prog_opts::parse_command_line(argc, argv, desc), var_map);

    if (var_map.count("help") || argc == 1)
    {
      std::cerr << desc << '\n';
      return 0;
    }

    if (var_map.count("testing"))
    {

      try
      {
        TESTING::SequenceTestingMattHello();
        // TESTING::SequenceTestingFloss();
      }
      catch(const std::exception& e)
      {
        std::cerr << e.what() << '\n';
        return 0;
      }

      return 0;
    }
    prog_opts::notify(var_map);
  }

  catch (std::exception &e)
  {
    const std::string err = e.what();
    BOOST_LOG_TRIVIAL(error) << err;
    std::cerr << "Run: " << argv[0] << " -h   for usage" << std::endl;
    return 1;
  }

  catch (...)
  {
    std::string e = "Unknown Exception encountered!";
    BOOST_LOG_TRIVIAL(error) << e;
    std::cerr << e << std::endl;
  }

  // Mesh tracking from user input
  try
  {

    // start timer
    time_t tic, toc; 
    tic = time(0);

    MeshTracking mt;    
    KeyFramer kf;

    // create result directory
    std::string res_path = "/result/";
    if (mkdir(res_path.c_str(), 0777) == -1)
    {
        std::cerr
        << "Couldn't create result directory at:\n" 
        << res_path
        << "\n path may already exist..."
        << std::endl;
    }

    // Setup logs
    std::string f_log = res_path + "output.log";
    init_logging(0,f_log); //log everything
    BOOST_LOG_TRIVIAL(debug) << "Tracking from user input...\n";

    // Create meshlist
    // std::string cmd = "ls /input/*{.obj,.ply} > /input/meshlist.txt"; 
    std::string cmd = "ls /input/*.ply /input/*.obj > /input/meshlist.txt"; 
    system(cmd.c_str());
    // Read meshlist and store
    mt.ReadMeshList("/input/meshlist.txt");

    // Generate Abs layers from meshlist and save in result folder
    // mt.GenAbslayersFromList("/result/");

    // Generate Spherical Harmonics descriptor for each abs layers using the 
    // third party exe
    // std::string msh2shd = "/MeshTracker/resource/thirdparty/msh2shd";
    // kf.GenerateSPHDescriptors("/result/abs/meshlist.txt",msh2shd);

    // Generate feasibility scores from abs layers
    // kf.GenerateFeasScore(kf.ReadMeshList("/result/abs/meshlist.txt"),path);
    // kf.GenerateFeasScore(kf.ReadMeshList("/result/abs/meshlist.txt"),"");

    // Execute matlab script to find and plot eligible keyframes and tracking
    // regions
    // std::string descriptors_p = "/result/abs/descriptors.csv";
    // std::string feas_p = "/result/abs/feas.txt";
    // cmd = "sh /MeshTracker/resource/scripts/gen_regions.sh "+descriptors_p+" "+feas_p;
    // system(cmd.c_str());
 
    // set variables and pass object to mt
    icp::IcpOptions icp_ops;
    icp_ops.use_adaptive = true;
    icp_ops.parallel = true;
    std::string regions_file = "/input/regions.txt";
    std::string tracking_out = "/result/";
    
    // perform tracking
    mt.ReadMeshList("/input/meshlist.txt");
    mt.ReadRegions(regions_file);
    mt.TrackMeshSequence(icp_ops,tracking_out);

    // Create meshlist for smoothing
    cmd = "ls /result/*.obj > /result/smooth/meshlist.txt"; 
    system(cmd.c_str());
    KalmanSmoother ks;
    ks.readMeshes("/result/smooth/meshlist.txt");
    ks.ReadRegions(regions_file);
    ks.smoothInRegions();
    ks.setupSmoothMeshes(true);

    toc = time(0); 
    cout << "Elapsed: "<< difftime(toc, tic) <<" second(s)."<< endl;

    // Calculate and export hdorf metrics 
    std::vector<double> errors = TESTING::GetSequenceError("/input/meshlist.txt",tracking_out);
    std::string err_path = tracking_out + "error.txt";
    TESTING::SaveSequenceError(errors,err_path);
  }
  catch (std::exception& e)
  {
    std::string dump =
      std::string(e.what()) + "\n\tCheck output.log in project folder for details.\n";
    std::cerr<<dump;
  }

  catch(...)
  {
    std::cerr<<'Unkown exception!\n\tCheck output.log in project folder for details.\n';
  }

  return 0;
}
