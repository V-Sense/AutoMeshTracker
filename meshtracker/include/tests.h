#ifndef _TESTING_H
#define _TESTING_H

#include <algorithm>

#include "mesh_tracking.h"
#include "kalman_smoothing.h"
#include "mesh_processing.h"
#include "cloud_processing.h"
#include "log.h"

#include "keyframer.h"

namespace TESTING
{
////////////////////////////////////////////////////////////////////////////////
// Fundamental tracking test for pairwise alignment

    static void OneShotTracking
    (
        const std::string & f_moving,
        const std::string &  f_fixed,
        const std::string & f_result,
        std::string       debug_path = "/MeshTracker/resource/tests/oneshot_test/debug/",
        bool              save_debug = true
    )
    {
        // start timer
        time_t tic, toc; 
        tic = time(0);

        MeshTracking mt;
        
        MeshProcessing mp;
        CloudProcessing cp;
        pcl::PolygonMesh moving, fixed, result;

        // setup debug output
        mt.SaveDebugOutput(save_debug,debug_path);

        BOOST_LOG_TRIVIAL(debug) << "loading input...";
        moving = mt.LoadMeshFromFile(f_moving);
        fixed = mt.LoadMeshFromFile(f_fixed);

        // Clean up input 
        moving = mp.RemoveSmallComponents(moving,0.01);
        fixed = mp.RemoveSmallComponents(fixed,0.01);
        moving = mp.FixManifoldness(moving);
        fixed = mp. FixManifoldness(fixed);

        result = mt.DoAdaptiveTracking(moving,fixed);

        if(pcl::io::saveOBJFile(f_result, result)<0)
        {
            std::string e = "\tIn 1Shot Test: Couldn't save output";
            BOOST_LOG_TRIVIAL(error) << e;
            throw std::runtime_error(e);
        }

        toc = time(0); 
        cout << "Elapsed: "<< difftime(toc, tic) <<" second(s)."<< endl;

        // Get error
        CloudPtr res_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
        CloudPtr gt_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
        pcl::fromPCLPointCloud2(result.cloud,*res_cloud);
        pcl::fromPCLPointCloud2(fixed.cloud,*gt_cloud);
        const double err = cp.ComputeHausdorffDistance(res_cloud,gt_cloud);
        BOOST_LOG_TRIVIAL(info) << "Error: "<<err;

        // Save error
        const std::string f_err = f_result.substr(0,f_result.size()-4) + "_error.txt";
        ofstream errfile(f_err);
        if (errfile.is_open())
        {
            errfile << std::to_string(err);
            errfile.close();
        }
        else cout << "Unable to open file: "<<f_err;

    }

////////////////////////////////////////////////////////////////////////////////



    static void OneShotDebug()
    {
        std::string f_log = "/MeshTracker/resource/tests/oneshot_test/debug/output.log";
        init_logging(0,f_log); //log everything
        BOOST_LOG_TRIVIAL(debug) << "Running test: OneShotDebug\n";

        std::string f_moving = "/MeshTracker/resource/tests/oneshot_test/debug/moving.obj";
        std::string f_fixed = "/MeshTracker/resource/tests/oneshot_test/debug/fixed.obj";
        std::string f_result = "/MeshTracker/resource/tests/oneshot_test/debug/result.obj";

        OneShotTracking(f_moving,f_fixed,f_result);
    }


 ///////////////////////////////////////////////////////////////////////////////
 // Error outputs

    static std::vector<double> GetSequenceError(
        const std::string meshlist,
        const std::string res_path)
    {
        CloudProcessing cp;
        std::vector<double> errors; 

        std::string line;
        ifstream myfile (meshlist);
        if (myfile.is_open())
        {
            while ( getline (myfile,line) )
            {
                pcl::PolygonMesh res_mesh, gt_mesh;
                int path_end = line.find_last_of("/\\")+1;
                std::string fend = line.substr(path_end,line.size()-path_end-4);
                std::string res_fname = res_path + fend + "_t.obj";

                // load result mesh
                pcl::io::loadPolygonFile(res_fname, res_mesh);

                // load target as ground truth
                pcl::io::loadPolygonFile(line, gt_mesh);

                CloudPtr res_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
                CloudPtr gt_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
                pcl::fromPCLPointCloud2(res_mesh.cloud,*res_cloud);
                pcl::fromPCLPointCloud2(gt_mesh.cloud,*gt_cloud);

                errors.push_back(cp.ComputeHausdorffDistance(gt_cloud,res_cloud));
            }
            myfile.close();


        } else
        {
            cout << "Unable to open file";
        } 
        return errors;
    }

 ///////////////////////////////////////////////////////////////////////////////

    static void SaveSequenceError(
        const std::vector<double> errors, 
        const std::string& fname
    )
    {
        // max error
        double max_e = *std::max_element(errors.begin(), errors.end());

        // median error
        double median_e = 0;

        std::vector<double> sorted_errors(errors);
        std::sort(sorted_errors.begin(),sorted_errors.end()); 
    
        // check for even case 
        if (errors.size() % 2 != 0) 
        {
            median_e = (double)sorted_errors[errors.size()/2]; 
        } else
        {
            median_e = 
                (double)(sorted_errors[(errors.size()-1)/2] + 
                        sorted_errors[errors.size()/2])/2.0; 
        }
        
        ofstream myfile(fname);
        if (myfile.is_open())
        {
            myfile << std::to_string(max_e) << " " ;
            myfile << std::to_string(median_e) << "\n";
            for (size_t e = 0; e < errors.size(); e++)
            {
                myfile << std::to_string(errors[e]) << "\n"; 
            }
            myfile.close();
        }
        else cout << "Unable to open file: "<<fname;
        
    }

///////////////////////////////////////////////////////////////////////////////

    static void PlotSequenceError(
        const std::string& error_fname
    )
    {   std::string cmd = "python3 /MeshTracker/resource/scripts/plot_error.py "+error_fname;
        system(cmd.c_str());
    }

 ///////////////////////////////////////////////////////////////////////////////
 // Extends pairwise tracking to test performance across a sequence 

    static void SequenceTesting(
        const std::string & path
    )
    {
        // start timer
        time_t tic, toc; 
        tic = time(0);

        MeshTracking mt;    
        KeyFramer kf;

        // create result directory
        std::string res_path = path + "result/";
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
        BOOST_LOG_TRIVIAL(debug) << "Running test from resource...\n";

        #if false
        {
        // Create meshlist
        std::string cmd = "ls -d $PWD/"+path+"input/*.ply > "+path+"input/meshlist.txt"; 
        system(cmd.c_str());

        // Read meshlist and store
        mt.ReadMeshList(path+"input/meshlist.txt");

        // Generate Abs layers from meshlist and save in result folder
        mt.GenAbslayersFromList(path+"result/");

        // Generate Spherical Harmonics descriptor for each abs layers using the 
        // third party exe
        std::string msh2shd = "/MeshTracker/resource/thirdparty/msh2shd";
        kf.GenerateSPHDescriptors(path+"result/abs/meshlist.txt",msh2shd);

        // Generate feasibility scores from abs layers
        // kf.GenerateFeasScore(kf.ReadMeshList(path+"result/abs/meshlist.txt"),path);
        kf.GenerateFeasScore(kf.ReadMeshList(path+"result/abs/meshlist.txt"),path);

        // Execute matlab script to find and plot eligible keyframes and tracking
        // regions
        std::string descriptors_p = path+"result/abs/descriptors.csv";
        std::string feas_p = path+"result/abs/feas.txt";
        cmd = "sh /MeshTracker/resource/scripts/gen_regions.sh "+descriptors_p+" "+feas_p;
        system(cmd.c_str());
        exit(1);
        }
        #endif

        // set variables and pass object to mt
        icp::IcpOptions icp_ops;
        icp_ops.use_adaptive = true;
        icp_ops.parallel = true;
        std::string regions_file = path+"input/regions.txt";
        std::string tracking_out = path+"result/";
       
        // perform tracking and detail synth
        mt.ReadMeshList(path+"input/meshlist.txt");
        mt.ReadRegions(regions_file);
        mt.TrackMeshSequence(icp_ops,tracking_out);
        
        // Kalman smoothing
        // Create meshlist
        // std::string cmd = "ls -d $PWD/"+path+"result/*_t.obj > "+path+"result/meshlist.txt"; 
        // system(cmd.c_str());
        // KalmanSmoother ks;
        // ks.readMeshes(path+"result/smooth/meshlist.txt");
        // ks.ReadRegions(regions_file);
        // ks.smoothInRegions();
        // std::cerr<<"done smoothing, saving"<<"\n";
        // ks.setupSmoothMeshes(true);

        toc = time(0); 
        cout << "Elapsed: "<< difftime(toc, tic) <<" second(s)."<< endl;

        // Calculate and export hdorf metrics 
        // std::vector<double> errors = GetSequenceError(path+"input/meshlist.txt",tracking_out);
        // std::string err_path = tracking_out + "error.txt";
        // SaveSequenceError(errors,err_path);
    }

    static void SequenceTestingMattHello()
    {
        std::string path = "/MeshTracker/resource/tests/sequence_test/matt_hello/";
        SequenceTesting(path);
    }

    static void SequenceTestingFloss()
    {
        std::string path = "/MeshTracker/resource/tests/sequence_test/floss/";
        SequenceTesting(path);
    }

} // end namespace
#endif