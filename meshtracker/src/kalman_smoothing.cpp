#include <iostream>
#include <fstream>
#include <algorithm>

#include <pcl/io/ply_io.h>
#include <pcl/io/obj_io.h>
#include <pcl/io/vtk_lib_io.h>

#include <opencv2/video/tracking.hpp>

#include "kalman_smoothing.h"

int KalmanSmoother::readMeshes(const std::string &_fileName)
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
            // std::cerr<<line<<"\n";
            pcl::PolygonMesh mesh;

            if (pcl::io::loadPolygonFile(line, mesh) == -1)
            {
                std::cerr << "Mesh loader: couldn't open file " << line << std::endl;
                return -2;
            }
            meshes_.push_back(mesh);

        } while (!inputfile.eof());
    }
    else
    {
        std::cerr << "Mesh sequence reader: couldn't open file " << _fileName << std::endl;
        return -1;
    }

    std::cerr << "Number of meshes loaded: " << meshes_.size() << std::endl;

    this->setupClouds();
    std::cerr << "Number of clouds loaded: " << clouds_.size() << std::endl;

    return 0;
}

int KalmanSmoother::ReadRegions(const std::string &file_name)
{
  if (meshNames_.empty())
  {
    std::cerr << "Load mesh sequence before finding keyframes!"<<std::endl;
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
      keyFrameIndices_.push_back(std::stoi(value));
    }

    // Print regions
    std::cerr << "Regions: \n";
    for (size_t region = 0; region < regions_.size(); region++)
    {
      std::cerr 
        <<regions_[region].first<< " "<<regions_[region].second<<", kf:"
        <<keyFrameIndices_[region]<<"\n";
    }
    

    // Safety checks
    if(regions_.begin()->first != 0)
    {
      std::cerr << "Warning: Check regions.txt.\n"
        << "\tFirst index must always be 0\n";
      regions_.begin()->first = 0;
    }

    if(regions_.back().second + 1 < meshNames_.size())
    {
      std::cerr
        << "Warning: End of specified regions is less than input meshes!"
        <<"["<<regions_.back().second + 1<<","<<meshNames_.size()<<"]\n"
        << "\tEnd of last region will be extended to end of sequence.\n";
        regions_.end()->second = meshNames_.size() - 1;
    }
  }
  else
  {

    std::cerr <<  "Keyframe reader: couldn't open file " << file_name<<"\n";
    return -1;
  }


  return 1;
}

void KalmanSmoother::setupClouds()
{

    for (unsigned int i = 0; i < meshes_.size(); i++)
    {

        pcl::PolygonMesh mesh = meshes_[i];
        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
        pcl::fromPCLPointCloud2(mesh.cloud, *cloud);

        clouds_.push_back(cloud);
    }
}

void KalmanSmoother::smoothLinearly()
{

    // Allocate space for the smooth clouds
    for (unsigned int i = 0; i < meshes_.size(); i++)
    {
        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
        smoothClouds_.push_back(cloud);
    }

    // First find the mesh interval
    std::vector<unsigned int>::iterator kfit = keyFrameIndices_.begin();
    for (; kfit != keyFrameIndices_.end(); ++kfit)
    {

        // First and last indices in the tracked interval
        unsigned int first = *kfit;
        unsigned int last;

        if (kfit + 1 != keyFrameIndices_.end())
        {
            last = *(kfit + 1) - 1; // The one previous to the next keyframe
        }
        else
        {
            last = meshes_.size() - 1; // The last one
        }

        std::cerr << first << "/" << last << std::endl;

        // We iterate through the meshes in this specific range
        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr keycloud = clouds_[first];
        unsigned int nvtx = keycloud->height * keycloud->width;

        // Store the keycloud in the smoothclouds vector, as this one does not change
        smoothClouds_[first] = keycloud;

        // This vector will store all equivalent vertices in the mesh sequence
        std::vector<pcl::PointXYZRGBNormal> vertexseq(last - first + 1);
        std::vector<pcl::PointXYZRGBNormal> smoothvertexseq(last - first + 1);

        for (unsigned int i = 0; i < nvtx; i++)
        {

            unsigned int k = 0;
            for (unsigned int j = first; j <= last; j++, k++)
            {

                const pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud = clouds_[j];
                const pcl::PointXYZRGBNormal &point = cloud->points[i];
                vertexseq[k] = point;

            } // end per cloud

            // vertexseq is already filled here
            // We use a 3D Kalman filter to smooth the positions
            this->filterKalman(vertexseq, smoothvertexseq);

            // Assign the new vertex to the whole sequence
            for (unsigned int j = first; j <= last; j++)
            {
                pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud = smoothClouds_[j];
                cloud->push_back(smoothvertexseq[j - first]);

            } // end per cloud

        } // end per vertex
    }
}

void KalmanSmoother::smoothInRegions()
{

    // Allocate space for the smooth clouds
    for (unsigned int i = 0; i < meshes_.size(); i++)
    {
        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
        smoothClouds_.push_back(cloud);
    }

    // First find the mesh interval
    for (unsigned int r = 0; r < regions_.size(); r++)
    {
        // We do two passes:
        // KF -> Rinit
        // KF -> Rend
        unsigned int first = regions_[r].first;
        unsigned int last = regions_[r].second;

        std::cerr << first << "/" << last << std::endl;

        // We iterate through the meshes in this specific range
        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr keycloud = clouds_[first];
        unsigned int nvtx = keycloud->height * keycloud->width;

        // Store the keycloud in the smoothclouds vector, as this one does not change
        smoothClouds_[first] = keycloud;

        // This vector will store all equivalent vertices in the mesh sequence
        std::vector<pcl::PointXYZRGBNormal> vertexseq(last - first + 1);
        std::vector<pcl::PointXYZRGBNormal> smoothvertexseq(last - first + 1);
        
        for (unsigned int i = 0; i < nvtx; i++)
        {

            unsigned int k = 0;
            for (unsigned int j = first; j <= last; j++, k++)
            {

                const pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud = clouds_[j];
                const pcl::PointXYZRGBNormal &point = cloud->points[i];
                vertexseq[k] = point;

            } // end per cloud

            // vertexseq is already filled here
            // We use a 3D Kalman filter to smooth the positions
            this->filterKalman(vertexseq, smoothvertexseq);

            // Assign the new vertex to the whole sequence
            for (unsigned int j = first + 1; j <= last; j++)
            {
                pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud = smoothClouds_[j];
                cloud->push_back(smoothvertexseq[j - first]);

            } // end per cloud
        } // end per vertex
    }
}

void KalmanSmoother::filterKalman(const std::vector<pcl::PointXYZRGBNormal> &_input, std::vector<pcl::PointXYZRGBNormal> &_output) const
{

    // Initial values

    double dt = 0.04; // time between measurements (1/FPS = 1/25)

    cv::KalmanFilter KF(9, 3, 0, CV_64F);

    double proc_nc_max = 1.0;
    double proc_nc_min = 1e-5; 

    double max_v = DBL_MIN;
    double min_v = DBL_MAX;

    std::vector<double> v_normd;

    // analyze displacement magnitude
    for (unsigned int i = 1; i < _input.size(); i++)
    {
        double v = 0;
        const pcl::PointXYZRGBNormal &pclinitial = _input[i];
        const pcl::PointXYZRGBNormal &pclprev = _input[i-1];

        v += std::abs(double(pclinitial.x - pclprev.x));
        v += std::abs(double(pclinitial.y - pclprev.y));
        v += std::abs(double(pclinitial.z - pclprev.z));

        max_v = v > max_v ? v : max_v;
        min_v = v < min_v ? v : min_v;

        v_normd.push_back(v);
    }

    // normalize displacements between [0, nc_max - nc_min]
    for (size_t i = 0; i < v_normd.size(); i++)
    {
        v_normd[i] = ((v_normd[i] - min_v)*(proc_nc_max-proc_nc_min))/(max_v-min_v);
    }
    // zero pad the front
    v_normd.insert(v_normd.begin(), 0);
    
    cv::setIdentity(KF.measurementNoiseCov, cv::Scalar::all(1e-3)); // set measurement noise

    cv::setIdentity(KF.errorCovPost, cv::Scalar::all(1)); // error covariance

    KF.transitionMatrix.at<double>(0, 3) = dt;
    KF.transitionMatrix.at<double>(1, 4) = dt;
    KF.transitionMatrix.at<double>(2, 5) = dt;
    KF.transitionMatrix.at<double>(3, 6) = dt;
    KF.transitionMatrix.at<double>(4, 7) = dt;
    KF.transitionMatrix.at<double>(5, 8) = dt;
    KF.transitionMatrix.at<double>(0, 6) = 0.5 * dt * dt;
    KF.transitionMatrix.at<double>(1, 7) = 0.5 * dt * dt;
    KF.transitionMatrix.at<double>(2, 8) = 0.5 * dt * dt;

    KF.measurementMatrix.at<double>(0, 0) = 1; // x
    KF.measurementMatrix.at<double>(1, 1) = 1; // y
    KF.measurementMatrix.at<double>(2, 2) = 1; // z

    for (unsigned int i = 0; i < _input.size(); i++)
    {
        
        
        cv::setIdentity(KF.processNoiseCov, cv::Scalar::all(proc_nc_min + v_normd[i]));     // set process noise

        cv::Mat measurement = cv::Mat::zeros(3, 1, CV_64F);

        const pcl::PointXYZRGBNormal &pclinitial = _input[i];
        cv::Point3f initial(pclinitial.x, pclinitial.y, pclinitial.z);

        measurement.at<double>(0) = initial.x; // Initial X observation
        measurement.at<double>(1) = initial.y; // Initial Y observation
        measurement.at<double>(2) = initial.z; // Initial Z observation

        // First predict, to update the internal statePre variable
        cv::Mat prediction = KF.predict();
        // The "correct" phase that is going to use the predicted value and our measurement
        cv::Mat estimated = KF.correct(measurement);
        // Estimated translation
        cv::Mat finalE(3, 1, CV_64F);
        finalE.at<double>(0) = estimated.at<double>(0);
        finalE.at<double>(1) = estimated.at<double>(1);
        finalE.at<double>(2) = estimated.at<double>(2);

        pcl::PointXYZRGBNormal finalP;
        finalP.x = finalE.at<double>(0);
        finalP.y = finalE.at<double>(1);
        finalP.z = finalE.at<double>(2);

        _output[i] = finalP;
    }
}

std::vector<pcl::PolygonMesh> KalmanSmoother::setupSmoothMeshes(bool _export) const
{

    std::vector<pcl::PolygonMesh> smoothmeshes;

    for (unsigned int i = 0; i < meshes_.size(); i++)
    {

        pcl::PCLPointCloud2 nucloud;
        pcl::toPCLPointCloud2(*smoothClouds_[i], nucloud);

        pcl::PolygonMesh mesh;

        mesh.cloud = nucloud;
        mesh.polygons = meshes_[i].polygons;

        if (_export)
        {
            std::string outname = meshNames_[i].substr(0, meshNames_[i].size() - 4);
            pcl::io::saveOBJFile(outname + "_s.obj", mesh);
        }

        smoothmeshes.push_back(mesh);
    }

    return smoothmeshes;
}
