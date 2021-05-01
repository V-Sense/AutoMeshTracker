#ifndef KALMAN_SMOOTHING_H
#define KALMAN_SMOOTHING_H

class KalmanSmoother
{

  public:
    KalmanSmoother(){};
    ~KalmanSmoother(){};

    inline void setMeshNames(const std::vector<std::string> & _meshNames)
    {
        meshNames_ = _meshNames;
    }

    inline const std::vector<std::string> & getMeshNames() const
    {
        return meshNames_;
    }

    inline void setMeshes(const std::vector<pcl::PolygonMesh> & _meshes)
    {
        meshes_ = _meshes;
    }

    inline const std::vector<pcl::PolygonMesh> & getMeshes() const
    {
        return meshes_;
    }

    inline void setClouds(const std::vector<pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr> & _clouds)
    {
        clouds_ = _clouds;
    }

    inline const std::vector<pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr> & getClouds() const
    {
        return clouds_;
    }

    inline void setKeyFrameIndices(const std::vector<unsigned int> & _keyFrameIndices)
    {
        keyFrameIndices_ = _keyFrameIndices;
    }
    
    inline const std::vector<unsigned int> &getKeyFrameIndices() const
    {
        return keyFrameIndices_;
    }

    // Read the mesh sequence and keyframes from a text file
    int readMeshes(const std::string & _fileName);

    // Finds the keyframe indices in the mesh sequence
    int ReadRegions(const std::string & _fileName);

    // Smooths the mesh sequence keyframe to keyframe
    // It is simpler but not the rigth way
    void smoothLinearly();

    //
    void smoothInRegions();

    // 3D Kalman filter to smooth the sequence
    void filterKalman(
        const std::vector<pcl::PointXYZRGBNormal> &  _input, 
        std::vector<pcl::PointXYZRGBNormal>       & _output) const;

    // Takes the smooth clouds and builds meshes
    std::vector<pcl::PolygonMesh> setupSmoothMeshes(bool _export = true) const;

  private:
    //
    void setupClouds();

    std::vector<std::string>                                     meshNames_;
    std::vector<pcl::PolygonMesh>                                   meshes_;
    std::vector<pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr>       clouds_;
    std::vector<unsigned int>                              keyFrameIndices_;
    std::vector<std::pair<unsigned int, unsigned int>>             regions_;
    std::vector<pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr> smoothClouds_;
};

#endif
