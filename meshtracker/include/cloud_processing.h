#ifndef CLOUD_PROCESSING_H
#define CLOUD_PROCESSING_H

#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/search/kdtree.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/filters/bilateral.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/random_sample.h>
#include <pcl/registration/icp.h>
#include <pcl/surface/mls.h>

#include <Eigen/Geometry>

#include "log.h"

typedef pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr CloudPtr;


namespace cloud_processing
{
    struct BBox
    {

        double min_x     = DBL_MAX;
        double min_y     = DBL_MAX;
        double min_z     = DBL_MAX;
        double max_x     = -DBL_MAX;
        double max_y     = -DBL_MAX;
        double max_z     = -DBL_MAX;
        double cx        = 0;
        double cy        = 0;
        double cz        = 0;
        double l_diag_sq = 0;

        BBox(const pcl::PolygonMesh& mesh){ this->get(mesh);};
        BBox(const CloudPtr& cloud){ this->get(cloud);};

        void get(const pcl::PolygonMesh& mesh)
        {
            CloudPtr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
            pcl::fromPCLPointCloud2(mesh.cloud,*cloud);

            for (size_t pt = 0; pt < cloud->points.size(); pt++)
            {
                pcl::PointXYZRGBNormal p = cloud->points[pt];

                if (p.x < this->min_x) this->min_x = p.x; 
                if (p.x > this->max_x) this->max_x = p.x; 
                if (p.y < this->min_y) this->min_y = p.y; 
                if (p.y > this->max_y) this->max_y = p.y; 
                if (p.z < this->min_z) this->min_z = p.z; 
                if (p.z > this->max_z) this->max_z = p.z; 
            }

            this->cx = this->min_x + (this->max_x - this->min_x)*0.5;
            this->cy = this->min_y + (this->max_y - this->min_y)*0.5;
            this->cz = this->min_z + (this->max_z - this->min_z)*0.5;

            this->l_diag_sq = (max_x - min_x) + (max_y - min_y);
        };

        void get(const CloudPtr& cloud)
        {
            for (size_t pt = 0; pt < cloud->points.size(); pt++)
            {
                pcl::PointXYZRGBNormal p = cloud->points[pt];

                if (p.x < this->min_x) this->min_x = p.x; 
                if (p.x > this->max_x) this->max_x = p.x; 
                if (p.y < this->min_y) this->min_y = p.y; 
                if (p.y > this->max_y) this->max_y = p.y; 
                if (p.z < this->min_z) this->min_z = p.z; 
                if (p.z > this->max_z) this->max_z = p.z; 
            }

            this->cx = this->min_x + (this->max_x - this->min_x)*0.5;
            this->cy = this->min_y + (this->max_y - this->min_y)*0.5;
            this->cz = this->min_z + (this->max_z - this->min_z)*0.5;

            this->l_diag_sq = (max_x - min_x) + (max_y - min_y);
        };

        void print()
        {
            std::cout<<"Bbox:"<<std::endl;
            std::cout<<"x: ["<<this->min_x<<", "<<this->max_x<<"] l:"<<this->max_x - this->min_x<<std::endl;
            std::cout<<"y: ["<<this->min_y<<", "<<this->max_y<<"] l:"<<this->max_y - this->min_y<<std::endl;
            std::cout<<"z: ["<<this->min_z<<", "<<this->max_z<<"] l:"<<this->max_z - this->min_z<<std::endl;
            std::cout<<"c: ["<<this->cx<<", "<<this->cy<<", "<<this->cz<<"]"<<std::endl;
        };
    };
}

class CloudProcessing
{
private:

    const int kSamplesUpperLimit = 600;
    const int kSamplesLowerLimit =  50;
    
public:
    CloudProcessing(){};
    ~CloudProcessing(){};

    // Use Moving Least Squares to smooth a point cloud
    CloudPtr SmoothCloudMLS(
        const CloudPtr & cloud);

    // returns avg distance between each point and it's K-nearest neighbour
    double ComputeCloudResolution(
        const CloudPtr & cloud);

    // Checks the first point to see if xyz normals are non-zero
    bool NormalsExist(
        const CloudPtr & cloud);

    // Calulcate and apply normals to input cloud
    void CalculatePointCloudNormals(
        CloudPtr & cloud);

    // Filters outliers using a neighbours-within-radius condition. Function returns
    // the indices of the outlier points
    std::vector<int> FindOutliers(
        const CloudPtr & cloud);

    // Calculates hdorf from a to b and returns avg dist per point
    double ComputeHausdorffDistance(
        const CloudPtr & cloud_a,
        const CloudPtr & cloud_b);

    // Calculates distance from point to cloud
    double ComputeHausdorffDistance(
        const pcl::PointXYZRGBNormal &    pt,
        const CloudPtr               & cloud);

    // Implements pcl VoxelGrid resampling
    CloudPtr ResamplePointCloud(
        CloudPtr & cloud);

    // Implements pcl VoxelGrid uniform sampling
    CloudPtr UniformsamplePointCloud(
        const CloudPtr &         cloud,
        std::optional<int> target_size);

    // Implements Random subsampling
    CloudPtr RandomsamplePointCloud(
        CloudPtr           &   cloud,
        const unsigned int & samples);

    // Translate moving cloud such that center of bbox aligns with 
    // centre of target bbox
    CloudPtr AlignCentres(
        const CloudPtr & moving,
        const CloudPtr &  fixed);

    // transform cloud to origin
    CloudPtr TransformToOrigin(
        const CloudPtr & cloud);

    // Perform rigid ICP and return aligned cloud
    CloudPtr RigidICP(
        const CloudPtr &          moving,
        const CloudPtr &           fixed,
        const double   &   max_pair_dist =  1.0,
        const int      &       max_iters =   50,
        const double   &   tform_epsilon = 1e-8,
        const double   & fitness_epsilon = 0.05);

    inline int GetSamplesUpperLimit() const
    {
        return kSamplesUpperLimit;
    }

    inline int GetSamplesLowerLimit() const
    {
        return kSamplesLowerLimit;
    }

    void ColorizeCloud(
        CloudPtr          & cloud, 
        std::array<int,3> &   rgb);
};


#endif