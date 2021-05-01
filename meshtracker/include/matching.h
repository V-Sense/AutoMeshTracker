#ifndef MATCHING_H
#define MATCHING_H

#include "cloud_processing.h"
#include "mesh_processing.h"
#include "log.h"

// <point index, confidence[0,1]>
typedef std::pair<unsigned int, float>                               match_conf;
typedef pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr                   CloudPtr;

// Use to find correspondence pairs between point clouds or meshes.
// Functions can return confidences or use KNN to assume matches where 
// none are found. 
class Matcher
{
private:
    // TODO: Getter setter for source target with template func
    // to accept either mesh or pt cloud

public:
    Matcher(){};
    ~Matcher(){};

    // Simple pcl nearest neighbour search that will match 1-to-1 for each input
    std::vector<int> GetKNNMatches(
        const CloudPtr &        source,
        const CloudPtr &         query,
        const int      & search_length = 10);

    // Simple pcl nearest neighbour search that will match 1-to-1 for each input   
    std::vector<int> GetKNNMatches(
        const pcl::PolygonMesh &        source,
        const pcl::PolygonMesh &         query,
        const int              & search_length = 10);

    // Simple pcl nearest neighbour search that will match 1-to-1 for each input
    std::vector<int> GetKNNMatches(
        const pcl::PointXYZRGBNormal &    source, 
        const CloudPtr               &     query,
        const int                    & tolerance,
        int                        search_length = 10);

    // Returns only the point correspondences which conform to strict matching rules.
    // i.e. sparse matching only
    std::vector<int> GetConstrainedMatches(
        const pcl::PolygonMesh & source,
        const pcl::PolygonMesh &  query,
        double          align_tolerance = 0.25);

    // Returns only the point correspondences which conform to strict matching rules.
    // i.e. sparse matching only
    std::vector<int> GetConstrainedMatches(
        const CloudPtr & source,
        const CloudPtr &  query,
        double  align_tolerance = 0.25,
        float            radius =    0); // flag to auto calculate max_radius

    // Returns only the point correspondences which conform to strict matching rules.
    // i.e. sparse matching only
    std::vector<int> GetConstrainedMatches(
        const pcl::PointXYZRGBNormal &    source, 
        const CloudPtr               &     query,
        const double                 & tolerance,
        int                        search_length = 10);

    // Performs KNN search from source to query, returns exact 1-to-1 matches as well
    // as confidence between [0, 1.0] based on normal alignment
    std::vector<match_conf> GetMatchesAndConfidence(
        CloudPtr &         source, 
        CloudPtr &          query, 
        const float    normal_tol = 0.5,
        const float  max_ray_dist = 2.0); 

    // // Debug method outputs a colorized cloud which illustrates matching
    // void VisualizeMatches(
    //     const std::vector<match_conf>& matches, 
    //     const CloudPtr& source,
    //     const CloudPtr& query,
    //     CloudPtr& _colored_source,
    //     CloudPtr& _colored_query);

    // override for polygonmesh input
    void VisualizeMatches(
        const std::vector<match_conf> &         matches, 
        const pcl::PolygonMesh        &          source,
        const pcl::PolygonMesh        &           query,
        pcl::PolygonMesh              & _colored_source,
        pcl::PolygonMesh              & _colored_query);
    
    // Override for basic KNNMatches
    void VisualizeMatches(
        const std::vector<match_conf> &         matches, 
        const CloudPtr                &          source,
        const CloudPtr                &           query,
        CloudPtr                      & _colored_source,
        CloudPtr                      & _colored_query);
        
    // Override for basic KNNMatches
    void VisualizeMatches(
        const std::vector<int> &         matches, 
        const CloudPtr         &          source,
        const CloudPtr         &           query,
        CloudPtr               & _colored_source,
        CloudPtr               &  _colored_query);
};



#endif