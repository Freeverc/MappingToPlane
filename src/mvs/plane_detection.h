#include "util/ply.h"
#include <cstdlib>
#include <ctime>
#include <float.h>
#include <iostream>

#include <pcl/common/transforms.h>
#include <pcl/features/normal_3d.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h>
#include <pcl/search/search.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/segmentation/organized_multi_plane_segmentation.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <vector>

namespace colmap {
namespace mvs {

struct PlaneDetectionOptions {
  int num_neighbor = 30;

  int min_num = 50;

  int max_num = 100000;

  double smooth_thresh = 10.0;

  double curve_thresh = 3.0;

  // The number of threads used for the Poisson reconstruction.
  int num_threads = -1;

  bool Check() const;
};

void down_sample(pcl::PointCloud<pcl::PointXYZ>::Ptr& point_cloud,
                 pcl::PointCloud<pcl::PointXYZ>::Ptr& sampled_point_cloud);

void filter(pcl::PointCloud<pcl::PointXYZ>::Ptr& point_cloud,
            pcl::PointCloud<pcl::PointXYZ>::Ptr& filtered_point_cloud);

void plane_detection(pcl::PointCloud<pcl::PointXYZ>::Ptr& point_cloud,
                     pcl::PointCloud<pcl::PointXYZ>::Ptr& inner_point_cloud);

void multi_plane_detection(
    pcl::PointCloud<pcl::PointXYZ>::Ptr& point_cloud,
    pcl::PointCloud<pcl::PointXYZL>::Ptr& inner_point_cloud,
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr& colored_point_cloud,
    std::vector<std::vector<float>>& plane_list);

void IndiceToClustered(
    pcl::PointCloud<pcl::PointXYZ>::Ptr& point_cloud,
    std::vector<pcl::PointIndices>& cluster_indices,
    pcl::PointCloud<pcl::PointXYZL>::Ptr& clustered_point_cloud);

// Perform plane detection.
bool PlaneDetection(const PlaneDetectionOptions& options,
                    const std::string& input_path,
                    const std::string& output_path,
                    std::vector<PlyPoint>& plane_points,
                    std::vector<std::vector<float>>& plane_list);
}  // namespace mvs
}  // namespace colmap