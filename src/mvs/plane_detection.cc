#include "mvs/plane_detection.h"

namespace colmap {
namespace mvs {

void down_sample(pcl::PointCloud<pcl::PointXYZ>::Ptr& point_cloud,
                 pcl::PointCloud<pcl::PointXYZ>::Ptr& sampled_point_cloud) {
  std::cout << "clustering by euclid..." << std::endl;

  pcl::VoxelGrid<pcl::PointXYZ> vg;
  vg.setInputCloud(point_cloud);
  vg.setLeafSize(0.05f, 0.05f, 0.05f);
  vg.filter(*sampled_point_cloud);
  std::cout << "PointCloud after sampling has: " << sampled_point_cloud->size()
            << " data points." << std::endl;  //*
}

void filter(pcl::PointCloud<pcl::PointXYZ>::Ptr& point_cloud,
            pcl::PointCloud<pcl::PointXYZ>::Ptr& filtered_point_cloud) {
  // Filtering.
  pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
  sor.setInputCloud(point_cloud);
  sor.setMeanK(50);
  sor.setStddevMulThresh(1);
  sor.filter(*filtered_point_cloud);
  std::cout << "PointCloud after filtering has: "
            << filtered_point_cloud->size() << " data points." << std::endl;
}

void plane_detection(pcl::PointCloud<pcl::PointXYZ>::Ptr& point_cloud,
                     pcl::PointCloud<pcl::PointXYZ>::Ptr& inner_point_cloud) {
  pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
  pcl::PointIndices::Ptr inliers(new pcl::PointIndices);

  // Create the segmentation object
  pcl::SACSegmentation<pcl::PointXYZ> seg;

  // Optional
  seg.setOptimizeCoefficients(true);
  // Mandatory
  seg.setModelType(pcl::SACMODEL_PLANE);
  seg.setMethodType(pcl::SAC_RANSAC);
  seg.setDistanceThreshold(0.2);

  seg.setInputCloud(point_cloud);
  seg.segment(*inliers, *coefficients);

  for (std::vector<int>::const_iterator pit = inliers->indices.begin();
       pit != inliers->indices.end(); ++pit) {
    pcl::PointXYZ p;
    p.x = (*point_cloud)[*pit].x;
    p.y = (*point_cloud)[*pit].y;
    p.z = (*point_cloud)[*pit].z;
    inner_point_cloud->push_back(p);
  }

  float a = coefficients->values[0];
  float b = coefficients->values[1];
  float c = coefficients->values[2];
  float d = coefficients->values[3];

  std::cout << "Inner num : " << inner_point_cloud->size() << std::endl;
  std::cout << "Model coefficients: " << a << " " << b << " " << c << " " << d
            << std::endl;
  pcl::io::savePLYFileASCII("../results/plane.ply", *inner_point_cloud);

  Eigen::Matrix4f rotate_x = Eigen::Matrix4f::Identity();
  float theta_x = std::acos(c / std::sqrt(c * c + b * b));
  rotate_x(1, 1) = std::cos(theta_x);
  rotate_x(1, 2) = std::sin(theta_x);
  rotate_x(2, 1) = -std::sin(theta_x);
  rotate_x(2, 2) = std::cos(theta_x);
  std::cout << "rotate x : " << std::endl;
  std::cout << rotate_x << std::endl;
  pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_x_cloud(
      new pcl::PointCloud<pcl::PointXYZ>());
  pcl::transformPointCloud(*inner_point_cloud, *transformed_x_cloud, rotate_x);
  pcl::io::savePLYFileASCII("../results/transformed_x.ply",
                            *transformed_x_cloud);

  Eigen::Matrix4f rotate_y = Eigen::Matrix4f::Identity();
  float theta_y = std::acos(c / std::sqrt(c * c + a * a));
  rotate_y(0, 0) = std::cos(theta_y);
  rotate_y(0, 2) = -std::sin(theta_y);
  rotate_y(2, 0) = std::sin(theta_y);
  rotate_y(2, 2) = std::cos(theta_y);
  std::cout << "rotate y : " << std::endl;
  std::cout << rotate_y << std::endl;
  pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud(
      new pcl::PointCloud<pcl::PointXYZ>());
  pcl::transformPointCloud(*transformed_x_cloud, *transformed_cloud, rotate_y);
  pcl::io::savePLYFileASCII("../results/transformed.ply", *transformed_cloud);

  float min_x, min_y, min_z, max_x, max_y, max_z;
  min_x = min_y = min_z = FLT_MAX;
  max_x = max_y = max_z = FLT_MIN;
  for (int i = 0; i < transformed_cloud->size(); ++i) {
    float x = (*inner_point_cloud)[i].x;
    float y = (*inner_point_cloud)[i].y;
    float z = (*inner_point_cloud)[i].z;
    if (x < min_x) min_x = x;
    if (y < min_y) min_y = y;
    if (z < min_z) min_z = z;
    if (x > max_x) max_x = x;
    if (y > max_y) max_y = y;
    if (z > max_z) max_z = z;
  }
  std::cout << "inner cloud : " << std::endl;
  std::cout << "Min : " << min_x << " " << min_y << " " << min_z << std::endl;
  std::cout << "Max : " << max_x << " " << max_y << " " << max_z << std::endl;
  min_x = min_y = min_z = FLT_MAX;
  max_x = max_y = max_z = FLT_MIN;
  for (int i = 0; i < transformed_cloud->size(); ++i) {
    float x = (*transformed_cloud)[i].x;
    float y = (*transformed_cloud)[i].y;
    float z = (*transformed_cloud)[i].z;
    if (x < min_x) min_x = x;
    if (y < min_y) min_y = y;
    if (z < min_z) min_z = z;
    if (x > max_x) max_x = x;
    if (y > max_y) max_y = y;
    if (z > max_z) max_z = z;
  }

  std::cout << "transformed cloud : " << std::endl;
  std::cout << "Min : " << min_x << " " << min_y << " " << min_z << std::endl;
  std::cout << "Max : " << max_x << " " << max_y << " " << max_z << std::endl;
}

void multi_plane_detection(
    pcl::PointCloud<pcl::PointXYZ>::Ptr& point_cloud,
    pcl::PointCloud<pcl::PointXYZL>::Ptr& inner_point_cloud,
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr& colored_point_cloud) {
  pcl::PassThrough<pcl::PointXYZ> pass;
  pcl::PointCloud<pcl::PointXYZ>::Ptr org_point_cloud(
      new pcl::PointCloud<pcl::PointXYZ>);
  pass.setInputCloud(point_cloud);
  pass.setFilterFieldName("z");
  pass.setFilterLimits(-1000, 1000);
  pass.setKeepOrganized(true);
  pass.filter(*org_point_cloud);

  pcl::PointCloud<pcl::Normal>::Ptr point_cloud_n(
      new pcl::PointCloud<pcl::Normal>);
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
  ne.setInputCloud(org_point_cloud);
  ne.setKSearch(20);
  ne.compute(*point_cloud_n);

  // util::displayPC(viewer, source_pc, 1, 1, 1);
  // viewer->addPointCloudNormals<PointXYZ, Normal>(source_pc, source_n, 64,
  // 0.05, "NN");

  std::vector<pcl::PlanarRegion<pcl::PointXYZ>,
              Eigen::aligned_allocator<pcl::PlanarRegion<pcl::PointXYZ>>>
      regions;
  std::vector<pcl::ModelCoefficients> model_coefficients;
  std::vector<pcl::PointIndices> inlier_indices;
  pcl::PointCloud<pcl::Label>::Ptr labels(new pcl::PointCloud<pcl::Label>);
  std::vector<pcl::PointIndices> label_indices;
  std::vector<pcl::PointIndices> boundary_indices;

  // Segment Planes
  pcl::search::Search<pcl::PointXYZ>::Ptr tree =
      boost::shared_ptr<pcl::search::Search<pcl::PointXYZ>>(
          new pcl::search::KdTree<pcl::PointXYZ>);
  pcl::RegionGrowing<pcl::PointXYZ, pcl::Normal> reg;
  reg.setMinClusterSize(100);
  reg.setMaxClusterSize(10000000);
  reg.setSearchMethod(tree);
  reg.setNumberOfNeighbours(50);
  reg.setInputCloud(point_cloud);
  // reg.setIndices (indices);
  reg.setInputNormals(point_cloud_n);
  reg.setSmoothnessThreshold(4.0 / 180.0 * M_PI);
  reg.setCurvatureThreshold(0.8);

  std::vector<pcl::PointIndices> clusters;
  reg.extract(clusters);
  // int64 mps_start = cv::getTickCount();

  // pcl::OrganizedMultiPlaneSegmentation<pcl::PointXYZ, pcl::Normal,
  // pcl::Label> mps; mps.setAngularThreshold(pcl::deg2rad(10.0));
  // mps.setMinInliers(1000);
  // mps.setInputNormals(point_cloud_n);
  // mps.setInputCloud(org_point_cloud);

  // mps.segmentAndRefine(regions, model_coefficients, inlier_indices, labels,
  //                      label_indices, boundary_indices);

  // double mps_end = cv::getTickCount();
  // PCL_WARN("MPS + Refine took: %f sec \n=== \n",
  //          double(mps_end - mps_start) / cv::getTickFrequency());

  IndiceToClustered(point_cloud, clusters, inner_point_cloud);
  colored_point_cloud = reg.getColoredCloud();
}

void IndiceToClustered(
    pcl::PointCloud<pcl::PointXYZ>::Ptr& point_cloud,
    std::vector<pcl::PointIndices>& cluster_indices,
    pcl::PointCloud<pcl::PointXYZL>::Ptr& clustered_point_cloud) {
  int j = 0;
  for (std::vector<pcl::PointIndices>::const_iterator it =
           cluster_indices.begin();
       it != cluster_indices.end(); ++it) {
    for (std::vector<int>::const_iterator pit = it->indices.begin();
         pit != it->indices.end(); ++pit) {
      pcl::PointXYZL p;
      p.x = (*point_cloud)[*pit].x;
      p.y = (*point_cloud)[*pit].y;
      p.z = (*point_cloud)[*pit].z;
      p.label = j;
      clustered_point_cloud->push_back(p);
    }
    j++;
  }
  clustered_point_cloud->width = clustered_point_cloud->size();
  clustered_point_cloud->height = 1;
  clustered_point_cloud->is_dense = true;
}

bool PlaneDetection(const PlaneDetectionOptions& options,
                    const std::string& input_path,
                    const std::string& output_path) {
  std::cout << "Point cloud measuring : " << std::endl;
  std::cout << "./point_cluster measuring or merge : "
               "(measure, merge)"
            << std::endl;

  std::string point_cloud_path = input_path + "/fused.ply";
  std::string sampled_path = output_path + "/sampled.ply";
  std::string filtered_path = output_path + "/filtered.ply";
  std::string inner_path = output_path + "/inner.ply";
  std::string color_path = output_path + "/colored.ply";

  // Read point cloud.
  pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud(
      new pcl::PointCloud<pcl::PointXYZ>);
  if (pcl::io::loadPLYFile<pcl::PointXYZ>(point_cloud_path, *point_cloud) ==
      -1) {
    PCL_ERROR("Couldn't read file\n");
    return (-1);
  }
  std::cout << "Size : " << point_cloud->width << " " << point_cloud->height
            << std::endl;

  // Do sampling.
  pcl::PointCloud<pcl::PointXYZ>::Ptr sampled_point_cloud(
      new pcl::PointCloud<pcl::PointXYZ>);
  down_sample(point_cloud, sampled_point_cloud);
  std::cout << "Down sampled : " << sampled_point_cloud->size() << std::endl;
  pcl::io::savePLYFileASCII(sampled_path, *sampled_point_cloud);

  // Do filtering.
  pcl::PointCloud<pcl::PointXYZ>::Ptr filtered_point_cloud(
      new pcl::PointCloud<pcl::PointXYZ>);
  filter(sampled_point_cloud, filtered_point_cloud);
  std::cout << "Filtered : " << filtered_point_cloud->size() << std::endl;
  pcl::io::savePLYFileASCII(filtered_path, *filtered_point_cloud);

  // Do plane detecting.
  pcl::PointCloud<pcl::PointXYZL>::Ptr inner_point_cloud(
      new pcl::PointCloud<pcl::PointXYZL>);
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr color_point_cloud(
      new pcl::PointCloud<pcl::PointXYZRGB>);
  multi_plane_detection(filtered_point_cloud, inner_point_cloud,
                        color_point_cloud);
  std::cout << "Inner : " << inner_point_cloud->size() << std::endl;
  std::cout << "Color : " << color_point_cloud->size() << std::endl;
  pcl::io::savePLYFileASCII(inner_path, *inner_point_cloud);
  pcl::io::savePLYFileASCII(color_path, *color_point_cloud);

  // Visualizing.
  // pcl::visualization::PCLVisualizer viewer("Cloud Viewer");
  // // viewer.setBackgroundColor(0, 0, 0);
  // viewer.addPointCloud(color_point_cloud, "inner");
  // while (!viewer.wasStopped()) {
  // }

  // pcl::visualization::CloudViewer viewer("Cluster viewer");
  // viewer.showCloud(color_point_cloud);
  // while (!viewer.wasStopped()) {
  // }
}

}  // namespace mvs
}  // namespace colmap