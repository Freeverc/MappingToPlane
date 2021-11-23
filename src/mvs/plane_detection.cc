#include "mvs/plane_detection.h"

namespace colmap {
namespace mvs {

void down_sample(pcl::PointCloud<pcl::PointXYZ>::Ptr& point_cloud,
                 pcl::PointCloud<pcl::PointXYZ>::Ptr& sampled_point_cloud) {
  std::cout << "clustering by euclid..." << std::endl;

  pcl::VoxelGrid<pcl::PointXYZ> vg;
  vg.setInputCloud(point_cloud);
  vg.setLeafSize(0.02f, 0.02f, 0.02f);
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
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr& colored_point_cloud,
    std::vector<Plane>& plane_list) {
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

  IndiceToClustered(point_cloud, clusters, inner_point_cloud);
  colored_point_cloud = reg.getColoredCloud();

  int c_id = 0;
  for (std::vector<pcl::PointIndices>::const_iterator it = clusters.begin();
       it != clusters.end(); ++it) {
    pcl::PointCloud<pcl::PointXYZ>::Ptr point_cluser(
        new pcl::PointCloud<pcl::PointXYZ>);
    std::vector<PlyPoint> ply_points;
    for (std::vector<int>::const_iterator pit = it->indices.begin();
         pit != it->indices.end(); ++pit) {
      pcl::PointXYZ p;
      p.x = (*point_cloud)[*pit].x;
      p.y = (*point_cloud)[*pit].y;
      p.z = (*point_cloud)[*pit].z;
      point_cluser->push_back(p);

      PlyPoint py;
      py.x = (*point_cloud)[*pit].x;
      py.y = (*point_cloud)[*pit].y;
      py.z = (*point_cloud)[*pit].z;
      ply_points.push_back(py);
    }

    pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
    pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);

    // Create the segmentation object
    pcl::SACSegmentation<pcl::PointXYZ> seg;
    seg.setOptimizeCoefficients(true);
    seg.setModelType(pcl::SACMODEL_PLANE);
    seg.setMethodType(pcl::SAC_RANSAC);
    seg.setDistanceThreshold(0.002);

    seg.setInputCloud(point_cluser);
    seg.segment(*inliers, *coefficients);

    float a = coefficients->values[0];
    float b = coefficients->values[1];
    float c = coefficients->values[2];
    float d = coefficients->values[3];

    float s2 = a * a + b * b + c * c;
    float s = std::sqrt(s2);

    Plane pl;
    pl.para = coefficients->values;
    pl.direction = std::atan2(b, a) / M_PI * 180;
    pl.inclination = std::atan2(s, c) / M_PI * 180;
    pl.area = it->indices.size();
    pl.inner_points = ply_points;

    plane_list.emplace_back(pl);

    std::cout << "Cluster id : " << c_id << std::endl;
    std::cout << "Cluster size : " << point_cluser->size() << std::endl;
    std::cout << "Inner size : " << inliers->indices.size() << std::endl;
    std::cout << "Model coefficients: " << a << " " << b << " " << c << " " << d
              << std::endl;
    ++c_id;
  }
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
                    const std::string& output_path,
                    std::vector<PlyPoint>& plane_points,
                    std::vector<Plane>& plane_list) {
  std::cout << "Point cloud measuring : " << std::endl;

  std::string point_cloud_path = input_path + "/fused.ply";
  std::string sampled_path = output_path + "/sampled.ply";
  std::string filtered_path = output_path + "/filtered.ply";
  std::string inner_path = output_path + "/plane_points.ply";
  std::string color_path = output_path + "/cluster_points.ply";

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
                        color_point_cloud, plane_list);

  std::cout << "Inner : " << inner_point_cloud->size() << std::endl;
  std::cout << "Color : " << color_point_cloud->size() << std::endl;
  pcl::io::savePLYFileASCII(inner_path, *inner_point_cloud);
  pcl::io::savePLYFileASCII(color_path, *color_point_cloud);

  for (int i = 0; i < color_point_cloud->points.size(); i++) {
    PlyPoint p;
    p.x = color_point_cloud->points[i].x;
    p.y = color_point_cloud->points[i].y;
    p.z = color_point_cloud->points[i].z;
    p.r = color_point_cloud->points[i].r;
    p.g = color_point_cloud->points[i].g;
    p.b = color_point_cloud->points[i].b;
    plane_points.push_back(p);
  }
}

bool GenerateDEM(const std::string& input_path,
                 const std::string& output_path) {
  std::cout << "Generating DEM : " << std::endl;

  std::string point_cloud_path = input_path + "/filtered.ply";
  std::string dem_path = output_path + "/dem.jpg";

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

  float a = coefficients->values[0];
  float b = coefficients->values[1];
  float c = coefficients->values[2];
  float d = coefficients->values[3];
  std::cout << "ground : " << a << " " << b << " " << c << " " << d
            << std::endl;

  float h = std::sqrt(a * a + b * b);
  float s = std::sqrt(a * a + b * b + c * c);
  float sin_ab = b / h;
  float cos_ab = a / h;
  float sin_z = c / s;
  float cos_z = h / s;

  Eigen::Matrix<float, 3, 3, Eigen::RowMajor> R;
  R(0, 0) = sin_ab;
  R(0, 1) = sin_z * cos_ab;
  R(0, 2) = a / s;
  R(1, 0) = -cos_ab;
  R(1, 1) = sin_z * sin_ab;
  R(1, 2) = b / s;
  R(2, 0) = 0;
  R(2, 1) = -cos_z;
  R(2, 2) = c / s;
  Eigen::Matrix<float, 3, 3, Eigen::RowMajor> R_inv = R.inverse();

  float min_d = FLT_MAX;
  float max_d = FLT_MIN;
  float min_xx = FLT_MAX;
  float max_xx = FLT_MIN;
  float min_yy = FLT_MAX;
  float max_yy = FLT_MIN;
  float min_zz = FLT_MAX;
  float max_zz = FLT_MIN;
  for (int i = 0; i < point_cloud->points.size(); ++i) {
    double x = (*point_cloud)[i].x;
    double y = (*point_cloud)[i].y;
    double z = (*point_cloud)[i].z;
    Eigen::Vector3f xyz(x, y, z);
    Eigen::Vector3f xyz_p = R_inv * xyz;
    float dis = (a * x + b * y + c * z + d) / s;
    min_d = std::min(dis, min_d);
    max_d = std::max(dis, max_d);
    min_xx = std::min(xyz_p(0), min_xx);
    max_xx = std::max(xyz_p(0), max_xx);
    min_yy = std::min(xyz_p(1), min_yy);
    max_yy = std::max(xyz_p(1), max_yy);
    min_zz = std::min(xyz_p(2), min_zz);
    max_zz = std::max(xyz_p(2), max_zz);
  }

  double reso = 0.04;
  Bitmap bitmap;
  int width = static_cast<int>((max_xx - min_xx) / reso + 1);
  int height = static_cast<int>((max_yy - min_yy) / reso + 1);
  int bar_width = width * 0.4;
  width += bar_width;
  int houdu = static_cast<int>((max_zz - min_zz) / reso + 1);
  bitmap.Allocate(width, height, true);
  std::vector<std::vector<float>> d_list(height,
                                         std::vector<float>(width, min_d));

  const BitmapColor<uint8_t> white_color(255, 255, 255);
  for (int x = 0; x < width; ++x) {
    for (int y = 0; y < height; ++y) {
      bitmap.SetPixel(x, y, white_color);
    }
  }

  std::cout << "x: " << min_xx << " " << max_xx << std::endl;
  std::cout << "y: " << min_yy << " " << max_yy << std::endl;
  std::cout << "z: " << min_zz << " " << max_zz << std::endl;
  std::cout << "d: " << min_d << " " << max_d << std::endl;
  std::cout << "w : " << width << " h :  " << height << " z : " << houdu
            << std::endl;

  int block_w = static_cast<int>(width * 0.06);
  int block_h = static_cast<int>(height * 0.06);

  // Draw bar
  cv::Mat image(height, bar_width, CV_8UC3, cv::Scalar(255, 255, 255));
  for (int i = 0; i < 10; i++) {
    float gray = static_cast<float>(i) / 10;
    std::string dis = std::to_string(gray * (max_d - min_d));
    cv::Rect r(10, 10 + block_h * i, block_w, block_h);
    cv::Scalar c(255 * JetColormap::Red(gray), 255 * JetColormap::Green(gray),
                 255 * JetColormap::Blue(gray));
    cv::rectangle(image, r, c, -1);
    cv::putText(image, dis, cv::Point(20 + block_w, 12 + block_h * (i + 1)),
                cv::FONT_HERSHEY_SIMPLEX, 0.02 / reso, cv::Scalar(0, 0, 0));
  }

  for (int x = 0; x < bar_width; x++) {
    for (int y = 0; y < height; y++) {
      cv::Vec3b c = image.at<cv::Vec3b>(y, x);
      const BitmapColor<float> color(c[0], c[1], c[2]);
      bitmap.SetPixel(width - bar_width + x, y, color.Cast<uint8_t>());
    }
  }

  for (int i = 0; i < point_cloud->points.size(); ++i) {
    double x = (*point_cloud)[i].x;
    double y = (*point_cloud)[i].y;
    double z = (*point_cloud)[i].z;
    float dis = a * x + b * y + c * z + d / s;

    Eigen::Vector3f xyz(x, y, z);
    Eigen::Vector3f xyz_p = R_inv * xyz;

    int x_id = static_cast<int>((xyz_p(0) - min_xx) / reso);
    int y_id = static_cast<int>((xyz_p(1) - min_yy) / reso);

    if (i < 10) {
      std::cout << "p1: " << xyz(0) << " " << xyz(1) << " " << xyz(2) << " "
                << dis << " - " << dis - xyz_p(2) << std::endl;
      std::cout << "p2: " << xyz_p(0) << " " << xyz_p(1) << " " << xyz_p(2)
                << std::endl;
      std::cout << "id : " << x_id << " " << y_id << std::endl;
    }

    if (dis > d_list[y_id][x_id]) {
      const float gray = 1 - (dis - min_d) / (max_d - min_d);
      const BitmapColor<float> color(255 * JetColormap::Red(gray),
                                     255 * JetColormap::Green(gray),
                                     255 * JetColormap::Blue(gray));
      bitmap.SetPixel(x_id, y_id, color.Cast<uint8_t>());
      d_list[y_id][x_id] = dis;
    }
  }

  bitmap.Write(dem_path);
}

}  // namespace mvs
}  // namespace colmap