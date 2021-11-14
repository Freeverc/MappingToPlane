
#ifndef COLMAP_SRC_BASE_GPS_H_
#define COLMAP_SRC_BASE_GPS_H_

#include <vector>

#include <Eigen/Core>
#include <float.h>

#include "util/alignment.h"
#include "util/types.h"

#define AngelToRadian(angle) (angle / 180.0 * M_PI)
#define RadianToAngel(radian) (radian * 180.0 / M_PI)

// GRS-80	6378137.0 m	 6356752.314140 298.257222100882711...
// WGS-84	6378137.0 m	 6356752.314245 298.257223563
// const double EARTH_A = 6371000.0;
// const double EARTH_F_r = 0;

// WGS84标准
const double EARTH_A = 6378137.0;
const double EARTH_F = 298.257223563;
const double EARTH_F_r = 1.0 / EARTH_F;

const double EARTH_RADIUS = 6371000.0;
// const double EARTH_B = EARTH_A * (1.0 - 1.0 / EARTH_F);  // ≈ 6356752.314245

namespace colmap {

struct GPSPoint {
  double latitude;
  double longitude;
  double altitude;

  GPSPoint() : latitude(0), longitude(0), altitude(0){};

  GPSPoint(double lat, double lon, double alt) {
    latitude = lat;
    longitude = lon;
    altitude = alt;
  };
};

struct Point {
  double x;
  double y;
  double z;

  Point() : x(0), y(0), z(0){};

  Point(double x_, double y_, double z_) {
    x = x_;
    y = y_;
    z = z_;
  };
};

struct GeoConvertor {
  GeoConvertor(GPSPoint& ref_gps);

  void SetRefGps(GPSPoint& ref_gps);

  // Refer : https://en.wikipedia.org/wiki/Geographic_coordinate_conversion
  static void GpsToECEF(GPSPoint gps, Point& ecef);

  void GPSToENU(GPSPoint ref_gps, GPSPoint cur_gps, Point& enu);

  // Refer : https://en.wikipedia.org/wiki/Azimuthal_equidistant_projection
  void AEP(GPSPoint ref_gps, GPSPoint gps, Point& p);

  void QuickGPSToENU(GPSPoint ref_gps, GPSPoint gps, Point& p);

  void QuickENUToGPS(GPSPoint ref_gps, Point p, GPSPoint& gps);

  GPSPoint ref_gps;

  Point ref_ecef;
  double ref_lat_rad;
  double ref_lon_rad;

  double ref_sin_lat;
  double ref_cos_lat;

  double ref_sin_lon;
  double ref_cos_lon;
  double e2;
  double r_m;
  double r_n;
  double r_n_lat;
};

double Distance(Point& p1, Point& p2) {
  return std::sqrt((p1.x - p2.x) * (p1.x - p2.x) +
                   (p1.y - p2.y) * (p1.y - p2.y) +
                   (p1.z - p2.z) * (p1.z - p2.z));
}

}  // namespace colmap

#endif  // COLMAP_SRC_BASE_GPS_H_
