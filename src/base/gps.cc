#include "base/gps.h"

#include "util/math.h"

namespace colmap {

GeoConvertor::GeoConvertor(GPSPoint& ref_gps_in) {
  SetRefGps(ref_gps_in);

  GpsToECEF(ref_gps, ref_ecef);
}

void GeoConvertor::SetRefGps(GPSPoint& ref_gps_in) {
  ref_gps = ref_gps_in;

  double ref_lat_rad = AngelToRadian(ref_gps.latitude);
  double ref_lon_rad = AngelToRadian(ref_gps.longitude);
  double ref_sin_lat = std::sin(ref_lat_rad);
  double ref_cos_lat = std::cos(ref_lat_rad);
  double ref_sin_lon = std::sin(ref_lon_rad);
  double ref_cos_lon = std::cos(ref_lon_rad);

  double e2 = EARTH_F_r * (2 - EARTH_F_r);
  double r_m =
      EARTH_A * (1 - e2) / std::pow((1 - e2 * ref_sin_lat * ref_sin_lat), 1.5) +
      ref_gps.altitude;
  double r_n = EARTH_A / std::sqrt(1 - e2 * ref_sin_lat * ref_sin_lat) +
               ref_gps.altitude;
  double r_n_lat = r_n * ref_cos_lat;
}

// Refer : https://en.wikipedia.org/wiki/Azimuthal_equidistant_projection
void GeoConvertor::AEP(GPSPoint ref_gps, GPSPoint gps, Point& p) {
  double lat_rad = AngelToRadian(gps.latitude);
  double lon_rad = AngelToRadian(gps.longitude);
  double sin_lat = std::sin(lat_rad);
  double cos_lat = std::cos(lat_rad);
  double cos_d_lon = std::cos(lon_rad - ref_lon_rad);

  double arg = ref_sin_lat * sin_lat + ref_cos_lat * cos_lat * cos_d_lon;
  arg = std::min(arg, 1.0);
  arg = std::max(arg, -1.0);
  double c = std::acos(arg);
  double k = (fabs(c) < DBL_EPSILON) ? 1.0 : (c / std::sin(c));
  p.y = k * (ref_cos_lat * sin_lat - ref_sin_lat * cos_lat * cos_d_lon) * r_m;
  p.x = k * cos_lat * std::sin(lon_rad - ref_lon_rad) * r_n;
  p.z = gps.altitude - ref_gps.altitude;
}

void GeoConvertor::QuickGPSToENU(GPSPoint ref_gps, GPSPoint gps, Point& p) {
  double lat_rad = AngelToRadian(gps.latitude);
  double lon_rad = AngelToRadian(gps.longitude);
  p.x = (lon_rad - ref_lon_rad) * r_n_lat;
  p.y = (lat_rad - ref_lat_rad) * r_m;
  p.z = gps.altitude - ref_gps.altitude;
}

void GeoConvertor::QuickENUToGPS(GPSPoint ref_gps, Point p, GPSPoint& gps) {
  gps.longitude = RadianToAngel(ref_lon_rad + p.x / r_n_lat);
  gps.latitude = RadianToAngel(ref_lat_rad + p.y / r_m);
  gps.altitude = p.z + ref_gps.altitude;
}

void GeoConvertor::GpsToECEF(GPSPoint gps, Point& ecef) {
  double lat_rad = AngelToRadian(gps.latitude);
  double lon_rad = AngelToRadian(gps.longitude);

  double cos_lat = std::cos(lat_rad);
  double cos_lon = std::cos(lon_rad);
  double sin_lat = std::sin(lat_rad);
  double sin_lon = std::sin(lon_rad);

  double e2_cur = EARTH_F_r * (2 - EARTH_F_r);
  double N = EARTH_A / (std::sqrt(1.0 - (e2_cur * sin_lat * sin_lat)));

  ecef.x = (N + gps.altitude) * cos_lat * cos_lon;
  ecef.y = (N + gps.altitude) * cos_lat * sin_lon;
  ecef.z = ((1 - e2_cur) * N + gps.altitude) * sin_lat;
}

// Refer : https://en.wikipedia.org/wiki/Geographic_coordinate_conversion
void GeoConvertor::GPSToENU(GPSPoint ref_gps, GPSPoint cur_gps, Point& enu) {
  Point cur_ecef;
  GpsToECEF(cur_gps, cur_ecef);

  double delta_x = cur_ecef.x - ref_ecef.x;
  double delta_y = cur_ecef.y - ref_ecef.y;
  double delta_z = cur_ecef.z - ref_ecef.z;

  enu.x = delta_x * (-ref_sin_lon) + delta_y * (ref_cos_lon);
  enu.y = delta_x * (-ref_sin_lat) * (ref_cos_lon) +
          delta_y * (-ref_sin_lat) * (ref_sin_lon) + delta_z * ref_cos_lat;
  enu.z = delta_x * (ref_cos_lat) * (ref_cos_lon) +
          delta_y * (ref_cos_lat) * (ref_sin_lon) + delta_z * ref_sin_lat;
}

}  // namespace colmap
