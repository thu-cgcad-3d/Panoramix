#include "algorithms.hpp"
#include "geometry.hpp"
#include "utility.hpp"

namespace pano {
namespace core {

template <class T> T AreaTemplated(const Polygon<T, 3> &polygon) {
  T area = 0.0;
  Vec<T, 3> x, y;
  std::tie(x, y) = ProposeXYDirectionsFromZDirection(polygon.normal);
  TriangulatePolygon(polygon.corners.begin(), polygon.corners.end(),
                     [&](const Point<T, 3> &v) {
                       auto proj = v - v.dot(polygon.normal) * polygon.normal;
                       return Vec<T, 2>(proj.dot(x), proj.dot(y));
                     },
                     [&area](const Point<T, 3> &a, const Point<T, 3> &b,
                             const Point<T, 3> &c) {
                       double aa = AreaOfTriangle<T, 3>(a, b, c);
                       area += IsInfOrNaN(aa) ? 0.0 : aa;
                     });
  return area;
}

float Area(const Polygon3f &polygon) { return AreaTemplated(polygon); }
double Area(const Polygon3 &polygon) { return AreaTemplated(polygon); }

double PointTest(const Polygon2f &poly, const Point2f &p) {
  return cv::pointPolygonTest(poly.corners, p, true);
}

bool Contains(const Polygon2f &poly, const Point2f &p) {
  return cv::pointPolygonTest(poly.corners, p, false) >= 0;
}
}
}