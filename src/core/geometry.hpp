#pragma once

#include <cstdint>
#include <numeric>
#include <vector>

#include <opencv2/opencv.hpp>

#include "meta.hpp"
#include "ratio.hpp"

namespace cv {
template <class T, int N>
inline bool operator<(const Vec<T, N> &a, const Vec<T, N> &b) {
  for (int i = 0; i < N; i++) {
    if (a[i] > b[i]) {
      return false;
    } else if (a[i] < b[i]) {
      return true;
    }
  }
  return false;
}
}

namespace pano {
namespace core {

// vectors/points
template <class T, int N> using Vec = cv::Vec<T, N>;
using Vec2 = Vec<double, 2>;
using Vec3 = Vec<double, 3>;
using Vec4 = Vec<double, 4>;
using Vec2f = Vec<float, 2>;
using Vec3f = Vec<float, 3>;
using Vec4f = Vec<float, 4>;
using Vec2i = Vec<int, 2>;
using Vec3i = Vec<int, 3>;
using Vec4i = Vec<int, 4>;
using Vec2ub = Vec<uint8_t, 2>;
using Vec3ub = Vec<uint8_t, 3>;
using Vec4ub = Vec<uint8_t, 4>;
template <class T, int N> using Point = cv::Vec<T, N>;
using Point2 = Point<double, 2>;
using Point3 = Point<double, 3>;
using Point4 = Point<double, 4>;
using Point2f = Point<float, 2>;
using Point3f = Point<float, 3>;
using Point4f = Point<float, 4>;
using Point2i = Point<int, 2>;
using Point3i = Point<int, 3>;
using Point4i = Point<int, 4>;

using Vec5 = Vec<double, 5>;
using Vec5f = Vec<float, 5>;
using Vec7 = Vec<double, 7>;
using Vec7f = Vec<float, 7>;

template <class T, int N> struct IsNotContainerByHand<Point<T, N>> : yes {};

// matrix
template <class T, int M, int N> using Mat = cv::Matx<T, M, N>;
using Mat3 = Mat<double, 3, 3>;
using Mat4 = Mat<double, 4, 4>;
using Mat3f = Mat<float, 3, 3>;
using Mat4f = Mat<float, 4, 4>;

template <class T, int M, int N>
struct IsNotContainerByHand<Mat<T, M, N>> : yes {};

using cv::norm;
template <class T> inline T normalize(const T &d) { return d / norm(d); }

template <int N = 3, class T = double> inline const Point<T, N> &Origin() {
  static const Point<T, N> _origin;
  return _origin;
}

template <int N = 3, class T = double> inline const Vec<T, N> &X() {
  static const Vec<T, N> _x(1);
  return _x;
}

template <int N = 3, class T = double> inline const Vec<T, N> &Y() {
  static const Vec<T, N> _y(0, 1);
  return _y;
}

template <int N = 3, class T = double> inline const Vec<T, N> &Z() {
  static const Vec<T, N> _z(0, 0, 1);
  return _z;
}

template <class To, class From, int N>
inline Vec<To, N> ecast(const Vec<From, N> &v) {
  return v;
}

template <class To, class From, int N>
inline std::vector<Vec<To, N>> ecast(const std::vector<Vec<From, N>> &v) {
  std::vector<Vec<To, N>> result(v.size());
  for (int i = 0; i < v.size(); i++) {
    result[i] = v[i];
  }
  return result;
}

template <class T, int M, int N>
inline Vec<T, M + N> cat(const Mat<T, M, 1> &a, const Mat<T, N, 1> &b) {
  Vec<T, M + N> ab;
  std::copy(a.val, a.val + M, ab.val);
  std::copy(b.val, b.val + N, ab.val + M);
  return ab;
}

template <class T, int M>
inline Vec<T, M + 1> cat(const Mat<T, M, 1> &a, const T &b) {
  Vec<T, M + 1> ab;
  std::copy(a.val, a.val + M, ab.val);
  ab[M] = b;
  return ab;
}

template <class T, int M>
inline Vec<T, M + 1> cat(const T &a, const Mat<T, M, 1> &b) {
  Vec<T, M + 1> ab;
  ab.val[0] = a;
  std::copy(b.val, b.val + M, ab.val + 1);
  return ab;
}

// homogeneous point
template <class T, int N> using HPoint = Ratio<Point<T, N>, T>;
template <class T, int N>
Vec<T, N + 1> VectorFromHPoint(const HPoint<T, N> &p, const T &scale = 1.0) {
  Vec<T, N + 1> v;
  std::copy(p.numerator.val, p.numerator.val + N, v.val);
  v[N] = p.denominator * scale;
  return v;
}
template <class T, int N>
HPoint<T, N - 1> HPointFromVector(const Vec<T, N> &v, const T &scale = 1.0) {
  HPoint<T, N - 1> hp;
  std::copy(v.val, v.val + N - 1, hp.numerator.val);
  hp.denominator = v[N - 1] / scale;
  return hp;
}
template <class T, int N> inline Ratio<T, T> norm(const HPoint<T, N> &p) {
  return Ratio<T, T>(norm(p.numerator), p.denominator);
}
template <class T, int N>
inline Ratio<T, T> dot(const HPoint<T, N> &a, const HPoint<T, N> &b) {
  return Ratio<T, T>(a.numerator.dot(b.numerator),
                     a.denominator * b.denominator);
}
using HPoint2 = HPoint<double, 2>;
using HPoint3 = HPoint<double, 3>;
using HPoint4 = HPoint<double, 4>;

// matrix transform
template <class T> Mat<T, 3, 3> MakeMat3Rotate(const Vec<T, 3> &axis, T angle) {
  auto a = core::normalize(axis);
  double l = a[0], m = a[1], n = a[2];
  double cosv = cos(angle), sinv = sin(angle);
  return Mat<T, 3, 3>(l * l * (1 - cosv) + cosv, m * l * (1 - cosv) - n * sinv,
                      n * l * (1 - cosv) + m * sinv,
                      l * m * (1 - cosv) + n * sinv, m * m * (1 - cosv) + cosv,
                      n * m * (1 - cosv) - l * sinv,
                      l * n * (1 - cosv) - m * sinv,
                      m * n * (1 - cosv) + l * sinv, n * n * (1 - cosv) + cosv);
}

template <class T> Mat<T, 4, 4> MakeMat4Rotate(const Vec<T, 3> &axis, T angle) {
  auto a = core::normalize(axis);
  double l = a[0], m = a[1], n = a[2];
  double cosv = cos(angle), sinv = sin(angle);
  return Mat<T, 4, 4>(
      l * l * (1 - cosv) + cosv, m * l * (1 - cosv) - n * sinv,
      n * l * (1 - cosv) + m * sinv, 0, l * m * (1 - cosv) + n * sinv,
      m * m * (1 - cosv) + cosv, n * m * (1 - cosv) - l * sinv, 0,
      l * n * (1 - cosv) - m * sinv, m * n * (1 - cosv) + l * sinv,
      n * n * (1 - cosv) + cosv, 0, 0, 0, 0, 1);
}

template <class T> Mat<T, 4, 4> MakeMat4Translate(const Vec<T, 3> &trans) {
  return Mat<T, 4, 4>(1, 0, 0, trans[0], 0, 1, 0, trans[1], 0, 0, 1, trans[2],
                      0, 0, 0, 1);
}

template <class T> Mat<T, 4, 4> MakeMat4Scale(const T &scale) {
  return Mat<T, 4, 4>(scale, 0, 0, 0, 0, scale, 0, 0, 0, 0, scale, 0, 0, 0, 0,
                      1);
}

template <class T> Mat<T, 4, 4> MakeMat4Scale(T sx, T sy, T sz) {
  return Mat<T, 4, 4>(sx, 0, 0, 0, 0, sy, 0, 0, 0, 0, sz, 0, 0, 0, 0, 1);
}

template <class T>
Mat<T, 4, 4>
MakeMat4LocalToWorld(const Vec<T, 3> &localx, const Vec<T, 3> &localy,
                     const Vec<T, 3> &localz, const Point<T, 3> &localo) {
  return Mat<T, 4, 4>(localx[0], localy[0], localz[0], localo[0], localx[1],
                      localy[1], localz[1], localo[1], localx[2], localy[2],
                      localz[2], localo[2], 0, 0, 0, 1)
      .t();
}

// camera functions with matrix
// make a lookat view matrix
template <class T>
Mat<T, 4, 4> MakeMat4LookAt(const Vec<T, 3> &eye, const Vec<T, 3> &center,
                            const Vec<T, 3> &up) {
  Vec<T, 3> zaxis = (center - eye);
  zaxis /= core::norm(zaxis);
  Vec<T, 3> xaxis = up.cross(zaxis);
  xaxis /= core::norm(xaxis);
  Vec<T, 3> yaxis = zaxis.cross(xaxis);
  Mat<T, 4, 4> m(-xaxis(0), yaxis(0), -zaxis(0), 0, -xaxis(1), yaxis(1),
                 -zaxis(1), 0, -xaxis(2), yaxis(2), -zaxis(2), 0,
                 xaxis.dot(eye), -yaxis.dot(eye), zaxis.dot(eye), 1);
  return m.t();
}

// make a perspective projection matrix
template <class T>
Mat<T, 4, 4> MakeMat4Perspective(const T &fovyRadians, const T &aspect,
                                 const T &nearZ, const T &farZ) {
  T cotan = T(1.0) / std::tan(fovyRadians / 2.0);
  Mat<T, 4, 4> m(cotan / aspect, 0, 0, 0, 0, cotan, 0, 0, 0, 0,
                 (farZ + nearZ) / (nearZ - farZ), -1, 0, 0,
                 (2 * farZ * nearZ) / (nearZ - farZ), 0);
  return m.t();
}

template <class T>
Mat<T, 4, 4> MakeMat4Perspective(const T &fx, const T &fy, const T &cx,
                                 const T &cy, const T &nearZ, const T &farZ) {
  Mat<T, 4, 4> m(fx / cx, 0, 0, 0, 0, fy / cy, 0, 0, 0, 0,
                 (farZ + nearZ) / (nearZ - farZ), -1, 0, 0,
                 (2 * farZ * nearZ) / (nearZ - farZ), 0);
  return m.t();
}

// geographic coordinate
template <class T> inline double Argument(const Vec<T, 2> &d) {
  return std::atan2(d(1), d(0));
}
template <class T> inline Vec<T, 2> Direction(double angle) {
  return Vec<T, 2>(cos(angle), sin(angle));
}
struct GeoCoord {
  inline explicit GeoCoord(double longi = 0.0, double lati = 0.0)
      : longitude(longi), latitude(lati) {}
  template <class T>
  inline GeoCoord(const Vec<T, 3> &d)
      : longitude(std::atan2(d(1), d(0))),
        latitude(std::atan(d(2) / std::sqrt((d(1) * d(1)) + (d(0) * d(0))))) {}
  template <class T = double> inline Vec<T, 3> toVector() const {
    return Vec<T, 3>(static_cast<T>(cos(longitude) * cos(latitude)),
                     static_cast<T>(sin(longitude) * cos(latitude)),
                     static_cast<T>(sin(latitude)));
  }
  double longitude; // - M_PI ~ + M_PI
  double latitude;  // - M_PI_2 ~ + M_PI_2
};
inline bool operator==(const GeoCoord &a, const GeoCoord &b) {
  return a.longitude == b.longitude && a.latitude == b.latitude;
}
template <class Archive> inline void serialize(Archive &ar, GeoCoord &gc) {
  ar(gc.longitude, gc.latitude);
}

// key point
using KeyPoint = cv::KeyPoint;

// size
using Size = cv::Size2f;
using Sizei = cv::Size2i;

// infinite line
template <class T, int N> struct Ray {
  inline Ray() {}
  inline Ray(const Point<T, N> &a, const Vec<T, N> &d)
      : anchor(a), direction(d) {}
  Point<T, N> anchor;
  Vec<T, N> direction;
};
template <class T, int N>
inline bool operator==(const Ray<T, N> &a, const Ray<T, N> &b) {
  return a.anchor == b.anchor && a.direction == b.direction;
}
template <class Archive, class T, int N>
inline void serialize(Archive &ar, Ray<T, N> &p) {
  ar(p.anchor, p.direction);
}
using Ray2 = Ray<double, 2>;
using Ray3 = Ray<double, 3>;
template <class T> inline Vec<T, 3> GetCoeffs(const Ray<T, 2> &line) {
  return Vec<T, 3>{line.direction[1], -line.direction[0],
                   -(line.direction[1] * line.anchor[0] -
                     line.direction[0] * line.anchor[1])};
}
template <class T> inline Ray<T, 2> Ray2FromCoeffs(const Vec<T, 3> &c) {
  T d = c[0] * c[0] + c[1] * c[1];
  Point<T, 2> anchor(-c[2] * c[0] / d, -c[2] * c[1] / d);
  Vec<T, 2> dir(c[1], -c[0]);
  return Ray<T, 2>(anchor, dir);
}

// plane
template <class T, int N> struct Plane {
  inline Plane() {}
  inline Plane(const Point<T, N> &a, const Vec<T, N> &n)
      : anchor(a), normal(n) {}
  inline Point<T, N> root() const {
    return normal * (anchor.dot(normal)) / norm(normal) / norm(normal);
  }
  inline T signedDistanceTo(const Point<T, N> &p) const {
    return (p - anchor).dot(normalize(normal));
  }
  inline T distanceTo(const Point<T, N> &p) const {
    return abs((p - anchor).dot(normalize(normal)));
  }
  Point<T, N> anchor;
  Vec<T, N> normal;
};
template <class T, int N>
inline bool operator==(const Plane<T, N> &a, const Plane<T, N> &b) {
  return a.anchor == b.anchor && a.normal == b.normal;
}
template <class Archive, class T, int N>
inline void serialize(Archive &ar, Plane<T, N> &p) {
  ar(p.anchor, p.normal);
}
using Plane3 = Plane<double, 3>;
template <class T>
inline Plane<T, 3> Plane3From3Points(const Point<T, 3> &a, const Point<T, 3> &b,
                                     const Point<T, 3> &c) {
  return Plane<T, 3>(a, normalize((b - a).cross(c - a)));
}
// ax + by + cz = 1
template <class T> inline Plane<T, 3> Plane3FromEquation(T a, T b, T c) {
  T k = (a * a + b * b + c * c);
  return Plane<T, 3>(Point<T, 3>(a, b, c) / k, normalize(Vec<T, 3>(a, b, c)));
}
template <class T> inline Plane<T, 3> Plane3FromEquation(const Vec<T, 3> &equ) {
  return Plane3FromEquation(equ[0], equ[1], equ[2]);
}
template <class T> inline Vec<T, 3> Plane3ToEquation(const Plane<T, 3> &p) {
  auto dotv = p.anchor.dot(p.normal);
  assert(dotv != 0.0);
  return p.normal / dotv;
}

// line
template <class T, int N> struct Line {
  static const int Dimension = N;
  inline Line() {}
  inline Line(const Point<T, N> &f, const Point<T, N> &s)
      : first(f), second(s) {}
  inline Point<T, N> center() const { return (first + second) / 2.0; }
  inline T length() const { return norm(first - second); }
  inline Vec<T, N> direction() const { return second - first; }
  inline Line reversed() const { return Line(second, first); }
  inline Ray<T, N> ray() const { return Ray<T, N>{first, second - first}; }
  inline Line &operator*=(const T &factor) {
    first *= factor;
    second *= factor;
    return *this;
  }
  inline Line &operator/=(const T &factor) {
    first /= factor;
    second /= factor;
    return *this;
  }
  inline Line &operator+=(const Vec<T, N> &trans) {
    first += trans;
    second += trans;
    return *this;
  }
  inline Line &operator-=(const Vec<T, N> &trans) {
    first -= trans;
    second -= trans;
    return *this;
  }
  Point<T, N> first, second;
};
template <class T, int N>
inline bool operator==(const Line<T, N> &a, const Line<T, N> &b) {
  return a.first == b.first && a.second == b.second;
}
template <class Archive, class T, int N>
inline void serialize(Archive &ar, Line<T, N> &l) {
  ar(l.first, l.second);
}
template <class T, int N>
inline Line<T, N> operator*(const Line<T, N> &line, const T &factor) {
  return Line<T, N>(line.first * factor, line.second * factor);
}
template <class T, int N>
inline Line<T, N> operator/(const Line<T, N> &line, const T &factor) {
  return Line<T, N>(line.first / factor, line.second / factor);
}
template <class T, int N>
inline Line<T, N> operator*(const T &factor, const Line<T, N> &line) {
  return Line<T, N>(line.first * factor, line.second * factor);
}
template <class T, int N>
inline Line<T, N> operator+(const Line<T, N> &line, const Vec<T, N> &trans) {
  return Line<T, N>(line.first + trans, line.second + trans);
}
template <class T, int N>
inline Line<T, N> operator-(const Line<T, N> &line, const Vec<T, N> &trans) {
  return Line<T, N>(line.first - trans, line.second - trans);
}
using Line2 = Line<double, 2>;
using Line3 = Line<double, 3>;

template <class T, int N> inline Line<T, N> normalize(const Line<T, N> &line) {
  return Line<T, N>(normalize(line.first), normalize(line.second));
}

template <class T> struct IsLine : std::false_type {};

template <class T, int N> struct IsLine<Line<T, N>> : std::true_type {};

// position on line/infline
template <class T, int N> struct PositionOnLine {
  inline PositionOnLine() {}
  inline PositionOnLine(const Line<T, N> &line, const T &r)
      : ratio(r), position(line.first + (line.second - line.first) * ratio) {}
  inline PositionOnLine(const Ray<T, N> &line, const T &r)
      : ratio(r), position(line.anchor + line.direction * ratio) {}
  T ratio; // [0 ~ 1] on line, or [-inf, +inf] on infinite line
  // position = line.first + (line.second - line.fist) * ratio
  // or       = line.anchor + line.direction * ratio
  Point<T, N> position;
};
template <class T, int N>
inline bool operator==(const PositionOnLine<T, N> &a,
                       const PositionOnLine<T, N> &b) {
  return a.ratio == b.ratio && a.position == b.position;
}
template <class Archive, class T, int N>
inline void serialize(Archive &ar, PositionOnLine<T, N> &p) {
  ar(p.ratio, p.position);
}
using PositionOnLine2 = PositionOnLine<double, 2>;
using PositionOnLine3 = PositionOnLine<double, 3>;

// chain
template <class T, int N> struct Chain {
  std::vector<Point<T, N>> points;
  bool closed;

  Chain() : closed(true) {}
  explicit Chain(const std::vector<Point<T, N>> &ps, bool c = true)
      : points(ps), closed(c) {}
  explicit Chain(std::vector<Point<T, N>> &&ps, bool c = true)
      : points(std::move(ps)), closed(c) {}

  const Point<T, N> &at(size_t i) const { return points.at(i % points.size()); }
  const Point<T, N> &prev(size_t i) const {
    return points.at((i - 1 + points.size()) % points.size());
  }
  const Point<T, N> &next(size_t i) const {
    return points.at((i + 1) % points.size());
  }

  Line<T, N> edge(size_t i) const { return Line<T, N>(at(i), next(i)); }

  const Point<T, N> &operator[](size_t i) const { return points[i]; }
  Point<T, N> &operator[](size_t i) { return points[i]; }

  size_t size() const { return points.size(); }
  void clear() { points.clear(); }
  bool empty() const { return points.empty(); }

  void append(const Point<T, N> &p) { points.push_back(p); }
  void insert(size_t i, const Point<T, N> &p) {
    points.insert(points.begin() + (i % points.size()), p);
  }

  T length() const {
    if (points.empty())
      return 0;
    T len = 0;
    for (int i = 0; i + 1 < points.size(); i++) {
      len += norm(points[i] - points[i + 1]);
    }
    if (closed) {
      len += norm(points.front() - points.back());
    }
    return len;
  }

  template <class K> operator Chain<K, N>() const {
    Chain<K, N> c;
    c.points.resize(points.size());
    for (int i = 0; i < points.size(); i++) {
      c.points[i] = points[i];
    }
    c.closed = closed;
    return c;
  }

  template <class FunT> void fixedStepSample(double stepLen, FunT &&fun) const {
    if (empty()) {
      return;
    }
    fun(points.front());
    auto last = points.front();
    for (int i = 1; i < points.size(); i++) {
      double dist = Distance(last, points[i]);
      auto dir = normalize(points[i] - last);
      while (dist >= stepLen) {
        auto p = last + dir * stepLen;
        fun(p);
        dist -= stepLen;
        last = p;
      }
    }
  }
};

template <class T, int N>
inline bool operator==(const Chain<T, N> &a, const Chain<T, N> &b) {
  return a.points == b.points && a.closed == b.closed;
}
template <class Archive, class T, int N>
inline void serialize(Archive &ar, Chain<T, N> &l) {
  ar(l.points, l.closed);
}
template <class T, int N> inline Chain<T, N> normalize(const Chain<T, N> &c) {
  auto r = c;
  for (auto &p : c.points) {
    p = normalize(p);
  }
  return r;
}
using Chain2i = Chain<int, 2>;
using Chain2 = Chain<double, 2>;
using Chain3 = Chain<double, 3>;

// sphere
template <class T, int N> struct Sphere {
  Point<T, N> center;
  T radius;
};
template <class T, int N>
inline bool operator==(const Sphere<T, N> &a, const Sphere<T, N> &b) {
  return a.center == b.center && a.radius == b.radius;
}
template <class Archive, class T, int N>
inline void serialize(Archive &ar, Sphere<T, N> &s) {
  ar(s.center, s.radius);
}

using Circle = Sphere<double, 2>;
using Sphere2 = Sphere<double, 2>;
using Sphere3 = Sphere<double, 3>;

template <int N = 3, class T = double> const Sphere<T, N> &UnitSphere() {
  static const Sphere<T, N> _us = {Point<T, N>(), static_cast<T>(1.0)};
  return _us;
}

namespace {
template <class T, int N>
Vec<T, N> AllMinOf(const Vec<T, N> &v1, const Vec<T, N> &v2) {
  Vec<T, N> v;
  for (int i = 0; i < N; i++)
    v[i] = v1[i] < v2[i] ? v1[i] : v2[i];
  return v;
}
template <class T, int N>
Vec<T, N> AllMaxOf(const Vec<T, N> &v1, const Vec<T, N> &v2) {
  Vec<T, N> v;
  for (int i = 0; i < N; i++)
    v[i] = v1[i] < v2[i] ? v2[i] : v1[i];
  return v;
}
}

// box
template <class T, int N> struct Box {
  using Type = T;
  static const int Dimension = N;

  Point<T, N> minCorner, maxCorner;
  bool isNull;

  Box() : isNull(true) {}
  Box(const Point<T, N> &c1, const Point<T, N> &c2)
      : minCorner(AllMinOf(c1, c2)), maxCorner(AllMaxOf(c1, c2)),
        isNull(false) {}
  template <class = std::enable_if_t<N == 1>>
  Box(const T &c1, const T &c2)
      : minCorner(std::min(c1, c2)), maxCorner(std::max(c1, c2)),
        isNull(false) {}

  Point<T, N> corner(std::initializer_list<bool> isMaxes) const {
    Point<T, N> c;
    auto it = isMaxes.begin();
    for (int i = 0; i < N; i++) {
      c[i] = (*it) ? maxCorner[i] : minCorner[i];
      ++it;
    }
    return c;
  }

  Vec<T, N> size() const { return maxCorner - minCorner; }
  T size(size_t i) const { return maxCorner[i] - minCorner[i]; }
  T volume() const {
    T v = 1.0;
    for (int i = 0; i < N; i++) {
      v *= (maxCorner[i] - minCorner[i]);
    }
    return v;
  }

  Point<T, N> center() const { return (maxCorner + minCorner) * (0.5); }
  Sphere<T, N> outerSphere() const {
    return Sphere<T, N>{center(),
                        static_cast<T>(norm(maxCorner - minCorner) / 2.0)};
  }
  Box &expand(const Vec<T, N> &s) {
    minCorner -= s;
    maxCorner += s;
    return *this;
  }
  Box &expand(const T &s) {
    for (int i = 0; i < N; i++) {
      minCorner[i] -= s;
      maxCorner[i] += s;
    }
    return *this;
  }

  bool contains(const Point<T, N> &p) const {
    if (isNull)
      return false;
    for (int i = 0; i < N; i++) {
      if (minCorner[i] > p[i] || maxCorner[i] < p[i])
        return false;
    }
    return true;
  }
  bool contains(const Box &b) const {
    return b.isNull ? true : contains(b.minCorner) && contains(b.maxCorner);
  }
  bool operator!=(const Box &b) const { return !(*this == b); }

  Box &operator|=(const Box &b) {
    if (isNull) {
      *this = b;
      return *this;
    }
    if (b.isNull)
      return *this;
    minCorner = AllMinOf(minCorner, b.minCorner);
    maxCorner = AllMaxOf(maxCorner, b.maxCorner);
    return *this;
  }
};
template <class T, int N>
inline bool operator==(const Box<T, N> &a, const Box<T, N> &b) {
  return a.isNull ? b.isNull : (!b.isNull && a.minCorner == b.minCorner &&
                                a.maxCorner == b.maxCorner);
}
template <class Archive, class T, int N>
inline void serialize(Archive &ar, Box<T, N> &b) {
  ar(b.isNull, b.minCorner, b.maxCorner);
}
template <class T, int N>
inline Box<T, N> operator|(const Box<T, N> &b1, const Box<T, N> &b2) {
  Box<T, N> b12 = b1;
  return b12 |= b2;
}
using Box1 = Box<double, 1>;
using Box2 = Box<double, 2>;
using Box3 = Box<double, 3>;

template <int N = 3, class T = double> const Box<T, N> &UnitBox() {
  static const Box<T, N> _ub(Point<T, N>(), Point<T, N>::ones());
  return _ub;
}

template <class T> struct IsBox : std::false_type {};

template <class T, int N> struct IsBox<Box<T, N>> : std::true_type {};

// polygon
template <class T, int N> struct Polygon {
  std::vector<Point<T, N>> corners;
  Vec<T, N> normal;

  Polygon() {}
  Polygon(const std::vector<Point<T, N>> &cs, const Vec<T, N> &n)
      : corners(cs), normal(n) {}
  Polygon(std::vector<Point<T, N>> &&cs, const Vec<T, N> &n)
      : corners(std::move(cs)), normal(n) {}
  Polygon(const Chain<T, N> &c)
      : corners(c.points),
        normal(normalize((corners.at(0) - corners.at(1))
                             .cross(corners.at(1) - corners.at(2)))) {}
  Polygon(Chain<T, N> &&c)
      : corners(std::move(c.points)),
        normal(normalize((corners.at(0) - corners.at(1))
                             .cross(corners.at(1) - corners.at(2)))) {}
  inline Plane<T, N> plane() const {
    return Plane<T, N>(corners.front(), normal);
  }
  inline Chain<T, N> boundary() const { return Chain<T, N>{corners, true}; }

  template <class K, std::enable_if_t<!std::is_same<K, T>::value>>
  operator Polygon<K, N>() const {
    std::vector<Point<K, N>> ps(corners.size());
    for (int i = 0; i < corners.size(); i++)
      ps[i] = corners[i];
    return Polygon<K, N>(std::move(ps), normal);
  }

  template <
      class = std::enable_if_t<std::is_floating_point<T>::value && N == 3>>
  inline T area() const {
    return Area(*this);
  }
};

template <class Archive, class T, int N>
inline void serialize(Archive &ar, Polygon<T, N> &p) {
  ar(p.corners, p.normal);
}

using Polygon2 = Polygon<double, 2>;
using Polygon2f = Polygon<float, 2>;
using Polygon3 = Polygon<double, 3>;
using Polygon3f = Polygon<float, 3>;

template <class T>
inline Polygon<T, 3> MakeTriangle(const Point<T, 3> &p1, const Point<T, 3> &p2,
                                  const Point<T, 3> &p3) {
  return Polygon<T, 3>{{p1, p2, p3}, (p2 - p1).cross(p3 - p1)};
}

float Area(const Polygon3f &polygon);
double Area(const Polygon3 &polygon);
// long double Area(const Polygon<long double, 3> & polygon);

double PointTest(const Polygon2f &poly, const Point2f &p);
bool Contains(const Polygon2f &poly, const Point2f &p);

// layered polygons
template <class T, int N> struct LayeredShape {
  std::vector<std::vector<Point<T, N>>> layers;
  Vec<T, N> normal;
  inline Polygon<T, N> layer(size_t i) const {
    return Polygon<T, N>(layers.at(i), normal);
  }
  inline size_t size() const { return layers.size(); }
  inline bool empty() const { return layers.empty(); }
};
template <class Archive, class T, int N>
inline void serialize(Archive &ar, LayeredShape<T, N> &p) {
  ar(p.layers, p.normal);
}
using LayeredShape3 = LayeredShape<double, 3>;

// transformed in 3d space
template <class T, class E = double> struct TransformedIn3D {
  T component;
  Mat<E, 4, 4> mat4;

  TransformedIn3D() : mat4(Mat<E, 4, 4>::eye()) {}
  TransformedIn3D(const T &c, const Mat<E, 4, 4> &m = Mat<E, 4, 4>::eye())
      : component(c), mat4(m) {}
  TransformedIn3D(T &&c, const Mat<E, 4, 4> &m = Mat<E, 4, 4>::eye())
      : component(std::move(c)), mat4(m) {}

  TransformedIn3D &translate(const Vec<E, 3> &t) {
    mat4 = MakeMat4Translate(t) * mat4;
    return *this;
  }
  TransformedIn3D &rotate(const Vec<E, 3> &axis, E angle) {
    mat4 = MakeMat4Rotate(axis, angle) * mat4;
    return *this;
  }
  TransformedIn3D &scale(E s) {
    mat4 = MakeMat4Scale(s) * mat4;
    return *this;
  }
  TransformedIn3D &scale(E sx, E sy, E sz) {
    mat4 = MakeMat4Scale(sx, sy, sz) * mat4;
    return *this;
  }
  TransformedIn3D &reset() {
    mat4 = Mat<E, 4, 4>::eye();
    return *this;
  }

  Point<E, 3> toWorld(const Point<E, 3> &p) const {
    Vec<E, 4> c = mat4 * cat(p, 1.0);
    return Point<E, 3>(c[0] / c[3], c[1] / c[3], c[2] / c[3]);
  }
  Point<E, 3> toLocal(const Point<E, 3> &p) const {
    Vec<E, 4> c = mat4 * cat(p, 1.0);
    Vec<E, 4> localc;
    bool solvable = cv::solve(mat4, c, localc);
    assert(solvable);
    return Point<E, 3>(localc[0] / localc[3], localc[1] / localc[3],
                       localc[2] / localc[3]);
  }

  template <class = std::enable_if_t<std::is_same<T, Line<E, 3>>::value>>
  Point<E, 3> first() const {
    return toWorld(component.first);
  }
  template <class = std::enable_if_t<std::is_same<T, Line<E, 3>>::value>>
  Point<E, 3> second() const {
    return toWorld(component.second);
  }
  template <class = std::enable_if_t<std::is_same<T, Box<E, 3>>::value>>
  Point<E, 3> corner(std::initializer_list<bool> isMaxes) const {
    return toWorld(component.corner(isMaxes));
  }
};

template <class E = double, class T>
TransformedIn3D<std::decay_t<T>, E> MakeTransformableIn3D(T &&c) {
  return TransformedIn3D<std::decay_t<T>, E>(std::forward<T>(c));
}

template <class E = double, class T>
TransformedIn3D<std::decay_t<T>, E>
AsInLocalCoordinates(T &&c, const Vec<E, 3> &x, const Vec<E, 3> &y,
                     const Vec<E, 3> &z, const Point<E, 3> &o) {
  return TransformedIn3D<std::decay_t<T>, E>(std::forward<T>(c),
                                             MakeMat4LocalToWorld(x, y, z, o));
}

template <class T, class E>
inline bool operator==(const TransformedIn3D<T, E> &a,
                       const TransformedIn3D<T, E> &b) {
  return a.component == b.component && a.mat4 == b.mat4;
}
template <class Archive, class T, class E>
inline void serialize(Archive &ar, TransformedIn3D<T, E> &b) {
  ar(b.component, b.mat4);
}
}
}

namespace std {

template <class T, int N, int M>
const T *begin(const pano::core::Mat<T, N, M> &v) {
  return v.val;
}

template <class T, int N, int M>
const T *end(const pano::core::Mat<T, N, M> &v) {
  return v.val + N * M;
}

template <class T, int N, int M>
const T *cbegin(const pano::core::Mat<T, N, M> &v) {
  return v.val;
}

template <class T, int N, int M>
const T *cend(const pano::core::Mat<T, N, M> &v) {
  return v.val + N * M;
}

template <class T, int N, int M> T *begin(pano::core::Mat<T, N, M> &v) {
  return v.val;
}

template <class T, int N, int M> T *end(pano::core::Mat<T, N, M> &v) {
  return v.val + N * M;
}
}
