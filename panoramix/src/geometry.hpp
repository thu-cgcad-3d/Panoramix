#pragma once

#include "math.hpp"

namespace pano {
namespace core {

// geographic coordinate
template <class T> double Argument(const Vec<T, 2> &d) {
  return std::atan2(d(1), d(0));
}
template <class T> Vec<T, 2> Direction(double angle) {
  return Vec<T, 2>(cos(angle), sin(angle));
}
struct GeoCoord {
  explicit GeoCoord(double longi = 0.0, double lati = 0.0)
      : longitude(longi), latitude(lati) {}
  template <class T>
  GeoCoord(const Vec<T, 3> &d)
      : longitude(std::atan2(d(1), d(0))),
        latitude(std::atan(d(2) / std::sqrt((d(1) * d(1)) + (d(0) * d(0))))) {}
  template <class T = double> Vec<T, 3> toVector() const {
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
template <class Archive> void serialize(Archive &ar, GeoCoord &gc) {
  ar(gc.longitude, gc.latitude);
}

// key point
using KeyPoint = cv::KeyPoint;

// size
template <class T> using Size_ = cv::Size_<T>;
using Size = Size_<float>;
using Sizei = Size_<int>;

// infinite line
template <class T, int N> struct Ray {
  Ray() {}
  Ray(const Point<T, N> &a, const Vec<T, N> &d) : anchor(a), direction(d) {}
  Point<T, N> anchor;
  Vec<T, N> direction;
};
template <class T, int N>
bool operator==(const Ray<T, N> &a, const Ray<T, N> &b) {
  return a.anchor == b.anchor && a.direction == b.direction;
}
template <class Archive, class T, int N>
void serialize(Archive &ar, Ray<T, N> &p) {
  ar(p.anchor, p.direction);
}
using Ray2 = Ray<double, 2>;
using Ray3 = Ray<double, 3>;
template <class T> Vec<T, 3> GetCoeffs(const Ray<T, 2> &line) {
  return Vec<T, 3>{line.direction[1], -line.direction[0],
                   -(line.direction[1] * line.anchor[0] -
                     line.direction[0] * line.anchor[1])};
}
template <class T> Ray<T, 2> Ray2FromCoeffs(const Vec<T, 3> &c) {
  T d = c[0] * c[0] + c[1] * c[1];
  Point<T, 2> anchor(-c[2] * c[0] / d, -c[2] * c[1] / d);
  Vec<T, 2> dir(c[1], -c[0]);
  return Ray<T, 2>(anchor, dir);
}

// plane
template <class T, int N> struct Plane {
  Plane() {}
  Plane(const Point<T, N> &a, const Vec<T, N> &n) : anchor(a), normal(n) {}
  Point<T, N> root() const {
    return normal * (anchor.dot(normal)) / norm(normal) / norm(normal);
  }
  T signedDistanceTo(const Point<T, N> &p) const {
    return (p - anchor).dot(normalize(normal));
  }
  T distanceTo(const Point<T, N> &p) const {
    return abs((p - anchor).dot(normalize(normal)));
  }
  Point<T, N> anchor;
  Vec<T, N> normal;
};
template <class T, int N>
bool operator==(const Plane<T, N> &a, const Plane<T, N> &b) {
  return a.anchor == b.anchor && a.normal == b.normal;
}
template <class Archive, class T, int N>
void serialize(Archive &ar, Plane<T, N> &p) {
  ar(p.anchor, p.normal);
}
using Plane3 = Plane<double, 3>;
template <class T>
Plane<T, 3> Plane3From3Points(const Point<T, 3> &a, const Point<T, 3> &b,
                              const Point<T, 3> &c) {
  return Plane<T, 3>(a, normalize((b - a).cross(c - a)));
}
// ax + by + cz = 1
template <class T> Plane<T, 3> Plane3FromEquation(T a, T b, T c) {
  T k = (a * a + b * b + c * c);
  return Plane<T, 3>(Point<T, 3>(a, b, c) / k, normalize(Vec<T, 3>(a, b, c)));
}
template <class T> Plane<T, 3> Plane3FromEquation(const Vec<T, 3> &equ) {
  return Plane3FromEquation(equ[0], equ[1], equ[2]);
}
template <class T> Vec<T, 3> Plane3ToEquation(const Plane<T, 3> &p) {
  auto dotv = p.anchor.dot(p.normal);
  assert(dotv != 0.0);
  return p.normal / dotv;
}

// line
template <class T, int N> struct Line {
  static const int Dimension = N;
  Line() {}
  Line(const Point<T, N> &f, const Point<T, N> &s) : first(f), second(s) {}
  Point<T, N> center() const { return (first + second) / 2.0; }
  T length() const { return norm(first - second); }
  Vec<T, N> direction() const { return second - first; }
  Line reversed() const { return Line(second, first); }
  Ray<T, N> ray() const { return Ray<T, N>{first, second - first}; }
  Line &operator*=(const T &factor) {
    first *= factor;
    second *= factor;
    return *this;
  }
  Line &operator/=(const T &factor) {
    first /= factor;
    second /= factor;
    return *this;
  }
  Line &operator+=(const Vec<T, N> &trans) {
    first += trans;
    second += trans;
    return *this;
  }
  Line &operator-=(const Vec<T, N> &trans) {
    first -= trans;
    second -= trans;
    return *this;
  }
  Point<T, N> first, second;
};
template <class T, int N>
bool operator==(const Line<T, N> &a, const Line<T, N> &b) {
  return a.first == b.first && a.second == b.second;
}
template <class Archive, class T, int N>
void serialize(Archive &ar, Line<T, N> &l) {
  ar(l.first, l.second);
}
template <class T, int N>
Line<T, N> operator*(const Line<T, N> &line, const T &factor) {
  return Line<T, N>(line.first * factor, line.second * factor);
}
template <class T, int N>
Line<T, N> operator/(const Line<T, N> &line, const T &factor) {
  return Line<T, N>(line.first / factor, line.second / factor);
}
template <class T, int N>
Line<T, N> operator*(const T &factor, const Line<T, N> &line) {
  return Line<T, N>(line.first * factor, line.second * factor);
}
template <class T, int N>
Line<T, N> operator+(const Line<T, N> &line, const Vec<T, N> &trans) {
  return Line<T, N>(line.first + trans, line.second + trans);
}
template <class T, int N>
Line<T, N> operator-(const Line<T, N> &line, const Vec<T, N> &trans) {
  return Line<T, N>(line.first - trans, line.second - trans);
}
using Line2 = Line<double, 2>;
using Line3 = Line<double, 3>;

template <class T, int N> Line<T, N> normalize(const Line<T, N> &line) {
  return Line<T, N>(normalize(line.first), normalize(line.second));
}

template <class T> struct IsLine : std::false_type {};

template <class T, int N> struct IsLine<Line<T, N>> : std::true_type {};

// position on line/infline
template <class T, int N> struct PositionOnLine {
  PositionOnLine() {}
  PositionOnLine(const Line<T, N> &line, const T &r)
      : ratio(r), position(line.first + (line.second - line.first) * ratio) {}
  PositionOnLine(const Ray<T, N> &line, const T &r)
      : ratio(r), position(line.anchor + line.direction * ratio) {}
  T ratio; // [0 ~ 1] on line, or [-inf, +inf] on infinite line
  // position = line.first + (line.second - line.fist) * ratio
  // or       = line.anchor + line.direction * ratio
  Point<T, N> position;
};
template <class T, int N>
bool operator==(const PositionOnLine<T, N> &a, const PositionOnLine<T, N> &b) {
  return a.ratio == b.ratio && a.position == b.position;
}
template <class Archive, class T, int N>
void serialize(Archive &ar, PositionOnLine<T, N> &p) {
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
bool operator==(const Chain<T, N> &a, const Chain<T, N> &b) {
  return a.points == b.points && a.closed == b.closed;
}
template <class Archive, class T, int N>
void serialize(Archive &ar, Chain<T, N> &l) {
  ar(l.points, l.closed);
}
template <class T, int N> Chain<T, N> normalize(const Chain<T, N> &c) {
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
bool operator==(const Sphere<T, N> &a, const Sphere<T, N> &b) {
  return a.center == b.center && a.radius == b.radius;
}
template <class Archive, class T, int N>
void serialize(Archive &ar, Sphere<T, N> &s) {
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
bool operator==(const Box<T, N> &a, const Box<T, N> &b) {
  return a.isNull ? b.isNull : (!b.isNull && a.minCorner == b.minCorner &&
                                a.maxCorner == b.maxCorner);
}
template <class Archive, class T, int N>
void serialize(Archive &ar, Box<T, N> &b) {
  ar(b.isNull, b.minCorner, b.maxCorner);
}
template <class T, int N>
Box<T, N> operator|(const Box<T, N> &b1, const Box<T, N> &b2) {
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
  Plane<T, N> plane() const { return Plane<T, N>(corners.front(), normal); }
  Chain<T, N> boundary() const { return Chain<T, N>{corners, true}; }

  template <class K, std::enable_if_t<!std::is_same<K, T>::value>>
  operator Polygon<K, N>() const {
    std::vector<Point<K, N>> ps(corners.size());
    for (int i = 0; i < corners.size(); i++)
      ps[i] = corners[i];
    return Polygon<K, N>(std::move(ps), normal);
  }

  template <
      class = std::enable_if_t<std::is_floating_point<T>::value && N == 3>>
  T area() const {
    return Area(*this);
  }
};

template <class Archive, class T, int N>
void serialize(Archive &ar, Polygon<T, N> &p) {
  ar(p.corners, p.normal);
}

using Polygon2 = Polygon<double, 2>;
using Polygon2f = Polygon<float, 2>;
using Polygon3 = Polygon<double, 3>;
using Polygon3f = Polygon<float, 3>;

template <class T>
Polygon<T, 3> MakeTriangle(const Point<T, 3> &p1, const Point<T, 3> &p2,
                           const Point<T, 3> &p3) {
  return Polygon<T, 3>{{p1, p2, p3}, (p2 - p1).cross(p3 - p1)};
}

inline double PointTest(const Polygon2f &poly, const Point2f &p) {
  return cv::pointPolygonTest(poly.corners, p, true);
}
inline bool Contains(const Polygon2f &poly, const Point2f &p) {
  return cv::pointPolygonTest(poly.corners, p, false) >= 0;
}

// layered polygons
template <class T, int N> struct LayeredShape {
  std::vector<std::vector<Point<T, N>>> layers;
  Vec<T, N> normal;
  Polygon<T, N> layer(size_t i) const {
    return Polygon<T, N>(layers.at(i), normal);
  }
  size_t size() const { return layers.size(); }
  bool empty() const { return layers.empty(); }
};
template <class Archive, class T, int N>
void serialize(Archive &ar, LayeredShape<T, N> &p) {
  ar(p.layers, p.normal);
}
using LayeredShape3 = LayeredShape<double, 3>;
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
