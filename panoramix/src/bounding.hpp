#include "math.hpp"

namespace pano {
namespace core {

// sphere
template <class T, int N> class Sphere {
public:
  Sphere() {}
  Sphere(const Point<T, N> &c, const T &r) : center(c), radius(r) {}

public:
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
template <class T, int N> class Box {
public:
  using Type = T;
  static const int Dimension = N;
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

public:
  Point<T, N> minCorner, maxCorner;
  bool isNull;
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
}
}