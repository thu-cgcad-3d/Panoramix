#pragma once

#include "../core/basic_types.hpp"

namespace pano {
namespace experimental {

using namespace pano::core;

// LineDrawing
template <class PointT> struct LineDrawing {
  std::vector<PointT> vertPositions;
  std::vector<std::pair<int, int>> line2verts;
  std::vector<std::vector<int>> face2verts;
  template <class ArchiverT> void serialize(ArchiverT &ar) {
    ar(vertPositions, line2verts, face2verts);
  }
};

// LoadLineDrawing
LineDrawing<Point2> LoadLineDrawing(const std::string &filename);
LineDrawing<Point3> LoadLineDrawing(const std::string &filename,
                                    const std::string &gtfilename);

// Transform
template <class T, class FunT>
auto Transform(const LineDrawing<T> &in, const FunT &fun) {
  LineDrawing<typename FunctionTraits<FunT>::ResultType> out;
  out.vertPositions.resize(in.vertPositions.size());
  std::transform(in.vertPositions.begin(), in.vertPositions.end(),
                 out.vertPositions.begin(), fun);
  out.line2verts = in.line2verts;
  out.face2verts = in.face2verts;
  return out;
}

template <class T, class FunT>
auto Transform(LineDrawing<T> &&in, const FunT &fun) {
  LineDrawing<typename FunctionTraits<FunT>::ResultType> out;
  out.vertPositions.resize(in.vertPositions.size());
  std::transform(in.vertPositions.begin(), in.vertPositions.end(),
                 out.vertPositions.begin(), fun);
  out.line2verts = std::move(in.line2verts);
  out.face2verts = std::move(in.face2verts);
  return out;
}

// SearchFace
void SearchFace(LineDrawing<Point2> &drawing);
}
}
