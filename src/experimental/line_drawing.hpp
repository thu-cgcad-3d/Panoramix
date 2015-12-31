#pragma once

#include "../core/basic_types.hpp"
#include "../core/mesh.hpp"

namespace pano {
namespace experimental {

using namespace pano::core;

// LineDrawing
template <class PointT> struct LineDrawing {
  std::vector<PointT> corners;
  std::vector<std::pair<int, int>> line2corners;
  std::vector<std::vector<int>> face2corners;
  template <class ArchiverT> void serialize(ArchiverT &ar) {
    ar(corners, line2corners, face2corners);
  }
};

// LoadLineDrawing
LineDrawing<Point2> LoadLineDrawing(const std::string &filename);
LineDrawing<Point3> LoadLineDrawing(const std::string &filename,
                                    const std::string &gtfilename);

// ToMesh
template <class T> Mesh<T> ToMesh(const LineDrawing<T> &ld) {
  Mesh<T> mesh;
  using VHandle = typename Mesh<T>::VertHandle;
  mesh.internalVertices().reserve(ld.corners.size());
  mesh.internalHalfEdges().reserve(ld.line2corners.size() * 2);
  mesh.internalFaces().reserve(ld.face2corners.size());
  for (int i = 0; i < ld.corners.size(); i++) {
    mesh.addVertex(ld.corners[i]);
  }
  for (int i = 0; i < ld.line2corners.size(); i++) {
    mesh.addEdge(VHandle(ld.line2corners[i].first),
                 VHandle(ld.line2corners[i].second));
  }
  for (int i = 0; i < ld.face2corners.size(); i++) {
    std::vector<VHandle> vhs(ld.face2corners[i].size());
    for (int j = 0; j < vhs.size(); j++) {
      vhs[j] = VHandle(ld.face2corners[i][j]);
    }
    mesh.addFace(vhs, true);
  }
  return mesh;
}

// ToLineDrawing
template <class T, class H, class F>
LineDrawing<T> ToLineDrawing(const Mesh<T, H, F> &mesh) {
  LineDrawing<T> ld;
  using VHandle = typename Mesh<T, H, F>::VertHandle;
  std::vector<int> vh2corner(mesh.internalVertices().size(), -1);
  for (auto &v : mesh.vertices()) {
    ld.corners.push_back(v.data);
    vh2corner[v.topo.hd.id] = ld.corners.size() - 1;
  }
  std::vector<bool> hh2used(mesh.internalHalfEdges().size(), false);
  for (auto &h : mesh.halfedges()) {
    if (hh2used[h.topo.hd.id] || hh2used[mesh.topo(h.topo.hd).opposite.id]) {
      continue;
    }
    ld.line2corners.emplace_back(vh2corner[mesh.topo(h.topo.hd).from().id],
                                 vh2corner[mesh.topo(h.topo.hd).to().id]);
    hh2used[h.topo.hd.id] = true;
  }
  for (auto &f : mesh.faces()) {
    std::vector<int> corners(f.topo.halfedges.size());
    for (int j = 0; j < corners.size(); j++) {
      corners[j] = vh2corner[mesh.topo(f.topo.halfedges[j]).to().id];
    }
    ld.face2corners.push_back(std::move(corners));
  }
  return ld;
}

// Transform
template <class T, class FunT>
auto Transform(const LineDrawing<T> &in, const FunT &fun) {
  LineDrawing<typename FunctionTraits<FunT>::ResultType> out;
  out.corners.resize(in.corners.size());
  std::transform(in.corners.begin(), in.corners.end(), out.corners.begin(),
                 fun);
  out.line2corners = in.line2corners;
  out.face2corners = in.face2corners;
  return out;
}

template <class T, class FunT>
auto Transform(LineDrawing<T> &&in, const FunT &fun) {
  LineDrawing<typename FunctionTraits<FunT>::ResultType> out;
  out.corners.resize(in.corners.size());
  std::transform(in.corners.begin(), in.corners.end(), out.corners.begin(),
                 fun);
  out.line2corners = std::move(in.line2corners);
  out.face2corners = std::move(in.face2corners);
  return out;
}

// SearchFace
inline void SearchAndAddFaces(LineDrawing<Point2> &drawing) {
  auto mesh = ToMesh(drawing);
  SearchAndAddFaces(mesh);
  drawing = ToLineDrawing(mesh);
}
}

namespace core {
template <class PointT>
inline auto
BoundingBox(const pano::experimental::LineDrawing<PointT> &drawing) {
  return BoundingBoxOfContainer(drawing.corners);
}
}
}
