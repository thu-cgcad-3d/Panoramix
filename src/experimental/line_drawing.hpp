#pragma once

#include <unordered_map>
#include <unordered_set>

#include "../core/basic_types.hpp"
#include "../core/mesh.hpp"
#include "../core/mesh_util.hpp"

#include "../misc/eigen.hpp"

#include "optimization.hpp"

namespace pano {
namespace experimental {

using namespace pano::core;

// LineDrawing
enum LineType { SolidLine, DashedLine };
template <class PointT> struct LineDrawing {
  std::vector<PointT> corners;
  std::vector<std::pair<int, int>> line2corners;
  std::vector<LineType> line2type;
  std::vector<std::vector<int>> face2corners;
  template <class ArchiverT> void serialize(ArchiverT &ar) {
    ar(corners, line2corners, line2type, face2corners);
  }
};

// LoadLineDrawing
LineDrawing<Point2> LoadLineDrawing(const std::string &filename);
LineDrawing<Point3> LoadLineDrawing(const std::string &filename,
                                    const std::string &gtfilename);

// ToMesh
template <class T, class PointConvertT = std::identity<T>>
Mesh<std::decay_t<typename FunctionTraits<PointConvertT>::ResultType>>
ToMesh(const LineDrawing<T> &ld, PointConvertT cvtFun = PointConvertT());

// ToLineDrawing
template <class T, class H, class F, class VertConvertT = std::identity<T>>
LineDrawing<std::decay_t<typename FunctionTraits<VertConvertT>::ResultType>>
ToLineDrawing(const Mesh<T, H, F> &mesh, VertConvertT cvtFun = VertConvertT());

// Transform
template <class T, class FunT>
auto Transform(const LineDrawing<T> &in, const FunT &fun);
template <class T, class FunT>
auto Transform(LineDrawing<T> &&in, const FunT &fun);
}

namespace core {
template <class PointT>
inline auto
BoundingBox(const pano::experimental::LineDrawing<PointT> &drawing) {
  return BoundingBoxOfContainer(drawing.corners);
}
}
}

////////////////////////////////////////////////
//// implementations
////////////////////////////////////////////////
namespace pano {
namespace experimental {
// ToMesh
template <class T, class PointConvertT>
Mesh<std::decay_t<typename FunctionTraits<PointConvertT>::ResultType>>
ToMesh(const LineDrawing<T> &ld, PointConvertT cvtFun) {
  Mesh<std::decay_t<typename FunctionTraits<PointConvertT>::ResultType>> mesh;
  using VHandle = VertHandle;
  mesh.internalVertices().reserve(ld.corners.size());
  mesh.internalHalfEdges().reserve(ld.line2corners.size() * 2);
  mesh.internalFaces().reserve(ld.face2corners.size());
  for (int i = 0; i < ld.corners.size(); i++) {
    mesh.addVertex(cvtFun(ld.corners[i]));
  }
  for (int i = 0; i < ld.line2corners.size(); i++) {
    mesh.addEdge(VHandle(ld.line2corners[i].first),
                 VHandle(ld.line2corners[i].second));
  }
  if (ld.face2corners.empty()) {
    return mesh;
  }
  std::vector<int> faceInsertOrder = {0};
  std::set<int> notInsertedFace;
  for (int i = 1; i < ld.face2corners.size(); i++) {
    notInsertedFace.insert(i);
  }
  while (!notInsertedFace.empty()) {
  LabelInsertFace:
    for (int candFace : notInsertedFace) {
      for (int insertedFace : faceInsertOrder) {
        auto &candCorners = ld.face2corners[candFace];
        auto &insertedCorners = ld.face2corners[insertedFace];
        for (int i1 = 0; i1 < candCorners.size(); i1++) {
          int from1 = candCorners[i1];
          int to1 = candCorners[(i1 + 1) % candCorners.size()];
          for (int i2 = 0; i2 < insertedCorners.size(); i2++) {
            int from2 = insertedCorners[i2];
            int to2 = insertedCorners[(i2 + 1) % insertedCorners.size()];
            if (from1 == from2 && to1 == to2 || from1 == to2 && to1 == from2) {
              // has common edge!
              notInsertedFace.erase(candFace);
              faceInsertOrder.push_back(candFace);
              goto LabelInsertFace;
            }
          }
        }
      }
    }
  }

  for (int i : faceInsertOrder) {
    std::vector<VHandle> vhs(ld.face2corners[i].size());
    for (int j = 0; j < vhs.size(); j++) {
      vhs[j] = VHandle(ld.face2corners[i][j]);
    }
    mesh.addFace(vhs, true);
  }

  AssertEdgesAreStiched(mesh);
  return mesh;
}

// ToLineDrawing
template <class T, class H, class F, class VertConvertT>
LineDrawing<std::decay_t<typename FunctionTraits<VertConvertT>::ResultType>>
ToLineDrawing(const Mesh<T, H, F> &mesh, VertConvertT cvtFun) {
  LineDrawing<std::decay_t<typename FunctionTraits<VertConvertT>::ResultType>>
      ld;
  using VHandle = VertHandle;
  std::vector<int> vh2corner(mesh.internalVertices().size(), -1);
  for (auto &v : mesh.vertices()) {
    ld.corners.push_back(cvtFun(v.data));
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
  ld.line2type.resize(ld.line2corners.size(), SolidLine);
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
  out.line2type = in.line2type;
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
  out.line2type = std::move(in.line2type);
  out.face2corners = std::move(in.face2corners);
  return out;
}
}
}
