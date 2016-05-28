#pragma once

#include <unordered_map>
#include <unordered_set>

#include "basic_types.hpp"
#include "mesh.hpp"
#include "mesh_util.hpp"

#include "eigen.hpp"

#include "optimization.hpp"

namespace pano {
namespace experimental {

using namespace pano::core;

class LineDrawingTopo {
public:
  LineDrawingTopo() {}
  LineDrawingTopo(const std::vector<std::pair<int, int>> &e2cs,
                  const std::vector<std::vector<int>> &f2cs);

  template <class ArchiverT> void serialize(ArchiverT &ar) {
    ar(face2corners, edge2corners, face2edges, edge2faces, corner2edges,
       corners2edge, corner2faces);
  }

  size_t ncorners() const { return corner2edges.size(); }
  size_t nedges() const { return edge2corners.size(); }
  size_t nfaces() const { return face2corners.size(); }

public:
  std::vector<std::vector<int>> face2corners;
  std::vector<std::pair<int, int>> edge2corners;
  std::vector<std::vector<int>> face2edges;
  std::vector<std::vector<int>> edge2faces;
  std::vector<std::vector<int>> corner2edges;
  std::map<std::pair<int, int>, int> corners2edge;
  std::vector<std::vector<int>> corner2faces;
};

// LineDrawing
template <class CornerDataT, class EdgeDataT = Dummy, class FaceDataT = Dummy>
class LineDrawing {
public:
  LineDrawing() {}
  explicit LineDrawing(const std::vector<std::pair<int, int>> &e2cs,
                       const std::vector<std::vector<int>> &f2cs,
                       const std::vector<CornerDataT> &cs = {},
                       const std::vector<EdgeDataT> &es = {},
                       const std::vector<FaceDataT> &fs = {})
      : corners(cs), edges(es), faces(fs), topo(e2cs, f2cs) {
    if (corners.empty()) {
      corners.resize(topo.ncorners());
    }
    if (edges.empty()) {
      edges.resize(topo.nedges());
    }
    if (faces.empty()) {
      faces.resize(topo.nfaces());
    }
  }

  size_t ncorners() const { return corners.size(); }
  size_t nedges() const { return edges.size(); }
  size_t nfaces() const { return faces.size(); }

  template <class ArchiverT> void serialize(ArchiverT &ar) {
    ar(corners, edges, faces, topo);
  }

public:
  std::vector<CornerDataT> corners;
  std::vector<EdgeDataT> edges;
  std::vector<FaceDataT> faces;
  LineDrawingTopo topo;
};

// LoadLineDrawingFromObjFile
LineDrawing<Point3> LoadLineDrawingFromObjFile(const std::string & fname);

// ToMesh
template <class CT, class ET, class FT,
          class CornerConvertT = std::identity<CT>,
          class EdgeConvertT = std::identity<ET>,
          class FaceConvertT = std::identity<FT>>
Mesh<std::decay_t<typename FunctionTraits<CornerConvertT>::ResultType>,
     std::decay_t<typename FunctionTraits<EdgeConvertT>::ResultType>,
     std::decay_t<typename FunctionTraits<FaceConvertT>::ResultType>>
ToMesh(const LineDrawing<CT, ET, FT> &ld,
       CornerConvertT cornerCvtFun = CornerConvertT(),
       EdgeConvertT edgeCvtFun = EdgeConvertT(),
       FaceConvertT faceCvtFun = FaceConvertT());



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
template <class CT, class ET, class FT, class CornerConvertT,
          class EdgeConvertT, class FaceConvertT>
Mesh<std::decay_t<typename FunctionTraits<CornerConvertT>::ResultType>,
     std::decay_t<typename FunctionTraits<EdgeConvertT>::ResultType>,
     std::decay_t<typename FunctionTraits<FaceConvertT>::ResultType>>
ToMesh(const LineDrawing<CT, ET, FT> &ld, CornerConvertT cornerCvtFun,
       EdgeConvertT edgeCvtFun, FaceConvertT faceCvtFun) {
  Mesh<std::decay_t<typename FunctionTraits<PointConvertT>::ResultType>> mesh;
  using VHandle = VertHandle;
  mesh.internalVertices().reserve(ld.ncorners());
  mesh.internalHalfEdges().reserve(ld.nedges() * 2);
  mesh.internalFaces().reserve(ld.nfaces());
  for (int i = 0; i < ld.corners.size(); i++) {
    mesh.addVertex(cornerCvtFun(ld.corners[i]));
  }
  for (int i = 0; i < ld.topo.edge2corners.size(); i++) {
    mesh.addEdge(VHandle(ld.topo.edge2corners[i].first),
                 VHandle(ld.topo.edge2corners[i].second),
                 edgeCvtFun(ld.edges[i]));
  }
  if (ld.topo.face2corners.empty()) {
    return mesh;
  }
  std::vector<int> faceInsertOrder = {0};
  std::set<int> notInsertedFace;
  for (int i = 1; i < ld.topo.face2corners.size(); i++) {
    notInsertedFace.insert(i);
  }
  while (!notInsertedFace.empty()) {
  LabelInsertFace:
    for (int candFace : notInsertedFace) {
      for (int insertedFace : faceInsertOrder) {
        auto &candCorners = ld.topo.face2corners[candFace];
        auto &insertedCorners = ld.topo.face2corners[insertedFace];
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
    std::vector<VHandle> vhs(ld.topo.face2corners[i].size());
    for (int j = 0; j < vhs.size(); j++) {
      vhs[j] = VHandle(ld.topo.face2corners[i][j]);
    }
    mesh.addFace(vhs, true, faceCvtFun(ld.faces[i]));
  }

  AssertEdgesAreStiched(mesh);
  return mesh;
}
}
}
