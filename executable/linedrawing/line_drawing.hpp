#pragma once

#include <unordered_map>
#include <unordered_set>

#include "basic_types.hpp"
#include "mesh.hpp"
#include "mesh_util.hpp"

#include "eigen.hpp"

#include "optimization.hpp"

#include "discretization.hpp"

namespace pano {
namespace experimental {
using namespace pano::core;

template <class VertDataT, class EdgeDataT = Dummy, class FaceDataT = Dummy>
class LineDrawing {
public:
  using VertData = VertDataT;
  using EdgeData = EdgeDataT;
  using FaceData = FaceDataT;

  struct VertTopo;
  struct EdgeTopo;
  struct FaceTopo;
  using VertHandle = Handle<VertTopo>;
  using EdgeHandle = Handle<EdgeTopo>;
  using FaceHandle = Handle<FaceTopo>;

  struct VertTopo {
    VertHandle hd;
    std::unordered_set<EdgeHandle> edges;
    std::unordered_set<FaceHandle> faces;
    template <class Archive> void serialize(Archive &ar) {
      ar(hd, edges, faces);
    }
  };
  struct EdgeTopo {
    EdgeHandle hd;
    VertHandle vertices[2];
    VertHandle &from() { return vertices[0]; }
    VertHandle &to() { return vertices[1]; }
    const VertHandle &from() const { return vertices[0]; }
    const VertHandle &to() const { return vertices[1]; }
    std::unordered_set<FaceHandle> faces;
    template <class Archive> void serialize(Archive &ar) {
      ar(hd, vertices, faces);
    }
  };
  struct FaceTopo {
    FaceHandle hd;
    std::vector<VertHandle> vertices;
    template <class Archive> void serialize(Archive &ar) { ar(hd, vertices); }
  };

public:
  decltype(auto) internalVertices() { return _verts; }
  decltype(auto) internalEdges() { return _edges; }
  decltype(auto) internalFaces() { return _faces; }
  decltype(auto) internalVertices() const { return _verts; }
  decltype(auto) internalEdges() const { return _edges; }
  decltype(auto) internalFaces() const { return _faces; }

  ArrayView<const VertDataT> internalVertsData() const {
    return ArrayView<const VertDataT>(&(_verts.front().data), _verts.size(),
                                      sizeof(Triplet<VertTopo, VertData>));
  }
  ArrayView<const EdgeDataT> internalEdgesData() const {
    return ArrayView<const EdgeDataT>(&(_edges.front().data), _edges.size(),
                                      sizeof(Triplet<EdgeTopo, EdgeData>));
  }
  ArrayView<const FaceDataT> internalFacesData() const {
    return ArrayView<const FaceDataT>(&(_faces.front().data), _faces.size(),
                                      sizeof(Triplet<FaceTopo, FaceData>));
  }

  auto vertices() {
    return MakeConditionalRange(_verts,
                                TripletExistsPred<VertTopo, VertDataT>());
  }
  auto edges() {
    return MakeConditionalRange(_edges,
                                TripletExistsPred<EdgeTopo, EdgeData>());
  }
  auto faces() {
    return MakeConditionalRange(_faces,
                                TripletExistsPred<FaceTopo, FaceDataT>());
  }
  auto vertices() const {
    return MakeConditionalRange(_verts,
                                TripletExistsPred<VertTopo, VertDataT>());
  }
  auto edges() const {
    return MakeConditionalRange(_edges,
                                TripletExistsPred<EdgeTopo, EdgeData>());
  }
  auto faces() const {
    return MakeConditionalRange(_faces,
                                TripletExistsPred<FaceTopo, FaceDataT>());
  }

  VertTopo &topo(VertHandle v) { return _verts[v.id].topo; }
  EdgeTopo &topo(EdgeHandle h) { return _edges[h.id].topo; }
  FaceTopo &topo(FaceHandle f) { return _faces[f.id].topo; }
  const VertTopo &topo(VertHandle v) const { return _verts[v.id].topo; }
  const EdgeTopo &topo(EdgeHandle h) const { return _edges[h.id].topo; }
  const FaceTopo &topo(FaceHandle f) const { return _faces[f.id].topo; }

  VertDataT &data(VertHandle v) { return _verts[v.id].data; }
  EdgeDataT &data(EdgeHandle h) { return _edges[h.id].data; }
  FaceDataT &data(FaceHandle f) { return _faces[f.id].data; }
  const VertDataT &data(VertHandle v) const { return _verts[v.id].data; }
  const EdgeDataT &data(EdgeHandle h) const { return _edges[h.id].data; }
  const FaceDataT &data(FaceHandle f) const { return _faces[f.id].data; }

  template <class VDT = VertDataT> VertHandle addVertex(VDT &&vd = VDT()) {
    VertTopo topo;
    topo.hd.id = _verts.size();
    _verts.emplace_back(std::move(topo), std::forward<VDT>(vd), true);
    return _verts.back().topo.hd;
  }

  template <class EDT = EdgeDataT>
  EdgeHandle addEdge(VertHandle from, VertHandle to, EDT &&ed = EDT(),
                     bool mergeDuplicateEdge = true) {
    if (from == to) {
      return EdgeHandle();
    }
    EdgeHandle eh;
    if (mergeDuplicateEdge) {
      eh = findEdge(from, to);
    }
    if (eh.valid()) {
      return eh;
    } else {
      eh = EdgeHandle(_edges.size());
      Triplet<EdgeTopo, EdgeDataT> et;
      et.topo.hd = eh;
      et.topo.from() = from;
      et.topo.to() = to;
      et.exists = true;
      et.data = std::forward<EDT>(ed);
      _edges.push_back(std::move(et));
      _verts[from.id].topo.edges.insert(eh);
      _verts[to.id].topo.edges.insert(eh);
      return eh;
    }
  }

  template <class FDT = FaceDataT>
  FaceHandle addFace(const std::vector<VertHandle> &verts, FDT &&fd = FDT(),
                     bool addBoundaryEdges = true) {
    assert(verts.size() >= 3);
    Triplet<FaceTopo, FaceDataT> ft;
    ft.topo.hd.id = _faces.size();
    ft.topo.vertices = verts;
    ft.exists = true;
    ft.data = std::forward<FDT>(fd);
    _faces.push_back(std::move(ft));

    for (int i = 0; i < verts.size(); i++) {
      VertHandle vh1 = verts[i];
      VertHandle vh2 = verts[(i + 1) % verts.size()];
      _verts[vh1.id].topo.faces.insert(_faces.back().topo.hd);
      if (addBoundaryEdges) {
        EdgeHandle eh = findEdge(vh1, vh2);
        if (eh.invalid()) {
          eh = addEdge(vh1, vh2);
        }
        _edges[eh.id].topo.faces.insert(_faces.back().topo.hd);
      }
    }
    return _faces.back().topo.hd;
  }

  EdgeHandle findEdge(VertHandle from, VertHandle to) const {
    for (EdgeHandle eh : _verts[from.id].topo.edges) {
      assert(_edges[eh.id].topo.vertices[0] == from ||
             _edges[eh.id].topo.vertices[1] == from);
      if (_edges[eh.id].topo.vertices[0] == to ||
          _edges[eh.id].topo.vertices[1] == to) {
        return eh;
      }
    }
    return EdgeHandle();
  }

  size_t degree(VertHandle v) const { return _verts[v.id].topo.edges.size(); }
  size_t degree(FaceHandle f) const {
    return _faces[f.id].topo.vertices.size();
  }

  void clear() {
    _verts.clear();
    _halfs.clear();
    _faces.clear();
  }

private:
  TripletArray<VertTopo, VertData> _verts;
  TripletArray<EdgeTopo, EdgeData> _edges;
  TripletArray<FaceTopo, FaceData> _faces;
};



class LineDrawingTopo {
public:
  LineDrawingTopo() {}
  LineDrawingTopo(const std::vector<std::pair<int, int>> &e2cs,
                  const std::vector<std::vector<int>> &f2cs);

  template <class ArchiverT> void serialize(ArchiverT &ar) {
    ar(face2corners, edge2corners, face2edges, edge2faces, corner2edges,
       corners2edge, corner2faces, edge_face2same_direction,
       adjecent_faces2same_direction);
  }

  size_t ncorners() const { return corner2edges.size(); }
  size_t nedges() const { return edge2corners.size(); }
  size_t nfaces() const { return face2corners.size(); }

  bool maybeManifold() const;

public:
  std::vector<std::vector<int>> face2corners;
  std::vector<std::pair<int, int>> edge2corners;
  std::vector<std::vector<int>> face2edges;
  std::vector<std::vector<int>> edge2faces;
  std::vector<std::vector<int>> corner2edges;
  std::map<std::pair<int, int>, int> corners2edge;
  std::vector<std::vector<int>> corner2faces;
  std::map<std::pair<int, int>, bool> edge_face2same_direction;
  std::map<std::pair<int, int>, bool> adjecent_faces2same_direction;
};

//// LineDrawing
//template <class CornerDataT, class EdgeDataT = Dummy, class FaceDataT = Dummy>
//class LineDrawingOld {
//public:
//  LineDrawingOld() {}
//  explicit LineDrawingOld(const std::vector<std::pair<int, int>> &e2cs,
//                          const std::vector<std::vector<int>> &f2cs,
//                          const std::vector<CornerDataT> &cs = {},
//                          const std::vector<EdgeDataT> &es = {},
//                          const std::vector<FaceDataT> &fs = {})
//      : corners(cs), edges(es), faces(fs), topo(e2cs, f2cs) {
//    if (corners.empty()) {
//      corners.resize(topo.ncorners());
//    }
//    if (edges.empty()) {
//      edges.resize(topo.nedges());
//    }
//    if (faces.empty()) {
//      faces.resize(topo.nfaces());
//    }
//  }
//
//  size_t ncorners() const { return corners.size(); }
//  size_t nedges() const { return edges.size(); }
//  size_t nfaces() const { return faces.size(); }
//
//  template <class ArchiverT> void serialize(ArchiverT &ar) {
//    ar(corners, edges, faces, topo);
//  }
//
//public:
//  std::vector<CornerDataT> corners;
//  std::vector<EdgeDataT> edges;
//  std::vector<FaceDataT> faces;
//  LineDrawingTopo topo;
//};

// LoadLineDrawingFromObjFile
//LineDrawingOld<Point3> LoadLineDrawingOldFromObjFile(const std::string &fname);

//// ToMesh
// template <class CT, class ET, class FT,
//          class CornerConvertT = std::identity<CT>,
//          class EdgeConvertT = std::identity<ET>,
//          class FaceConvertT = std::identity<FT>>
// Mesh<std::decay_t<typename FunctionTraits<CornerConvertT>::ResultType>,
//     std::decay_t<typename FunctionTraits<EdgeConvertT>::ResultType>,
//     std::decay_t<typename FunctionTraits<FaceConvertT>::ResultType>>
// ToMesh(const LineDrawing<CT, ET, FT> &ld,
//       CornerConvertT cornerCvtFun = CornerConvertT(),
//       EdgeConvertT edgeCvtFun = EdgeConvertT(),
//       FaceConvertT faceCvtFun = FaceConvertT());
}

// namespace core {
// template <class PointT>
// inline auto
// BoundingBox(const pano::experimental::LineDrawing<PointT> &drawing) {
//  return BoundingBoxOfContainer(drawing.corners);
//}
//}
//
namespace gui {
//template <class T, class E, class F>
//void Discretize(TriMesh &mesh,
//                const experimental::LineDrawingOld<Point<T, 3>, E, F> &ld,
//                const DiscretizeOptions &o);
}
}

//////////////////////////////////////////////////
////// implementations
//////////////////////////////////////////////////
namespace pano {
// namespace experimental {
//
//// ToMesh
// template <class CT, class ET, class FT, class CornerConvertT,
//          class EdgeConvertT, class FaceConvertT>
// Mesh<std::decay_t<typename FunctionTraits<CornerConvertT>::ResultType>,
//     std::decay_t<typename FunctionTraits<EdgeConvertT>::ResultType>,
//     std::decay_t<typename FunctionTraits<FaceConvertT>::ResultType>>
// ToMesh(const LineDrawing<CT, ET, FT> &ld, CornerConvertT cornerCvtFun,
//       EdgeConvertT edgeCvtFun, FaceConvertT faceCvtFun) {
//  Mesh<std::decay_t<typename FunctionTraits<PointConvertT>::ResultType>> mesh;
//  using VHandle = VertHandle;
//  mesh.internalVertices().reserve(ld.ncorners());
//  mesh.internalHalfEdges().reserve(ld.nedges() * 2);
//  mesh.internalFaces().reserve(ld.nfaces());
//  for (int i = 0; i < ld.corners.size(); i++) {
//    mesh.addVertex(cornerCvtFun(ld.corners[i]));
//  }
//  for (int i = 0; i < ld.topo.edge2corners.size(); i++) {
//    mesh.addEdge(VHandle(ld.topo.edge2corners[i].first),
//                 VHandle(ld.topo.edge2corners[i].second),
//                 edgeCvtFun(ld.edges[i]));
//  }
//  if (ld.topo.face2corners.empty()) {
//    return mesh;
//  }
//  std::vector<int> faceInsertOrder = {0};
//  std::set<int> notInsertedFace;
//  for (int i = 1; i < ld.topo.face2corners.size(); i++) {
//    notInsertedFace.insert(i);
//  }
//  while (!notInsertedFace.empty()) {
//  LabelInsertFace:
//    for (int candFace : notInsertedFace) {
//      for (int insertedFace : faceInsertOrder) {
//        auto &candCorners = ld.topo.face2corners[candFace];
//        auto &insertedCorners = ld.topo.face2corners[insertedFace];
//        for (int i1 = 0; i1 < candCorners.size(); i1++) {
//          int from1 = candCorners[i1];
//          int to1 = candCorners[(i1 + 1) % candCorners.size()];
//          for (int i2 = 0; i2 < insertedCorners.size(); i2++) {
//            int from2 = insertedCorners[i2];
//            int to2 = insertedCorners[(i2 + 1) % insertedCorners.size()];
//            if (from1 == from2 && to1 == to2 || from1 == to2 && to1 == from2)
//            {
//              // has common edge!
//              notInsertedFace.erase(candFace);
//              faceInsertOrder.push_back(candFace);
//              goto LabelInsertFace;
//            }
//          }
//        }
//      }
//    }
//  }
//
//  for (int i : faceInsertOrder) {
//    std::vector<VHandle> vhs(ld.topo.face2corners[i].size());
//    for (int j = 0; j < vhs.size(); j++) {
//      vhs[j] = VHandle(ld.topo.face2corners[i][j]);
//    }
//    mesh.addFace(vhs, true, faceCvtFun(ld.faces[i]));
//  }
//
//  AssertEdgesAreStiched(mesh);
//  return mesh;
//}
//}
//
namespace gui {
//template <class T, class E, class F>
//void Discretize(TriMesh &mesh,
//                const experimental::LineDrawingOld<Point<T, 3>, E, F> &ld,
//                const DiscretizeOptions &o) {
//  std::vector<TriMesh::VertHandle> vhandles;
//  for (const Point<T, 3> &c : ld.corners) {
//    TriMesh::Vertex v;
//    v.position = cat(c, (T)1.0);
//    v.color = o.color();
//    vhandles.push_back(mesh.addVertex(v, true, o.entity()));
//  }
//  for (int edge = 0; edge < ld.nedges(); edge++) {
//    mesh.addLine(vhandles[ld.topo.edge2corners[edge].first],
//                 vhandles[ld.topo.edge2corners[edge].second], o.entity());
//  }
//  for (int face = 0; face < ld.nfaces(); face++) {
//    std::vector<TriMesh::VertHandle> vhs;
//    vhs.reserve(ld.topo.face2corners[face].size());
//    for (int c : ld.topo.face2corners[face]) {
//      vhs.push_back(vhandles[c]);
//    }
//    mesh.addPolygon(vhs, o.entity());
//  }
//}
}
}
