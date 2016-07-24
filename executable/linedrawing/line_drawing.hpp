#pragma once

#include "cameras.hpp"
#include "geometry.hpp"
#include "image.hpp"
#include "discretization.hpp"

namespace pano {
namespace experimental {
using namespace ::pano::core;

// LineDrawingTopo
struct LineDrawingTopo {
  std::vector<std::pair<int, int>> edges;
  std::vector<std::vector<int>> coplanar_points;
  std::vector<std::vector<int>> colinear_points;
  template <class ArchiveT> void serialize(ArchiveT &ar) {
    ar(edges, coplanar_points, colinear_points);
  }
};

// LineDrawing
template <class PointT> struct LineDrawing {
  std::vector<PointT> points;
  LineDrawingTopo topo;
  template <class ArchiveT> void serialize(ArchiveT &ar) { ar(points, topo); }
};
using LineDrawing2 = LineDrawing<Point2>;
using LineDrawing3 = LineDrawing<Point3>;
LineDrawing3 LineDrawing3FromObjFile(const std::string &obj_file);
void SaveLineDrawing3ToObjFile(const std::string &obj_file,
                               const LineDrawing3 &line_drawing);

// AuxiliaryData
template <class EssentialDataT> struct AuxiliaryData {};

// AuxiliaryData<LineDrawing<PointT>>
template <class PointT> struct AuxiliaryData<LineDrawing<PointT>> {
  std::vector<std::set<int>> point2edges;
  std::map<std::pair<int, int>, int> points2edge;
  std::vector<std::set<int>> edge2faces;
  std::vector<std::set<int>> face2edges;
  std::vector<Line<PointT>> edge2line;
  template <class ArchiveT> void serialize(ArchiveT &ar) {
    ar(point2edges, points2edge, edge2faces, face2edges, edge2line, lines_box);
  }
};

// MakeAuxiliary
template <class PointT>
AuxiliaryData<LineDrawing<PointT>>
MakeAuxiliary(const LineDrawing<PointT> &line_drawing);

// ComputeFacesOverlap
std::map<std::pair<int, int>, bool>
ComputeFacesOverlap(const LineDrawing2 &line_drawing,
                    const AuxiliaryData<LineDrawing2> &aux);

// LineDrawingReconstruction
struct LineDrawingReconstruction {
  std::vector<double> point_depths;
  PerspectiveCamera camera;
  template <class Archive> void serialize(Archive &ar) {
    ar(point_depths, camera);
  }
};
// DecomposeProjectionAndDepths
std::pair<LineDrawing2, std::vector<double>>
DecomposeProjectionAndDepths(const LineDrawing3 &line_drawing,
                             const PerspectiveCamera &cam);
}

namespace core {
template <class PointT>
auto BoundingBox(const experimental::LineDrawing<PointT> &line_drawing) {
  return BoundingBoxOfContainer(line_drawing.points);
}
}

namespace gui {
template <class PointT>
void Discretize(TriMesh &mesh,
                const experimental::LineDrawing<PointT> &line_drawing,
                const DiscretizeOptions &o) {
  for (auto &e : line_drawing.topo.edges) {
    Discretize(mesh, Line<PointT>(line_drawing.points[e.first],
                                  line_drawing.points[e.second]),
               o);
  }
}
}
}

////////////////////////////////////////////////
//// implementations
////////////////////////////////////////////////
namespace pano {
namespace experimental {
template <class PointT>
AuxiliaryData<LineDrawing<PointT>>
MakeAuxiliary(const LineDrawing<PointT> &line_drawing) {
  const size_t npoints = line_drawing.points.size();
  const size_t nedges = line_drawing.topo.edges.size();
  const size_t nfaces = line_drawing.topo.coplanar_points.size();

  AuxiliaryData<LineDrawing<PointT>> aux;
  auto &point2edges = aux.point2edges;
  auto &points2edge = aux.points2edge;
  auto &edge2faces = aux.edge2faces;
  auto &face2edges = aux.face2edges;
  auto &edge2line = aux.edge2line;

  // compute additional data
  point2edges.resize(npoints);
  edge2faces.resize(nedges);
  face2edges.resize(nfaces);
  edge2line.resize(nedges);
  for (int edge = 0; edge < nedges; edge++) {
    auto &e = line_drawing.topo.edges[edge];
    points2edge[e] = edge;
    point2edges[e.first].insert(edge);
    point2edges[e.second].insert(edge);
  }
  for (int face = 0; face < nfaces; face++) {
    auto &ps = line_drawing.topo.coplanar_points[face];
    assert(ps.size() >= 3);
    for (int i = 0; i < ps.size(); i++) {
      int p1 = ps[i];
      int p2 = ps[(i + 1) % ps.size()];
      auto e = MakeOrderedPair(p1, p2);
      if (Contains(points2edge, e)) {
        int edge = points2edge.at(e);
        edge2faces[edge].insert(face);
        face2edges[face].insert(edge);
      }
    }
  }

  // construct 2d lines for edges
  for (int edge = 0; edge < nedges; edge++) {
    edge2line[edge].first =
        line_drawing.points[line_drawing.topo.edges[edge].first];
    edge2line[edge].second =
        line_drawing.points[line_drawing.topo.edges[edge].second];
  }

  return aux;
}
}
}
