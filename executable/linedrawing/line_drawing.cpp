#include "line_drawing.hpp"
#include "iterators.hpp"
#include "utility.hpp"

namespace pano {
namespace experimental {
// LineDrawing3FromObjFile
LineDrawing3 LineDrawing3FromObjFile(const std::string &obj_file) {
  LineDrawing3 ld;
  MakeMeshFromObjFile(
      [&ld](float x, float y, float z) { ld.points.emplace_back(x, y, z); },
      [&ld](auto &&face_corners) {
        ld.topo.coplanar_points.push_back(face_corners);
      },
      obj_file);
  std::set<std::pair<int, int>> edges;
  for (auto &f : ld.topo.coplanar_points) {
    assert(f.size() >= 3);
    for (int i = 0; i < f.size(); i++) {
      int c1 = f[i];
      int c2 = f[(i + 1) % f.size()];
      edges.insert(MakeOrderedPair(c1, c2));
    }
  }
  ld.topo.edges = std::vector<std::pair<int, int>>(edges.begin(), edges.end());
  return ld;
}

void SaveLineDrawing3ToObjFile(const std::string &obj_file,
                               const LineDrawing3 &line_drawing) {
  std::ofstream ofs(obj_file);
  if (ofs.is_open()) {
    for (int i = 0; i < line_drawing.points.size(); i++) {
      ofs << 'v';
      for (int k = 0; k < 3; k++) {
        ofs << ' ' << line_drawing.points[i][k];
      }
      ofs << '\n';
    }
    for (auto &f : line_drawing.topo.coplanar_points) {
      ofs << 'f';
      for (auto v : f) {
        ofs << ' ' << (v + 1);
      }
      ofs << '\n';
    }
  }
}

std::map<std::pair<int, int>, bool>
ComputeFacesOverlap(const LineDrawing2 &line_drawing,
                    const AuxiliaryData<LineDrawing2> &aux) {

  std::map<std::pair<int, int>, bool> faces_overlap;
  const size_t npoints = line_drawing.points.size();
  const size_t nedges = line_drawing.topo.edges.size();
  const size_t nfaces = line_drawing.topo.coplanar_points.size();

  // compute faces_overlap : whether two faces sharing a common edge
  // overlap?
  // first for each edge of a face, we select a corner of the face which
  // satisfies:
  // 1. the corner is not on the edge
  // 2. the triangle spanned by the corner and the edge is covered by the
  // face
  // polygon in 2d drawing
  std::vector<std::map<int, int>> face2far_corners_for_edges(nfaces);
  for (int face = 0; face < nfaces; face++) {
    auto &cs = line_drawing.topo.coplanar_points[face];
    auto &far_corners_for_edges = face2far_corners_for_edges[face];
    TriangulatePolygon(
        cs.begin(), cs.end(),
        [&line_drawing, &aux](int corner) -> const Point2 & {
          return line_drawing.points.at(corner);
        },
        [&far_corners_for_edges, &aux](int c1, int c2, int c3) {
          if (Contains(aux.points2edge, MakeOrderedPair(c1, c2))) {
            int edge = aux.points2edge.at(MakeOrderedPair(c1, c2));
            far_corners_for_edges[edge] = c3;
          }
          if (Contains(aux.points2edge, MakeOrderedPair(c2, c3))) {
            int edge = aux.points2edge.at(MakeOrderedPair(c2, c3));
            far_corners_for_edges[edge] = c1;
          }
          if (Contains(aux.points2edge, MakeOrderedPair(c3, c1))) {
            int edge = aux.points2edge.at(MakeOrderedPair(c3, c1));
            far_corners_for_edges[edge] = c2;
          }
        });
  }
  // then we compute face overlapnesss
  for (int edge = 0; edge < line_drawing.topo.edges.size(); edge++) {
    Line2 line(line_drawing.points[line_drawing.topo.edges[edge].first],
               line_drawing.points[line_drawing.topo.edges[edge].second]);
    auto &fs = aux.edge2faces.at(edge);
    if (fs.size() < 2) {
      continue;
    }
    for (auto iti = fs.begin(); iti != fs.end(); ++iti) {
      int far_corner_i = face2far_corners_for_edges.at(*iti).at(edge);
      bool on_left_side_of_edge_i = IsOnLeftSide(
          line_drawing.points.at(far_corner_i), line.first, line.second);
      for (auto itj = std::next(iti); itj != fs.end(); ++itj) {
        int far_corner_j = face2far_corners_for_edges.at(*itj).at(edge);
        bool on_left_side_of_edge_j = IsOnLeftSide(
            line_drawing.points.at(far_corner_j), line.first, line.second);
        faces_overlap[MakeOrderedPair(*iti, *itj)] =
            on_left_side_of_edge_i == on_left_side_of_edge_j;
      }
    }
  }

  return faces_overlap;
}

std::pair<LineDrawing2, std::vector<double>>
DecomposeProjectionAndDepths(const LineDrawing3 &line_drawing,
                             const PerspectiveCamera &cam) {
  std::vector<double> point_depths =
      MakeRange(line_drawing.points)
          .transform([&cam](const Point3 &p3) -> double {
            return Distance(p3, cam.eye());
          })
          .evalAsStdVector();
  double mean_depth =
      std::accumulate(point_depths.begin(), point_depths.end(), 0.0) /
      point_depths.size();
  for (auto &d : point_depths) {
    d /= mean_depth;
  }

  LineDrawing2 projection;
  projection.topo = line_drawing.topo;
  projection.points = MakeRange(line_drawing.points)
                          .transform([&cam](const Point3 &p3) -> Point2 {
                            return cam.toScreen(p3);
                          })
                          .evalAsStdVector();
  return std::make_pair(projection, point_depths);
}


}
}