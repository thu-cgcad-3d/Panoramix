#include "cache.hpp"
#include "cameras.hpp"
#include "canvas.hpp"
#include "clock.hpp"
#include "line_drawing_widget.hpp"
#include "matlab_api.hpp"
#include "optimization.hpp"
#include "parallel.hpp"
#include "scene.hpp"
#include "singleton.hpp"

#include "line_drawing_reconstruction.hpp"
#include "line_drawing_tools.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

void RoutineTrainEnergyWeights() {
  int gui_level = 0;

  std::vector<LineDrawingInput> inputs = {
      LineDrawingInput::FromObjFile("hex", "cam1"),          //
      LineDrawingInput::FromObjFile("triangle", "cam1"),     //
      LineDrawingInput::FromObjFile("twotriangles", "cam1"), //
      LineDrawingInput::FromObjFile("stool", "cam1"),        //
      LineDrawingInput::FromObjFile("plane", "cam1"),        //
      LineDrawingInput::FromObjFile("towers", "cam1"),       // skewed a bit
      LineDrawingInput::FromObjFile("tower", "cam1"),        //
      LineDrawingInput::FromObjFile("towerx", "cam1"),       //
      LineDrawingInput::FromObjFile("car", "cam1"),          //
      LineDrawingInput::FromObjFile("gate", "cam1"),         // skewed a bit
      LineDrawingInput::FromObjFile("smallgate", "cam1"),    //
      LineDrawingInput::FromObjFile("castle", "cam1"),       //
      LineDrawingInput::FromObjFile("desk", "cam1"),         //
      // LineDrawingInput::FromImageAnnotation("house_sketch.jpg"), // <-
      // LineDrawingInput::FromImageAnnotation("steps.jpg")
      //
  };

  for (auto &input : inputs) {
    core::Println("#### input.model_name = [", input.model_name,
                  "], input.cam_name = [", input.cam_name, "]");
    const size_t npoints = input.projection.line_drawing.points.size();
    const size_t nedges = input.projection.line_drawing.topo.edges.size();
    const size_t nfaces = input.projection.line_drawing.topo.coplanar_points.size();

    ////////////////////////////////////////////////
    // STEP 1: compute additional data
    std::vector<std::set<int>> point2edges(npoints);
    std::map<std::pair<int, int>, int> points2edge;
    std::vector<std::set<int>> edge2faces(nedges);
    std::vector<std::set<int>> face2edges(nfaces);
    std::map<std::pair<int, int>, bool> faces_overlap;
    std::vector<Line2> edge2line(nedges);
    Box2 lines_box;
    {
      for (int edge = 0; edge < nedges; edge++) {
        auto &e = input.projection.line_drawing.topo.edges[edge];
        points2edge[e] = edge;
        point2edges[e.first].insert(edge);
        point2edges[e.second].insert(edge);
      }
      for (int face = 0; face < nfaces; face++) {
        auto &ps = input.projection.line_drawing.topo.coplanar_points[face];
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
        auto &cs = input.projection.line_drawing.topo.coplanar_points[face];
        auto &far_corners_for_edges = face2far_corners_for_edges[face];
        TriangulatePolygon(
            cs.begin(), cs.end(),
            [&input](int corner) -> const Point2 & {
              return input.projection.line_drawing.points.at(corner);
            },
            [&far_corners_for_edges, &input, &points2edge](int c1, int c2,
                                                           int c3) {
              if (Contains(points2edge, MakeOrderedPair(c1, c2))) {
                int edge = points2edge.at(MakeOrderedPair(c1, c2));
                far_corners_for_edges[edge] = c3;
              }
              if (Contains(points2edge, MakeOrderedPair(c2, c3))) {
                int edge = points2edge.at(MakeOrderedPair(c2, c3));
                far_corners_for_edges[edge] = c1;
              }
              if (Contains(points2edge, MakeOrderedPair(c3, c1))) {
                int edge = points2edge.at(MakeOrderedPair(c3, c1));
                far_corners_for_edges[edge] = c2;
              }
            });
      }
      // then we compute face overlapnesss
      for (int edge = 0; edge < nedges; edge++) {
        Line2 line(
            input.projection.line_drawing.points[input.projection.line_drawing.topo.edges[edge].first],
            input.projection.line_drawing.points[input.projection.line_drawing.topo.edges[edge].second]);
        auto &fs = edge2faces.at(edge);
        if (fs.size() < 2) {
          continue;
        }
        for (auto iti = fs.begin(); iti != fs.end(); ++iti) {
          int far_corner_i = face2far_corners_for_edges.at(*iti).at(edge);
          bool on_left_side_of_edge_i =
              IsOnLeftSide(input.projection.line_drawing.points.at(far_corner_i), line.first,
                           line.second);
          for (auto itj = std::next(iti); itj != fs.end(); ++itj) {
            int far_corner_j = face2far_corners_for_edges.at(*itj).at(edge);
            bool on_left_side_of_edge_j =
                IsOnLeftSide(input.projection.line_drawing.points.at(far_corner_j),
                             line.first, line.second);
            faces_overlap[MakeOrderedPair(*iti, *itj)] =
                on_left_side_of_edge_i == on_left_side_of_edge_j;
          }
        }
      }

      // construct 2d lines for edges
      for (int edge = 0; edge < nedges; edge++) {
        edge2line[edge].first =
            input.projection.line_drawing.points[input.projection.line_drawing.topo.edges[edge].first];
        edge2line[edge].second =
            input.projection.line_drawing.points[input.projection.line_drawing.topo.edges[edge].second];
      }
      lines_box = BoundingBoxOfContainer(edge2line);
    }

    ////////////////////////////////////////////////
    // STEP 2: define energy terms
    auto energy_terms = [&](const PerspectiveCamera &cam,
                            const std::vector<double> &point_depths) {
      assert(input.projection.line_drawing.points.size() ==
                 input.groundtruth->point_depths.size() &&
             input.projection.line_drawing.points.size() ==
                 point_depths.size());
      // compute 3d points from the given depths
      std::vector<Point3> points(input.projection.line_drawing.points.size());
      for (int i = 0; i < points.size(); i++) {
        points[i] = normalize(cam.direction(input.projection.line_drawing.points[i])) *
                    point_depths[i];
      }
      // compute face equations from the given depths
      std::vector<Vec3> face_equations(nfaces);
      for (int face = 0; face < nfaces; face++) {
        // todo
      }

      // angles between adjacent edges
      std::vector<double> angles_between_adjacent_edges;
      for (int face = 0; face < input.projection.line_drawing.topo.coplanar_points.size();
           face++) {
        auto &vs = input.projection.line_drawing.topo.coplanar_points[face];
        assert(vs.size() >= 3);
        for (int i = 0; i < vs.size(); i++) {
          int v1 = vs[i];
          int v2 = vs[(i + 1) % vs.size()];
          int v3 = vs[(i + 2) % vs.size()];
          if (!Contains(points2edge, MakeOrderedPair(v1, v2)) ||
              !Contains(points2edge, MakeOrderedPair(v2, v3))) {
            continue;
          }
          double angle = AngleBetweenDirected(points[v1] - points[v2],
                                              points[v3] - points[v2]);
          angles_between_adjacent_edges.push_back(angle);
        }
      }
      // [msd adjacent edges]
      double msd_adjacent_edges =
          MeanSquaredDeviationOfContainer(angles_between_adjacent_edges);
      // [ortho adjacent edges]
      double ortho_adjacent_edges = MeanSquaredDeviationOfContainer(
          angles_between_adjacent_edges, {M_PI_2, M_PI});

      // angles between adjacent faces
      std::vector<double> angles_between_adjacent_faces;
      // for(int edge = 0)

    };
  }
}