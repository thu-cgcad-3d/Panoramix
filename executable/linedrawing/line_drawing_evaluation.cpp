#include "line_drawing_evaluation.hpp"
#include "line_drawing_tools.hpp"
#include "utility.hpp"

namespace pano {
namespace experimental {

std::function<std::vector<double>(const std::vector<double> &point_depths)>
MakeLineDrawingFeatureExtractor(
    const LineDrawing2 &line_drawing, const AuxiliaryData<LineDrawing2> &aux,
    const std::map<std::pair<int, int>, bool> &faces_overlap,
    const std::vector<std::set<int>> &face_sets, const PerspectiveCamera &cam) {

  const size_t npoints = line_drawing.points.size();
  const size_t nedges = line_drawing.topo.edges.size();
  const size_t nfaces = line_drawing.topo.coplanar_points.size();

  static const double clamp_thres_for_angle_ortho = DegreesToRadians(20);

  return [&line_drawing, &aux, &faces_overlap, &face_sets, &cam, npoints,
          nedges, nfaces](
             const std::vector<double> &point_depths) -> std::vector<double> {

    // compute 3d points from the given depths
    std::vector<Point3> points(line_drawing.points.size());
    for (int i = 0; i < points.size(); i++) {
      points[i] =
          normalize(cam.direction(line_drawing.points[i])) * point_depths[i] +
          cam.eye();
    }
    // compute face equations from the given depths
    std::vector<Plane3> face_planes(nfaces);
    for (int face = 0; face < nfaces; face++) {
      auto &vs = line_drawing.topo.coplanar_points[face];
      assert(vs.size() >= 3);
      bool found = false;
      for (int i = 0; i < vs.size() && !found; i++) {
        for (int j = i + 1; j < vs.size() && !found; j++) {
          for (int k = j + 1; k < vs.size() && !found; k++) {
            Plane3 plane = Plane3From3Points(points[i], points[j], points[k]);
            if (plane.normal != Origin()) {
              if (plane.normal.dot(plane.anchor - cam.eye()) <
                  0) { // make sure that all plane normals are toward the
                       // forward direction
                plane.normal = -plane.normal;
              }
              face_planes[face] = plane;
              found = true;
              break;
            }
          }
        }
      }
    }

    // angles between adjacent edges
    std::vector<double> angles_between_adjacent_edges;
    for (int face = 0; face < line_drawing.topo.coplanar_points.size();
         face++) {
      auto &vs = line_drawing.topo.coplanar_points[face];
      assert(vs.size() >= 3);
      for (int i = 0; i < vs.size(); i++) {
        int v1 = vs[i];
        int v2 = vs[(i + 1) % vs.size()];
        int v3 = vs[(i + 2) % vs.size()];
        if (!Contains(aux.points2edge, MakeOrderedPair(v1, v2)) ||
            !Contains(aux.points2edge, MakeOrderedPair(v2, v3))) {
          continue;
        }
        double angle = AngleBetweenDirected(points[v1] - points[v2],
                                            points[v3] - points[v2]);
        angles_between_adjacent_edges.push_back(angle);
      }
    }
    // [msda adjacent edges]
    double msda_adjacent_edges =
        MeanSquaredDeviationOfContainer(angles_between_adjacent_edges);
    // [ortho adjacent edges]
    double ortho_adjacent_edges = MeanSquaredDeviationOfContainer(
        angles_between_adjacent_edges, {M_PI_2, M_PI},
        clamp_thres_for_angle_ortho);

    // angles between adjacent faces
    std::vector<double> angles_between_adjacent_faces;
    for (auto &fs : aux.edge2faces) {
      if (fs.size() < 2) {
        continue;
      }
      for (auto it1 = fs.begin(); it1 != fs.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != fs.end(); ++it2) {
          double angle =
              AngleBetweenDirected(face_planes[*it1].normal,
                                   faces_overlap.at(MakeOrderedPair(*it1, *it2))
                                       ? face_planes[*it2].normal
                                       : -face_planes[*it2].normal);
          angles_between_adjacent_faces.push_back(angle);
        }
      }
    }
    // [msda adjacent faces]
    double msda_adjacent_faces =
        MeanSquaredDeviationOfContainer(angles_between_adjacent_faces);
    // [ortho adjacent faces]
    double ortho_adjacent_faces = MeanSquaredDeviationOfContainer(
        angles_between_adjacent_faces, {M_PI_2, M_PI, M_PI_2 * 3},
        clamp_thres_for_angle_ortho);


    // for each face set
    double mean_msda_adjacent_edges_face_sets = 0.0;
    double mean_ortho_adjacent_edges_face_sets = 0.0;
    double mean_msda_adjacent_faces_face_sets = 0.0;
    double mean_ortho_adjacent_faces_face_sets = 0.0;
    for (int ii = 0; ii < face_sets.size(); ii++) {
      // angles_between_adjacent_edges_this_face_set
      std::vector<double> angles_between_adjacent_edges_this_face_set;
      for (int face : face_sets[ii]) {
        auto &vs = line_drawing.topo.coplanar_points[face];
        assert(vs.size() >= 3);
        for (int i = 0; i < vs.size(); i++) {
          int v1 = vs[i];
          int v2 = vs[(i + 1) % vs.size()];
          int v3 = vs[(i + 2) % vs.size()];
          if (!Contains(aux.points2edge, MakeOrderedPair(v1, v2)) ||
              !Contains(aux.points2edge, MakeOrderedPair(v2, v3))) {
            continue;
          }
          double angle = AngleBetweenDirected(points[v1] - points[v2],
                                              points[v3] - points[v2]);
          angles_between_adjacent_edges_this_face_set.push_back(angle);
        }
      }
      mean_msda_adjacent_edges_face_sets += MeanSquaredDeviationOfContainer(
          angles_between_adjacent_edges_this_face_set);
      mean_ortho_adjacent_edges_face_sets += MeanSquaredDeviationOfContainer(
          angles_between_adjacent_edges_this_face_set, {M_PI_2, M_PI},
          clamp_thres_for_angle_ortho);

      // angles_between_adjacent_faces_this_face_set
      std::vector<double> angles_between_adjacent_faces_this_face_set;
      for (auto &fs : aux.edge2faces) {
        if (fs.size() < 2) {
          continue;
        }
        for (auto it1 = fs.begin(); it1 != fs.end(); ++it1) {
          for (auto it2 = std::next(it1); it2 != fs.end(); ++it2) {
            if (!Contains(face_sets[ii], *it1) ||
                !Contains(face_sets[ii], *it2)) {
              continue;
            }
            double angle = AngleBetweenDirected(
                face_planes[*it1].normal,
                faces_overlap.at(MakeOrderedPair(*it1, *it2))
                    ? face_planes[*it2].normal
                    : -face_planes[*it2].normal);
            angles_between_adjacent_faces_this_face_set.push_back(angle);
          }
        }
      }
      mean_msda_adjacent_faces_face_sets += MeanSquaredDeviationOfContainer(
          angles_between_adjacent_faces_this_face_set);
      mean_ortho_adjacent_faces_face_sets += MeanSquaredDeviationOfContainer(
          angles_between_adjacent_faces_this_face_set,
          {M_PI_2, M_PI, M_PI_2 * 3}, clamp_thres_for_angle_ortho);
    }
    mean_msda_adjacent_edges_face_sets /= face_sets.size();
    mean_ortho_adjacent_edges_face_sets /= face_sets.size();

    mean_msda_adjacent_faces_face_sets /= face_sets.size();
    mean_ortho_adjacent_faces_face_sets /= face_sets.size();

    // edge skew ratio
    double max_edge_skew_ratio_score = 0.0; // max ratio / median ratio
    {
      std::vector<double> edge_proj_ratios(nedges);
      for (int edge = 0; edge < nedges; edge++) {
        double length2d = aux.edge2line[edge].length();
        auto &edge_corners = line_drawing.topo.edges[edge];
        double length3d =
            Distance(points[edge_corners.first], points[edge_corners.second]);
        edge_proj_ratios[edge] = length3d / length2d;
      }
      std::sort(edge_proj_ratios.begin(), edge_proj_ratios.end());
      max_edge_skew_ratio_score = edge_proj_ratios.back() /
                                  edge_proj_ratios[edge_proj_ratios.size() / 2];
    }

    return {msda_adjacent_edges,
            ortho_adjacent_edges,
            msda_adjacent_faces,
            ortho_adjacent_faces,
            mean_msda_adjacent_edges_face_sets,
            mean_ortho_adjacent_edges_face_sets,
            mean_msda_adjacent_faces_face_sets,
            mean_ortho_adjacent_faces_face_sets,
            max_edge_skew_ratio_score};
  };
}
}
}