#include "cache.hpp"
#include "cameras.hpp"
#include "canvas.hpp"
#include "line_drawing.hpp"
#include "scene.hpp"
#include "singleton.hpp"
#include "optimization.hpp"
#include "parallel.hpp"
#include "clock.hpp"

#include "tools.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

struct Configuration {
  std::string model_name, cam_name;
  LineDrawingTopo topo;
  std::vector<Point2> corners2d;
  std::map<std::pair<int, int>, bool> faces_overlap;
  std::vector<Point3> cornersGT;
  PerspectiveCamera cameraGT;

  template <class... Ts>
  bool loadCache(const std::string &what, Ts &... ts) const {
    return misc::LoadCache(model_name + "_" + cam_name, what, ts...);
  }
  template <class... Ts>
  bool saveCache(const std::string &what, Ts &... ts) const {
    return misc::SaveCache(model_name + "_" + cam_name, what, ts...);
  }
};

// ParseConfig
Configuration ParseConfig(const std::string &modelName,
                          const std::string &camName) {
  std::string folder = "F:\\LineDrawings\\dataset\\" + modelName + "\\";

  auto line_drawing_gt =
      LoadLineDrawingFromObjFile(folder + modelName + "_w_intf.obj");
  assert(line_drawing_gt.ncorners() > 0);
  PerspectiveCamera cam;
  std::string camPath = folder + modelName + ".obj." + camName + ".cereal";
  bool succ = LoadFromDisk(camPath, cam);
  if (!succ) {
    gui::SceneBuilder sb;
    sb.add(line_drawing_gt);
    cam = sb.show(true, true, gui::RenderOptions()
                                  .renderMode(gui::Lines)
                                  .fixUpDirectionInCameraMove(false))
              .camera();
    SaveToDisk(camPath, cam);
  }

  core::Println("gt focal = ", cam.focal(), " gt pp = ", cam.principlePoint());
  std::vector<Point2> corners2d(line_drawing_gt.corners.size());
  for (int i = 0; i < corners2d.size(); i++) {
    corners2d[i] = cam.toScreen(line_drawing_gt.corners[i]);
  }
  // compute face overlapness
  // first for each edge of a face, we select a corner of the face which satisfies:
  // 1. the corner is not on the edge
  // 2. the triangle spanned by the corner and the edge is covered by the face
  // polygon in 2d drawing
  auto &topo = line_drawing_gt.topo;
  std::vector<std::map<int, int>> face2far_corners_for_edges(topo.nfaces());
  for (int face = 0; face < face2far_corners_for_edges.size(); face++) {
    auto &cs = topo.face2corners[face];
    auto &far_corners_for_edges = face2far_corners_for_edges[face];
    TriangulatePolygon(
        cs.begin(), cs.end(),
        [&corners2d](int corner) { return corners2d.at(corner); },
        [&far_corners_for_edges, &topo](int c1, int c2, int c3) {
          if (Contains(topo.corners2edge, MakeOrderedPair(c1, c2))) {
            int edge = topo.corners2edge.at(MakeOrderedPair(c1, c2));
            far_corners_for_edges[edge] = c3;
          }
          if (Contains(topo.corners2edge, MakeOrderedPair(c2, c3))) {
            int edge = topo.corners2edge.at(MakeOrderedPair(c2, c3));
            far_corners_for_edges[edge] = c1;
          }
          if (Contains(topo.corners2edge, MakeOrderedPair(c3, c1))) {
            int edge = topo.corners2edge.at(MakeOrderedPair(c3, c1));
            far_corners_for_edges[edge] = c2;
          }
        });
  }
  // the we compute face overlapnesss
  std::map<std::pair<int, int>, bool> faces_overlap;
  for (int edge = 0; edge < topo.nedges(); edge++) {
    Line2 line(corners2d[topo.edge2corners[edge].first],
               corners2d[topo.edge2corners[edge].second]);
    auto &fs = topo.edge2faces.at(edge);
    if (fs.size() < 2) {
      continue;
    }
    for (int i = 0; i < fs.size(); i++) {
      int far_corner_i = face2far_corners_for_edges.at(fs[i]).at(edge);
      bool on_left_side_of_edge_i =
          IsOnLeftSide(corners2d.at(far_corner_i), line.first, line.second);
      for (int j = i + 1; j < fs.size(); j++) {
        int far_corner_j = face2far_corners_for_edges.at(fs[j]).at(edge);
        bool on_left_side_of_edge_j =
            IsOnLeftSide(corners2d.at(far_corner_j), line.first, line.second);
        faces_overlap[MakeOrderedPair(fs[i], fs[j])] =
            on_left_side_of_edge_i == on_left_side_of_edge_j;
      }
    }
  }

  return Configuration{modelName,
                       camName,
                       std::move(line_drawing_gt.topo),
                       std::move(corners2d),
                       std::move(faces_overlap),
                       std::move(line_drawing_gt.corners),
                       std::move(cam)};
}

int main(int argc, char **argv, char **env) {
  gui::Singleton::SetCmdArgs(argc, argv, env);
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");
  misc::Matlab matlab;

  auto config = ParseConfig("towers", "cam1");
  core::Println("nedges: ", config.topo.nedges());

  bool rerun_preprocess = false;
  bool rerun_vps = true;
  bool rerun_orientation_estimation = true;

  std::vector<std::set<int>> face_sets;
  std::vector<CameraParam> pp_focals;

  if (rerun_preprocess ||
      !config.loadCache("preprocess", face_sets, pp_focals)) {
    // decompose the drawing
    face_sets = DecomposeFaces(config.topo.face2corners, config.corners2d);
    // calibrate camera
    pp_focals =
        CalibrateCamera(BoundingBoxOfContainer(config.corners2d), face_sets,
                        [&config](int face) {
                          Chain2 face_loop;
                          auto &vs = config.topo.face2corners[face];
                          face_loop.points.reserve(vs.size());
                          for (int v : vs) {
                            face_loop.append(config.corners2d[v]);
                          }
                          return face_loop;
                        },
                        5);
    config.saveCache("preprocess", face_sets, pp_focals);
  }

  std::vector<Line2> edge2line(config.topo.nedges());
  for (int edge = 0; edge < config.topo.nedges(); edge++) {
    edge2line[edge].first =
        config.corners2d[config.topo.edge2corners[edge].first];
    edge2line[edge].second =
        config.corners2d[config.topo.edge2corners[edge].second];
  }
  auto lines_box = BoundingBoxOfContainer(edge2line);

  // pick the best!
  assert(!pp_focals.empty());
  auto &pp_focal = pp_focals.front();
  core::Println("current focal = ", pp_focal.focal, " pp = ", pp_focal.pp);
  PerspectiveCamera cam(lines_box.size(0), lines_box.size(1), pp_focal.pp,
                        pp_focal.focal);

  // merge lines
  std::vector<int> edge2merged_line(edge2line.size());
  std::vector<Line2> merged_lines = MergeColinearLines(
      edge2line, pp_focal, DegreesToRadians(0.1), &edge2merged_line);

  // show merged lines
  if (true) {
    Image3ub im(config.cameraGT.screenSize(), Vec3ub(255, 255, 255));
    auto canvas = gui::MakeCanvas(im);
    canvas.color(gui::LightGray);
    canvas.colorTable(gui::CreateRandomColorTableWithSize(merged_lines.size()));
    canvas.thickness(2);
    for (int i = 0; i < merged_lines.size(); i++) {
      canvas.add(ClassifyAs(merged_lines[i], i));
    }
    // show noted corners
    canvas.color(gui::Black);
    canvas.paintingOptions().fontScale = 1.0;
    for (int i = 0; i < config.corners2d.size(); i++) {
      canvas.add(NoteAs(config.corners2d[i], std::to_string(i)));
    }
    canvas.color(gui::Blue);
    canvas.paintingOptions().fontScale = 0.9;
    for (int i = 0; i < config.topo.nedges(); i++) {
      canvas.add(NoteAs((config.corners2d[config.topo.edge2corners[i].first] +
                         config.corners2d[config.topo.edge2corners[i].second]) /
                            2.0,
                        "[" + std::to_string(i) + "]"));
    }
    canvas.show(0, "merged lines");
  }

  // collect vanishing points using merged lines
  std::vector<Point2> vps;
  if (rerun_vps || !config.loadCache("vps", vps)) {
    CollectVanishingPointsParam param;
    // param.angle_thres_phase1 = DegreesToRadians(2);
    // param.angle_thres_phase2 = DegreesToRadians(0.2);
    // param.angle_thres_phase3 = DegreesToRadians(1);
    param.use_mean_shift_merge_phase1 = true;
    for (int k = 0; k < 1; k++) {
      core::Println(k, "-th call of CollectVanishingPoints");
      vps = CollectVanishingPoints(merged_lines, pp_focal.focal, pp_focal.pp,
                                   param);
      if (vps.size() < 500) {
        break;
      }
      param.angle_thres_phase1 *= 1.2;
      param.angle_thres_phase3 *= 0.8;
    }
    config.saveCache("vps", vps);
  }

  // show vanishing points
  if (false) {
    auto vp2lines = BindPointsToLines(vps, edge2line, DegreesToRadians(8));
    for (int i = 0; i < vps.size(); i++) {
      if (vp2lines[i].empty()) {
        continue;
      }
      Image3ub im(config.cameraGT.screenSize(), Vec3ub(255, 255, 255));
      auto canvas = gui::MakeCanvas(im);
      canvas.color(gui::LightGray);
      canvas.thickness(2);
      for (auto &line : edge2line) {
        canvas.add(line);
      }
      canvas.color(gui::Gray);
      canvas.thickness(2);
      for (int edge : vp2lines[i]) {
        canvas.add(edge2line[edge].ray());
      }
      canvas.color(gui::Black);
      for (int edge : vp2lines[i]) {
        canvas.add(edge2line[edge]);
      }
      canvas.show(0, "raw vp_" + std::to_string(i));
    }
  }

  // build relation of face to merged lines
  std::vector<std::vector<int>> face2merged_lines(config.topo.nfaces());
  for (int face = 0; face < config.topo.nfaces(); face++) {
    auto &edges = config.topo.face2edges[face];
    auto &mlines = face2merged_lines[face];
    for (int edge : config.topo.face2edges[face]) {
      int mline = edge2merged_line[edge];
      if (!mlines.empty() && mlines.back() == mline) {
        continue;
      }
      mlines.push_back(mline);
    }
  }

  // estimate orientations of merged lines using vps
  std::vector<int> merged_line2vp;
  if (rerun_orientation_estimation ||
      !config.loadCache("merged_line2vp", merged_line2vp)) {
    EstimateEdgeOrientationsParam param;
    //param.angle_thres_allowed_vp_line_deviation = DegreesToRadians(5);
    merged_line2vp =
        EstimateEdgeOrientations(merged_lines, vps, face2merged_lines,
                                 pp_focal.focal, pp_focal.pp, param);
    config.saveCache("merged_line2vp", merged_line2vp);
  }

  std::vector<int> edge2vp(config.topo.nedges(), -1);

  // collect corner directions
  std::vector<Vec3> corner2dir(config.corners2d.size());
  for (int i = 0; i < config.corners2d.size(); i++) {
    corner2dir[i] = normalize(cam.direction(config.corners2d[i]));
  }

  std::vector<int> vps_in_use;
  for (int edge = 0; edge < edge2vp.size(); edge++) {
    edge2vp[edge] = merged_line2vp[edge2merged_line[edge]];
    if (edge2vp[edge] != -1 && !Contains(vps_in_use, edge2vp[edge])) {
      vps_in_use.push_back(edge2vp[edge]);
    }
  }

  // show line orientation estimation
  if (true) {
    std::vector<std::set<int>> vp2lines(vps.size());
    for (int l = 0; l < edge2vp.size(); l++) {
      int vp = edge2vp[l];
      if (vp != -1) {
        vp2lines[vp].insert(l);
      }
    }
    for (int i = 0; i < vps.size(); i++) {
      if (vp2lines[i].empty()) {
        continue;
      }
      Image3ub im(config.cameraGT.screenSize(), Vec3ub(255, 255, 255));
      auto canvas = gui::MakeCanvas(im);
      canvas.color(gui::LightGray);
      canvas.thickness(2);
      for (auto &line : edge2line) {
        canvas.add(line);
      }
      canvas.color(gui::Gray);
      canvas.thickness(2);
      for (int edge : vp2lines[i]) {
        canvas.add(edge2line[edge].ray());
      }
      canvas.color(gui::Black);
      for (int edge : vp2lines[i]) {
        canvas.add(edge2line[edge]);
      }
      canvas.show(0, "optimized vp_" + std::to_string(i));
    }
  }

  // basic_constraints
  std::vector<PlaneConstraint> basic_constraints;
  {
    // add face coplanarities
    for (int face = 0; face < config.topo.nfaces(); face++) {
      basic_constraints.push_back(
          PlaneConstraint{config.topo.face2corners[face], MakePlaneMatrix()});
    }
    // add colinear edges
    std::vector<std::set<int>> colinear_edge_groups(edge2vp.size());
    {
      colinear_edge_groups.resize(ConnectedComponents(
          MakeIotaIterator<int>(0), MakeIotaIterator<int>(edge2vp.size()),
          [&config, &edge2merged_line](int edge) {
            std::set<int> adj_colinear_edges;
            if (edge2merged_line[edge] == -1) {
              return adj_colinear_edges;
            }
            int corners[] = {config.topo.edge2corners[edge].first,
                             config.topo.edge2corners[edge].second};
            for (int corner : corners) {
              for (int adj_edge : config.topo.corner2edges[corner]) {
                if (adj_edge == edge) {
                  continue;
                }
                if (edge2merged_line[adj_edge] == edge2merged_line[edge]) {
                  adj_colinear_edges.insert(adj_edge);
                }
              }
            }
            return adj_colinear_edges;
          },
          [&colinear_edge_groups](int edge, int group) {
            colinear_edge_groups[group].insert(edge);
          }));
    }
    for (auto &colinear_edges : colinear_edge_groups) {
      if (colinear_edges.size() <= 1) {
        continue;
      }
      // get colinear corners
      std::set<int> colinear_corners;
      for (int edge : colinear_edges) {
        colinear_corners.insert({config.topo.edge2corners[edge].first,
                                 config.topo.edge2corners[edge].second});
      }
      if (colinear_corners.size() < 3) {
        continue;
      }
      auto &merged_line =
          merged_lines[edge2merged_line[*colinear_edges.begin()]];
      Vec3 line_perp_dir = cam.direction(merged_line.first)
                               .cross(cam.direction(merged_line.second));
      basic_constraints.push_back(PlaneConstraint{
          std::vector<int>(colinear_corners.begin(), colinear_corners.end()),
          MakePlaneMatrixAlongDirection(line_perp_dir)});
    }
  }

  // BeamSearch
  Weighted<std::vector<Point3>> best_result;
  best_result.weight() = std::numeric_limits<double>::infinity();

  std::default_random_engine rng;

  auto beam_search_energy_fun =
      [&basic_constraints, &edge2vp, &config, &face_sets, &corner2dir,
       &edge2line, &cam, &vps, &vps_in_use, &best_result,
       &rng](const std::vector<bool> &vps_in_use_disable_flags) -> double {
    
    std::set<int> enabled_vps_in_use;
    for (int i = 0; i < vps_in_use.size(); i++) {
      if (!vps_in_use_disable_flags[i]) {
        enabled_vps_in_use.insert(vps_in_use[i]);
      }
    }

    std::string msg = "searching with enabled_vps_in_use = {";
    for (int vp : enabled_vps_in_use) {
      msg += std::to_string(vp) + " ";
    }
    msg += "}";
    misc::Clock clock(msg);

    // show current configured line orientations
    if (false) {
      std::vector<std::set<int>> vp2lines(vps.size());
      for (int l = 0; l < edge2vp.size(); l++) {
        int vp = edge2vp[l];
        if (Contains(enabled_vps_in_use, vp)) {
          vp2lines[vp].insert(l);
        }
      }
      for (int i = 0; i < vps.size(); i++) {
        if (vp2lines[i].empty()) {
          continue;
        }
        Image3ub im(config.cameraGT.screenSize(), Vec3ub(255, 255, 255));
        auto canvas = gui::MakeCanvas(im);
        canvas.color(gui::LightGray);
        canvas.thickness(2);
        for (auto &line : edge2line) {
          canvas.add(line);
        }
        canvas.color(gui::Gray);
        canvas.thickness(2);
        for (int edge : vp2lines[i]) {
          canvas.add(edge2line[edge].ray());
        }
        canvas.color(gui::Black);
        for (int edge : vp2lines[i]) {
          canvas.add(edge2line[edge]);
        }
        canvas.show(0, "enabled vp_" + std::to_string(i));
      }
    }

    std::vector<PlaneConstraint> constraints = basic_constraints;
    // add oriented edges
    for (int edge = 0; edge < edge2vp.size(); edge++) {
      if (Contains(enabled_vps_in_use, edge2vp[edge])) { //
        auto &line = edge2line[edge];
        Vec3 vp = cam.direction(vps[edge2vp[edge]]);
        Vec3 line_perp_dir =
            cam.direction(line.first).cross(cam.direction(line.second));
        Vec3 normal = normalize(vp.cross(line_perp_dir));
        constraints.push_back(
            PlaneConstraint{{config.topo.edge2corners[edge].first,
                             config.topo.edge2corners[edge].second},
                            MakePlaneMatrixTowardDirection(normal)});
      }
    }

    // analyze the dependency matrices of verts and constraints
    static const int repeat_num = 3;
    std::vector<Weighted<std::vector<Point3>>>
        results_table_under_cur_constraints(repeat_num);

    // energy terms
    auto energy_fun = [&config, &face_sets, &edge2line](
        const Inferencer &infer, const DenseMatd &vars,
        const std::vector<Vec3> &vert2dir) -> double {
      DenseMatd variables = cv::abs(vars) / cv::norm(vars);
      // edge msda & edge ortho
      auto edge_angles = AnglesBetweenAdjacentEdges(
          vert2dir, config.topo.face2corners, variables, infer);
      double edge_msda = MeanSquaredDeviationOfContainer(edge_angles);
      double edge_ortho =
          MeanSquaredDeviationOfContainer(edge_angles, {0.0, M_PI_2});

      // face msda & face ortho
      double face_msda =
          MeanSquaredDeviationOfContainer(AnglesBetweenAdjacentFaces(
              config.topo.nfaces(), config.topo.edge2faces, variables, infer));
      double face_ortho = MeanSquaredDeviationOfContainer(
          AnglesBetweenAdjacentFaces(config.topo.nfaces(),
                                     config.topo.edge2faces, variables, infer,
                                     config.faces_overlap),
          {M_PI_2, M_PI, M_PI_2 * 3});

      // edge_msda_each_face_set
      double edge_msda_each_face_set = 0.0;
      {
        std::vector<double> each_face_set(face_sets.size(), 0.0);
        for (int i = 0; i < face_sets.size(); i++) {
          each_face_set[i] =
              MeanSquaredDeviationOfContainer(AnglesBetweenAdjacentEdges(
                  vert2dir, config.topo.face2corners, variables, infer,
                  [i, &face_sets](int face) {
                    return Contains(face_sets[i], face);
                  }));
        }
        edge_msda_each_face_set =
            std::accumulate(each_face_set.begin(), each_face_set.end(), 0.0) /
            face_sets.size();
      }

      // face_ortho_each_face_set
      double face_ortho_each_face_set = 0.0;
      {
        std::vector<double> each_face_set(face_sets.size(), 0.0);
        for (int i = 0; i < face_sets.size(); i++) {
          each_face_set[i] = MeanSquaredDeviationOfContainer(
              AnglesBetweenAdjacentFaces(config.topo.nfaces(),
                                         config.topo.edge2faces, variables,
                                         infer, config.faces_overlap,
                                         [i, &face_sets](int face) {
                                           return Contains(face_sets[i], face);
                                         }),
              {M_PI_2, M_PI, M_PI_2 * 3});
        }
        face_ortho_each_face_set =
            std::accumulate(each_face_set.begin(), each_face_set.end(), 0.0) /
            face_sets.size();
      }

      // edge skew ratio
      double max_edge_skew_ratio_score = 0.0;
      {
        std::vector<double> edge_proj_ratios(config.topo.nedges());
        for (int edge = 0; edge < config.topo.nedges(); edge++) {
          double length2d = edge2line[edge].length();
          auto &edge_corners = config.topo.edge2corners[edge];
          double length3d =
              Distance(normalize(vert2dir[edge_corners.first]) /
                           infer.getInversedDepth(edge_corners.first, vars),
                       normalize(vert2dir[edge_corners.second]) /
                           infer.getInversedDepth(edge_corners.second, vars));
          edge_proj_ratios[edge] = length3d / length2d;
        }
        std::sort(edge_proj_ratios.begin(), edge_proj_ratios.end());
        max_edge_skew_ratio_score =
            edge_proj_ratios.back() /
            edge_proj_ratios[edge_proj_ratios.size() / 2];
      }

      double e = 0.0;
      //{ // classic params
      //  e = 100.0 * edge_msda + 0.1 * face_msda + 10 * edge_ortho +
      //      20 * face_ortho + 100.0 * edge_msda_each_face_set +
      //      100 * face_ortho_each_face_set;
      //}
      {
        e = 100.0 * edge_msda + 0 * face_msda + 0 * edge_ortho +
            100 * face_ortho + 100.0 * edge_msda_each_face_set +
            100 * face_ortho_each_face_set +
            std::max(max_edge_skew_ratio_score - 2.5, 0.0) * 100;
      }
      /* core::Println("e = ", e, " edge_msda = ", edge_msda, " face_msda = ",
                     face_msda, " edge_ortho = ", edge_ortho, " face_ortho =
         ",
                     face_ortho, " edge_msda_each_face_set = ",
         edge_msda_each_face_set,
                     " face_ortho_each_face_set = ",
         face_ortho_each_face_set);*/
      return e;
    };
    size_t disabled_vps_num = vps_in_use.size() - enabled_vps_in_use.size();
    ParallelRun(repeat_num, std::thread::hardware_concurrency() - 1,
                [&constraints, &corner2dir, &rng, &energy_fun, disabled_vps_num,
                 &results_table_under_cur_constraints](int repeat_id) {
                  std::vector<int> fundamental_verts;
                  std::vector<Point3> cur_corner2position;
                  PerformReconstructionParam param;
                  param.max_iters = 100;
                  double cur_energy = PerformReconstruction(
                                          constraints, corner2dir,
                                          std::uniform_int_distribution<int>(
                                              0, corner2dir.size() - 1)(rng),
                                          energy_fun, rng, cur_corner2position,
                                          &fundamental_verts, param) +
                                      disabled_vps_num * 0.1; // compensation
                  // fill cur_energy and cur_corner2position to the tables
                  results_table_under_cur_constraints[repeat_id].component =
                      std::move(cur_corner2position);
                  results_table_under_cur_constraints[repeat_id].weight() =
                      cur_energy;
                });

    auto &best_result_under_cur_constraints =
        *std::min_element(results_table_under_cur_constraints.begin(),
                          results_table_under_cur_constraints.end());
    core::Println("best energy in this beam tree node = ",
                  best_result_under_cur_constraints.weight());
    if (best_result_under_cur_constraints < best_result) {
      best_result = best_result_under_cur_constraints;
    }

    if (false) { // show current node result
      gui::SceneBuilder sb;
      for (int e = 0; e < config.topo.nedges(); e++) {
        auto &p1 = best_result_under_cur_constraints
                       .component[config.topo.edge2corners[e].first];
        auto &p2 = best_result_under_cur_constraints
                       .component[config.topo.edge2corners[e].second];
        sb.add(Line3(p1, p2));
      }
      sb.show(true, true,
              gui::RenderOptions()
                  .renderMode(gui::Lines)
                  .fixUpDirectionInCameraMove(false)
                  .winName("best result in current searching node"));
    }

    return best_result_under_cur_constraints.weight();
  };

  auto best_vps_in_use_disable_flags =
      BeamSearch(vps_in_use.size(), beam_search_energy_fun, 1);

  core::Println("final energy after beam searching = ", best_result.weight());

  if (true) { // show final result
    gui::SceneBuilder sb;
    for (int e = 0; e < config.topo.nedges(); e++) {
      auto &p1 = best_result.component[config.topo.edge2corners[e].first];
      auto &p2 = best_result.component[config.topo.edge2corners[e].second];
      sb.add(Line3(p1, p2));
    }
    sb.show(true, true, gui::RenderOptions()
                            .renderMode(gui::Lines)
                            .fixUpDirectionInCameraMove(false)
                            .winName("final result"));
  }

  return 0;
}