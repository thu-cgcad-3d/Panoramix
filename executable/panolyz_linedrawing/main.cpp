#include "cache.hpp"
#include "cameras.hpp"
#include "canvas.hpp"
#include "line_drawing.hpp"
#include "scene.hpp"
#include "singleton.hpp"

#include "tools.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

struct Configuration {
  std::string model_name, cam_name;
  LineDrawingTopo topo;
  std::vector<Point2> corners2d;
  std::vector<Point3> cornersGT;
  PerspectiveCamera cameraGT;

  template <class ... Ts>
  bool loadCache(const std::string & what, Ts & ... ts) const {
    return misc::LoadCache(model_name + "_" + cam_name, what, ts ...);
  }
  template <class ... Ts>
  bool saveCache(const std::string & what, Ts & ... ts) const {
    return misc::SaveCache(model_name + "_" + cam_name, what, ts ...);
  }
};

// ParseConfig
Configuration ParseConfig(const std::string &modelName,
                            const std::string &camName) {
  std::string folder = "F:\\LineDrawings\\nonmanifold\\" + modelName + "\\";

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

  Println("gt focal = ", cam.focal(), " gt pp = ", cam.principlePoint());
  std::vector<Point2> corners2d(line_drawing_gt.corners.size());
  for (int i = 0; i < corners2d.size(); i++) {
    corners2d[i] = cam.toScreen(line_drawing_gt.corners[i]);
  }
  return Configuration{modelName,
                          camName,
                          std::move(line_drawing_gt.topo),
                          std::move(corners2d),
                          std::move(line_drawing_gt.corners),
                          std::move(cam)};
}

int main(int argc, char **argv, char **env) {
  gui::Singleton::SetCmdArgs(argc, argv, env);
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");
  misc::Matlab matlab;

  auto config = ParseConfig("2", "cam1");
  Println("nedges: ", config.topo.nedges());


  bool rerun_preprocess = false;
  bool rerun_vps = false;
  bool rerun_orientation_estimation = false;

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
  Println("current focal = ", pp_focal.focal, " pp = ", pp_focal.pp);
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
    canvas.show(0, "merged lines");
  }

  // collect vanishing points using merged lines
  std::vector<Point2> vps;
  if (rerun_vps || !config.loadCache("vps", vps)) {
    vps = CollectVanishingPoints(merged_lines, pp_focal.focal, pp_focal.pp);
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
    merged_line2vp = EstimateEdgeOrientations(
        merged_lines, vps, face2merged_lines, pp_focal.focal, pp_focal.pp);
    config.saveCache("merged_line2vp", merged_line2vp);
  }

  std::vector<int> edge2vp(config.topo.nedges(), -1);
  for (int edge = 0; edge < edge2vp.size(); edge++) {
    edge2vp[edge] = merged_line2vp[edge2merged_line[edge]];
  }

  // show line orientation results
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

  // collect corner directions
  std::vector<Vec3> corner2dir(config.corners2d.size());
  for (int i = 0; i < config.corners2d.size(); i++) {
    corner2dir[i] = normalize(cam.direction(config.corners2d[i]));
  }

  while (true) {
    // collect constraints
    std::vector<PlaneConstraint> constraints;
    // add face coplanarities
    for (int face = 0; face < config.topo.nfaces(); face++) {
      constraints.push_back(
          PlaneConstraint{config.topo.face2corners[face], MakePlaneMatrix()});
    }
    // add oriented edges
    for (int edge = 0; edge < edge2vp.size(); edge++) {
      if (edge2vp[edge] != -1) {
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
      auto & merged_line = merged_lines[edge2merged_line[*colinear_edges.begin()]];
      Vec3 line_perp_dir = cam.direction(merged_line.first)
                               .cross(cam.direction(merged_line.second));
      constraints.push_back(PlaneConstraint{
          std::vector<int>(colinear_corners.begin(), colinear_corners.end()),
          MakePlaneMatrixAlongDirection(line_perp_dir)});
    }

    // analyze the dependency matrices of verts and constraints
    std::vector<int> fundamental_verts;
    auto infer =
        GenerateInferenceFunctors(constraints, corner2dir, &fundamental_verts);

    if (true) { // show fundamental verts
      Image3ub im(config.cameraGT.screenSize(), Vec3ub(255, 255, 255));
      auto canvas = gui::MakeCanvas(im);
      canvas.color(gui::LightGray);
      canvas.thickness(2);
      canvas.add(edge2line);
      canvas.color(gui::Red);
      for (int v : fundamental_verts) {
        canvas.add(Circle{config.corners2d[v], 2});
      }
      canvas.show(0, "fundamental verts");
    }

    // energy terms
    auto energy_fun = [&infer, &corner2dir, &config,
                       &face_sets](const std::vector<double> &vars) -> double {
      DenseMatd variables = cv::abs(DenseMatd(vars, false)) / cv::norm(vars);
      auto edge_angles = AnglesBetweenAdjacentEdges(
          corner2dir, config.topo.face2corners, variables, *infer);
      auto face_angles = AnglesBetweenAdjacentFaces(
          config.topo.nfaces(), config.topo.edge2faces, variables, *infer);

      double edge_msda = MeanSquaredDeviationOfContainer(edge_angles);
      double face_msda = MeanSquaredDeviationOfContainer(face_angles);
      double edge_ortho =
          MeanSquaredDeviationOfContainer(edge_angles, {M_PI_2, 0.0, M_PI});
      double face_ortho =
          MeanSquaredDeviationOfContainer(face_angles, {M_PI_2, 0.0, M_PI});

      double edge_msda_face_set = 0.0;
      {
        std::vector<double> each_face_set(face_sets.size(), 0.0);
        for (int i = 0; i < face_sets.size(); i++) {
          each_face_set[i] =
              MeanSquaredDeviationOfContainer(AnglesBetweenAdjacentEdges(
                  corner2dir, config.topo.face2corners, variables, *infer,
                  [i, &face_sets](int face) {
                    return Contains(face_sets[i], face);
                  }));
        }
        edge_msda_face_set =
            std::accumulate(each_face_set.begin(), each_face_set.end(), 0.0) /
            face_sets.size();
      }

      double face_ortho_face_set = 0.0;
      {
        std::vector<double> each_face_set(face_sets.size(), 0.0);
        for (int i = 0; i < face_sets.size(); i++) {
          each_face_set[i] = MeanSquaredDeviationOfContainer(
              AnglesBetweenAdjacentFaces(config.topo.nfaces(),
                                         config.topo.edge2faces, variables,
                                         *infer,
                                         [i, &face_sets](int face) {
                                           return Contains(face_sets[i], face);
                                         }),
              {M_PI_2, 0.0, M_PI});
        }
        face_ortho_face_set =
            std::accumulate(each_face_set.begin(), each_face_set.end(), 0.0) /
            face_sets.size();
      }

      double e = 100.0 * edge_msda + 0.1 * face_msda + 10 * edge_ortho +
                 20 * face_ortho + 100.0 * edge_msda_face_set +
                 100 * face_ortho_face_set;
      Println("e = ", e, " edge_msda = ", edge_msda, " face_msda = ", face_msda,
              " edge_ortho = ", edge_ortho, " face_ortho = ", face_ortho,
              " edge_msda_face_set = ", edge_msda_face_set,
              " face_ortho_face_set = ", face_ortho_face_set);
      return e;
    };

    // reconstruct
    std::default_random_engine rng;
    std::vector<double> vars(infer->nvars(), 1);
    SimulatedAnnealing(
        vars, energy_fun,
        [](int iter) { return std::max(1.0 / log(log(log(iter + 2))), 1e-10); },
        [](std::vector<double> curX, int iter, auto &&forEachNeighborFun) {
          if (iter >= 100) {
            return;
          }
          double step = std::max(1.0 / (iter + 2), 1e-10);
          for (int i = 0; i < curX.size(); i++) {
            double curXHere = curX[i];
            curX[i] = curXHere + step;
            forEachNeighborFun(curX);
            curX[i] = curXHere - step;
            forEachNeighborFun(curX);
            curX[i] = curXHere;
          }
        },
        rng);

    std::vector<double> corner2gt_inv_depths(config.corners2d.size(), 1.0);
    for (int c = 0; c < config.corners2d.size(); c++) {
      corner2gt_inv_depths[c] = 1.0 / norm(config.corners2d[c]);
    }
    DenseMatd ground_truth_vars_mat =
        infer->recoverVariables(corner2gt_inv_depths);
    Println("groundtruth energy:");
    energy_fun(ground_truth_vars_mat);

    if (true) { // show results
      DenseMatd variables = cv::abs(DenseMatd(vars, false)) / cv::norm(vars);
      std::vector<Point3> final_corner2position(corner2dir.size());
      for (int c = 0; c < corner2dir.size(); c++) {
        final_corner2position[c] =
            normalize(corner2dir[c]) / infer->getInversedDepth(c, variables);
      }
      gui::SceneBuilder sb;
      for (int e = 0; e < config.topo.nedges(); e++) {
        auto &p1 = final_corner2position[config.topo.edge2corners[e].first];
        auto &p2 = final_corner2position[config.topo.edge2corners[e].second];
        sb.add(Line3(p1, p2));
      }
      sb.show(true, true, gui::RenderOptions()
                              .renderMode(gui::Lines)
                              .fixUpDirectionInCameraMove(false));
    }
  }
  return 0;
}