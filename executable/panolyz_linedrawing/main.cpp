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

struct LineDrawingInput {
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

// ParseInput
LineDrawingInput ParseInput(const std::string &modelName,
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
  return LineDrawingInput{modelName,
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

  auto input = ParseInput("1", "cam1");
  Println("nedges: ", input.topo.nedges());


  bool rerun_preprocess = false;
  bool rerun_vps = false;
  bool rerun_orientation_estimation = false;

  std::vector<std::set<int>> face_sets;
  std::vector<CameraParam> pp_focals;

  if (rerun_preprocess ||
      !input.loadCache("preprocess", face_sets, pp_focals)) {
    // decompose the drawing
    face_sets = DecomposeFaces(input.topo.face2corners, input.corners2d);

    // calibrate camera
    pp_focals =
        CalibrateCamera(BoundingBoxOfContainer(input.corners2d), face_sets,
                        [&input](int face) {
                          Chain2 face_loop;
                          auto &vs = input.topo.face2corners[face];
                          face_loop.points.reserve(vs.size());
                          for (int v : vs) {
                            face_loop.append(input.corners2d[v]);
                          }
                          return face_loop;
                        },
                        5);
    input.saveCache("preprocess", face_sets, pp_focals);
  }

  std::vector<Line2> edge2line(input.topo.nedges());
  for (int edge = 0; edge < input.topo.nedges(); edge++) {
    edge2line[edge].first =
        input.corners2d[input.topo.edge2corners[edge].first];
    edge2line[edge].second =
        input.corners2d[input.topo.edge2corners[edge].second];
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
    Image3ub im(input.cameraGT.screenSize(), Vec3ub(255, 255, 255));
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
  if (rerun_vps || !input.loadCache("vps", vps)) {
    vps = CollectVanishingPoints(merged_lines, pp_focal.focal, pp_focal.pp);
    input.saveCache("vps", vps);
  }

  // show vanishing points
  if (false) {
    auto vp2lines = BindPointsToLines(vps, edge2line, DegreesToRadians(8));
    for (int i = 0; i < vps.size(); i++) {
      if (vp2lines[i].empty()) {
        continue;
      }
      Image3ub im(input.cameraGT.screenSize(), Vec3ub(255, 255, 255));
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
  std::vector<std::vector<int>> face2merged_lines(input.topo.nfaces());
  for (int face = 0; face < input.topo.nfaces(); face++) {
    auto &edges = input.topo.face2edges[face];
    auto &mlines = face2merged_lines[face];
    for (int edge : input.topo.face2edges[face]) {
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
      !input.loadCache("merged_line2vp", merged_line2vp)) {
    auto merged_line2vp = EstimateEdgeOrientations(
        merged_lines, vps, face2merged_lines, pp_focal.focal, pp_focal.pp);
    input.saveCache("merged_line2vp", merged_line2vp);
  }

  std::vector<int> edge2vp(input.topo.nedges(), -1);
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
      Image3ub im(input.cameraGT.screenSize(), Vec3ub(255, 255, 255));
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
  std::vector<Vec3> corner2dir(input.corners2d.size());
  for (int i = 0; i < input.corners2d.size(); i++) {
    corner2dir[i] = normalize(cam.direction(input.corners2d[i]));
  }

  while (true) {
    // collect constraints
    std::vector<PlaneConstraint> constraints;
    // add face coplanarities
    for (int face = 0; face < input.topo.nfaces(); face++) {
      constraints.push_back(
          PlaneConstraint{input.topo.face2corners[face], MakePlaneMatrix()});
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
            PlaneConstraint{{input.topo.edge2corners[edge].first,
                             input.topo.edge2corners[edge].second},
                            MakePlaneMatrixTowardDirection(normal)});
      }
    }

    // analyze the dependency matrices of verts and constraints
    auto infer = GenerateInferenceFunctors(constraints, corner2dir);

    // reconstruct
  }
  return 0;
}