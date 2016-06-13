#include "cache.hpp"
#include "cameras.hpp"
#include "canvas.hpp"
#include "clock.hpp"
#include "line_drawing.hpp"
#include "line_drawing_widget.hpp"
#include "optimization.hpp"
#include "parallel.hpp"
#include "scene.hpp"
#include "singleton.hpp"

#include "tools.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

struct Configuration {
  std::string model_name, cam_name;
  LineDrawingAnnotation anno;

  struct GroundTruth {
    std::vector<Point3> points;
    PerspectiveCamera camera;
  };
  std::unique_ptr<GroundTruth> gt;

  template <class... Ts>
  bool loadCache(const std::string &what, Ts &... ts) const {
    return misc::LoadCache(model_name + "_" + cam_name, what, ts...);
  }
  template <class... Ts>
  bool saveCache(const std::string &what, Ts &... ts) const {
    return misc::SaveCache(model_name + "_" + cam_name, what, ts...);
  }

  // FromObjFile
  static Configuration FromObjFile(const std::string &model_name,
                                   const std::string &cam_name) {
    std::string folder =
        PANORAMIX_TEST_DATA_DIR_STR "\\linedrawing\\objs\\" + model_name + "\\";

    // load 3d points and faces from obj file
    std::vector<Point3> point3ds;
    std::vector<std::vector<int>> face2corners;
    MakeMeshFromObjFile(
        [&point3ds](float x, float y, float z) {
          point3ds.emplace_back(x, y, z);
        },
        [&face2corners](auto &&face_corners) {
          face2corners.push_back(face_corners);
        },
        folder + model_name + "_w_intf.obj");

    assert(point3ds.size() > 0);

    // setup a camera
    PerspectiveCamera cam;
    std::string cam_path = folder + model_name + ".obj." + cam_name + ".cereal";
    if (!LoadFromDisk(cam_path, cam)) {
      gui::SceneBuilder sb;
      for (auto &f : face2corners) {
        assert(f.size() >= 3);
        for (int i = 0; i < f.size(); i++) {
          sb.add(Line3(point3ds[f[i]], point3ds[f[(i + 1) % f.size()]]));
        }
      }
      cam = sb.show(true, true, gui::RenderOptions()
                                    .renderMode(gui::Lines)
                                    .fixUpDirectionInCameraMove(false))
                .camera();
      SaveToDisk(cam_path, cam);
    }
    core::Println("gt focal = ", cam.focal(), " gt pp = ",
                  cam.principlePoint());

    // build anno data
    LineDrawingAnnotation anno;
    // points
    anno.points.resize(point3ds.size());
    for (int i = 0; i < anno.points.size(); i++) {
      anno.points[i] = cam.toScreen(point3ds[i]);
    }
    // edges
    std::set<std::pair<int, int>> edges;
    for (auto &f : face2corners) {
      assert(f.size() >= 3);
      for (int i = 0; i < f.size(); i++) {
        int c1 = f[i];
        int c2 = f[(i + 1) % f.size()];
        edges.insert(MakeOrderedPair(c1, c2));
      }
    }
    anno.edges = std::vector<std::pair<int, int>>(edges.begin(), edges.end());
    // coplanar_points
    anno.coplanar_points = face2corners;
    // image (white board)
    anno.image = Image3ub(cam.screenSize(), Vec3ub(255, 255, 255));

    // groudtruth data
    auto gt_ptr = std::make_unique<Configuration::GroundTruth>();
    gt_ptr->camera = std::move(cam);
    gt_ptr->points = std::move(point3ds);

    return Configuration{model_name, cam_name, std::move(anno),
                         std::move(gt_ptr)};
  }

  // FromImageAnnotation
  static Configuration FromImageAnnotation(const std::string &image_name) {
    static const std::string folder =
        PANORAMIX_TEST_DATA_DIR_STR "\\linedrawing\\images\\";
    std::string anno_path = folder + image_name + ".ldanno.cereal";
    LineDrawingAnnotation anno;
    bool succ = LoadFromDisk(anno_path, anno);
    if (!succ) {
      core::Println("Failed to load ", image_name, "!");
      return Configuration();
    }
    return Configuration{"image_" + image_name, "no_cam", std::move(anno),
                         nullptr};
  }
};

void RoutineReconstruct() {

  //auto config = Configuration::FromObjFile("hex", "cam1");
	auto config = Configuration::FromImageAnnotation("house_sketch.jpg");
  //auto config = Configuration::FromImageAnnotation("washington.jpg");

  const size_t npoints = config.anno.points.size();
  const size_t nedges = config.anno.edges.size();
  const size_t nfaces = config.anno.coplanar_points.size();
  core::Println("nedges: ", config.anno.edges.size());

  bool rerun_preprocess = false;
  bool rerun_vps = false;
  bool rerun_orientation_estimation = false;
	bool rerun_reconstruction = false;

  // compute additional data
  std::vector<std::set<int>> point2edges(npoints);
  std::map<std::pair<int, int>, int> points2edge;
  std::vector<std::set<int>> edge2faces(nedges);
  std::vector<std::set<int>> face2edges(nfaces);
  std::map<std::pair<int, int>, bool> faces_overlap;
  std::vector<Line2> edge2line(nedges);
  {
    for (int edge = 0; edge < nedges; edge++) {
      auto &e = config.anno.edges[edge];
      points2edge[e] = edge;
      point2edges[e.first].insert(edge);
      point2edges[e.second].insert(edge);
    }
    for (int face = 0; face < nfaces; face++) {
      auto &ps = config.anno.coplanar_points[face];
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

    // compute faces_overlap : whether two faces sharing a common edge overlap?
    // first for each edge of a face, we select a corner of the face which
    // satisfies:
    // 1. the corner is not on the edge
    // 2. the triangle spanned by the corner and the edge is covered by the face
    // polygon in 2d drawing
    std::vector<std::map<int, int>> face2far_corners_for_edges(nfaces);
    for (int face = 0; face < nfaces; face++) {
      auto &cs = config.anno.coplanar_points[face];
      auto &far_corners_for_edges = face2far_corners_for_edges[face];
      TriangulatePolygon(cs.begin(), cs.end(),
                         [&config](int corner) -> const Point2 & {
                           return config.anno.points.at(corner);
                         },
                         [&far_corners_for_edges, &config,
                          &points2edge](int c1, int c2, int c3) {
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
    // the we compute face overlapnesss
    for (int edge = 0; edge < nedges; edge++) {
      Line2 line(config.anno.points[config.anno.edges[edge].first],
                 config.anno.points[config.anno.edges[edge].second]);
      auto &fs = edge2faces.at(edge);
      if (fs.size() < 2) {
        continue;
      }
      for (auto iti = fs.begin(); iti != fs.end(); ++iti) {
        int far_corner_i = face2far_corners_for_edges.at(*iti).at(edge);
        bool on_left_side_of_edge_i = IsOnLeftSide(
            config.anno.points.at(far_corner_i), line.first, line.second);
        for (auto itj = std::next(iti); itj != fs.end(); ++itj) {
          int far_corner_j = face2far_corners_for_edges.at(*itj).at(edge);
          bool on_left_side_of_edge_j = IsOnLeftSide(
              config.anno.points.at(far_corner_j), line.first, line.second);
          faces_overlap[MakeOrderedPair(*iti, *itj)] =
              on_left_side_of_edge_i == on_left_side_of_edge_j;
        }
      }
    }

    // construct 2d lines for edges
    for (int edge = 0; edge < nedges; edge++) {
      edge2line[edge].first = config.anno.points[config.anno.edges[edge].first];
      edge2line[edge].second =
          config.anno.points[config.anno.edges[edge].second];
    }
  }

  auto lines_box = BoundingBoxOfContainer(edge2line);

  // decompose into face sets
  // and estimate pps & focals
  std::vector<std::set<int>> face_sets;
  std::vector<CameraParam> pp_focals;
  if (rerun_preprocess ||
      !config.loadCache("preprocess", face_sets, pp_focals)) {
    // decompose the drawing
    face_sets = DecomposeFaces(config.anno.coplanar_points, config.anno.points);
    // calibrate camera
    pp_focals = CalibrateCamera(
        BoundingBoxOfContainer(config.anno.points), face_sets,
        [&config, &points2edge](int face) {
          std::vector<std::vector<int>> pids = {{}};
          auto &vs = config.anno.coplanar_points[face];
          for (int v : vs) {
            if (pids.back().empty() ||
                Contains(points2edge, MakeOrderedPair(pids.back().back(), v))) {
              pids.back().push_back(v);
            } else {
              pids.push_back({v});
            }
          }
          bool first_last_connected =
              Contains(points2edge, MakeOrderedPair(vs.back(), vs.front()));
          if (pids.size() > 1 &&
              first_last_connected) { // merge the first and the last chain
            pids.back().insert(pids.back().end(), pids.front().begin(),
                               pids.front().end());
            pids.front().clear();
          }
          std::vector<Chain2> face_loops;
          for (auto &ps : pids) {
            if (ps.size() <= 2) {
              continue;
            }
            Chain2 chain;
            chain.closed = false;
            for (int p : ps) {
              chain.append(config.anno.points[p]);
            }
            face_loops.push_back(std::move(chain));
          }
          if (face_loops.size() == 1 && first_last_connected) {
            face_loops.front().closed = true;
          }
          for (auto &loop : face_loops) {
            if (!loop.closed) {
              core::Println("a non-closed loop of size ", loop.size());
            } else {
              assert(face_loops.size() == 1 &&
                     face_loops.front().size() == vs.size());
            }
          }
          return face_loops;
        },
        5);
    config.saveCache("preprocess", face_sets, pp_focals);
  }

  if (true) { // show individual face sets
    for (auto & fs : face_sets) {
      auto canvas = gui::MakeCanvas(config.anno.image);
      canvas.colorTable(gui::CreateRandomColorTableWithSize(nfaces));

    }
  }


  // pick the top pp focal as the camera parameter
  assert(!pp_focals.empty());
  auto &pp_focal = pp_focals.front();
  core::Println("current focal = ", pp_focal.focal, " pp = ", pp_focal.pp);
  PerspectiveCamera cam(config.anno.image.cols, config.anno.image.rows,
                        pp_focal.pp, pp_focal.focal);

  // merge lines
  std::vector<int> edge2merged_line(edge2line.size());
  std::vector<Line2> merged_lines = MergeColinearLines(
      edge2line, pp_focal, DegreesToRadians(0.1), &edge2merged_line);

  // show merged lines
  if (true) {
    auto canvas = gui::MakeCanvas(config.anno.image);
    canvas.color(gui::LightGray);
    canvas.colorTable(gui::CreateRandomColorTableWithSize(merged_lines.size()));
    canvas.thickness(2);
    for (int i = 0; i < merged_lines.size(); i++) {
      canvas.add(ClassifyAs(merged_lines[i], i));
    }
    // show noted corners
    canvas.color(gui::Black);
    canvas.paintingOptions().fontScale = 1.0;
    for (int i = 0; i < npoints; i++) {
      canvas.add(NoteAs(config.anno.points[i], std::to_string(i)));
    }
    // show noted edges
    canvas.color(gui::Blue);
    canvas.paintingOptions().fontScale = 0.9;
    for (int i = 0; i < nedges; i++) {
      canvas.add(NoteAs((config.anno.points[config.anno.edges[i].first] +
                         config.anno.points[config.anno.edges[i].second]) /
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
      auto canvas = gui::MakeCanvas(config.anno.image);
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
  std::vector<std::vector<int>> face2merged_lines(nfaces);
  for (int face = 0; face < nfaces; face++) {
    auto &edges = face2edges[face];
    auto &mlines = face2merged_lines[face];
    for (int edge : edges) {
      int mline = edge2merged_line[edge];
      if (Contains(mlines, mline)) {
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
    // param.angle_thres_allowed_vp_line_deviation = DegreesToRadians(5);
    merged_line2vp =
        EstimateEdgeOrientations(merged_lines, vps, face2merged_lines,
                                 pp_focal.focal, pp_focal.pp, param);
    config.saveCache("merged_line2vp", merged_line2vp);
  }

  // collect point directions
  std::vector<Vec3> point2dir(npoints);
  for (int i = 0; i < npoints; i++) {
    point2dir[i] = normalize(cam.direction(config.anno.points[i]));
  }

  // edge2vp and vp2edges
  std::vector<int> edge2vp(nedges, -1);
  std::map<int, std::set<int>> vp2edges;
  for (int edge = 0; edge < edge2vp.size(); edge++) {
    int vp = merged_line2vp[edge2merged_line[edge]];
    edge2vp[edge] = vp;
    if (vp != -1) {
      vp2edges[vp].insert(edge);
    }
  }

  // show line orientation estimation
  if (true) {
    for (auto &vp_edges : vp2edges) {
      auto canvas = gui::MakeCanvas(config.anno.image);
      canvas.color(gui::LightGray);
      canvas.thickness(2);
      for (auto &line : edge2line) {
        canvas.add(line);
      }
      canvas.color(gui::Gray);
      canvas.thickness(2);
      for (int edge : vp_edges.second) {
        canvas.add(edge2line[edge].ray());
      }
      canvas.color(gui::Black);
      for (int edge : vp_edges.second) {
        canvas.add(edge2line[edge]);
      }
      canvas.show(0, "optimized vp_" + std::to_string(vp_edges.first));
    }
  }

  // basic_constraints
  std::vector<PlaneConstraint> basic_constraints;
  {
    // add face coplanarities
    for (int face = 0; face < nfaces; face++) {
      basic_constraints.push_back(PlaneConstraint{
          config.anno.coplanar_points[face], MakePlaneMatrix()});
    }
    // add colinearity from assignment
    for (auto &ps : config.anno.colinear_points) {
      assert(ps.size() >= 3);
      auto dir1 = cam.direction(config.anno.points[ps.front()]);
      auto dir2 = cam.direction(config.anno.points[ps.back()]);
      Vec3 line_perp_dir = dir1.cross(dir2);
      assert((!HasValue(line_perp_dir, [](auto e){return IsInfOrNaN(e);})));
      basic_constraints.push_back(
          PlaneConstraint{std::vector<int>(ps.begin(), ps.end()),
                          MakePlaneMatrixAlongDirection(line_perp_dir)});
    }
    // add colinearity from detection
    std::vector<std::set<int>> colinear_edge_groups(edge2vp.size());
    colinear_edge_groups.resize(ConnectedComponents(
        MakeIotaIterator<int>(0), MakeIotaIterator<int>(edge2vp.size()),
        [&config, &edge2merged_line, &point2edges](int edge) {
          std::set<int> adj_colinear_edges;
          if (edge2merged_line[edge] == -1) {
            return adj_colinear_edges;
          }
          int corners[] = {config.anno.edges[edge].first,
                           config.anno.edges[edge].second};
          for (int corner : corners) {
            for (int adj_edge : point2edges[corner]) {
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
    for (auto &colinear_edges : colinear_edge_groups) {
      if (colinear_edges.size() <= 1) {
        continue;
      }
      // get colinear corners
      std::set<int> colinear_corners;
      for (int edge : colinear_edges) {
        colinear_corners.insert(
            {config.anno.edges[edge].first, config.anno.edges[edge].second});
      }
      if (colinear_corners.size() < 3) {
        continue;
      }
      auto &merged_line =
          merged_lines[edge2merged_line[*colinear_edges.begin()]];
      auto dir1 = cam.direction(merged_line.first);
      auto dir2 = cam.direction(merged_line.second);
      Vec3 line_perp_dir = dir1.cross(dir2);
      basic_constraints.push_back(PlaneConstraint{
          std::vector<int>(colinear_corners.begin(), colinear_corners.end()),
          MakePlaneMatrixAlongDirection(line_perp_dir)});
    }
  }

  // BeamSearch
  Weighted<std::vector<Point3>> best_result;
  best_result.weight() = std::numeric_limits<double>::infinity();

  if (rerun_reconstruction || !config.loadCache("best_result", best_result)) {
    std::default_random_engine rng;
    std::vector<int> vps_table;
    for (auto &vp_edges : vp2edges) {
      vps_table.push_back(vp_edges.first);
    }

    auto beam_search_energy_fun =
        [&basic_constraints, &edge2vp, &config, &faces_overlap, &face_sets,
         &points2edge, &edge2faces, &point2dir, &edge2line, &cam, &vps,
         &vps_table, &best_result, &rng, nfaces, nedges,
         npoints](const std::vector<bool> &vps_in_use_disable_flags) -> double {

      std::set<int> enabled_vps_in_use;
      for (int i = 0; i < vps_table.size(); i++) {
        if (!vps_in_use_disable_flags[i]) {
          enabled_vps_in_use.insert(vps_table[i]);
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
          auto canvas = gui::MakeCanvas(config.anno.image);
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
          constraints.push_back(PlaneConstraint{
              {config.anno.edges[edge].first, config.anno.edges[edge].second},
              MakePlaneMatrixTowardDirection(normal)});
        }
      }

      // analyze the dependency matrices of verts and constraints
      static const int repeat_num = 3;
      std::vector<Weighted<std::vector<Point3>>>
          results_table_under_cur_constraints(repeat_num);

      // energy terms
      auto energy_fun = [&config, &face_sets, &edge2line, &faces_overlap,
                         &points2edge, &edge2faces, nfaces, nedges, npoints](
          const Inferencer &infer, const DenseMatd &vars,
          const std::vector<Vec3> &vert2dir) -> double {
        DenseMatd variables = cv::abs(vars) / cv::norm(vars);
        // edge msda & edge ortho
        auto edge_angles = AnglesBetweenAdjacentEdges(
            vert2dir, config.anno.coplanar_points, variables, infer,
            [&points2edge](int p1, int p2) -> bool {
              return Contains(points2edge, MakeOrderedPair(p1, p2));
            }); ///
        double edge_msda = MeanSquaredDeviationOfContainer(edge_angles);
        double edge_ortho =
            MeanSquaredDeviationOfContainer(edge_angles, {0.0, M_PI_2});

        // face msda & face ortho
        double face_msda = MeanSquaredDeviationOfContainer(
            AnglesBetweenAdjacentFaces(nfaces, edge2faces, variables, infer));
        double face_ortho = MeanSquaredDeviationOfContainer(
            AnglesBetweenAdjacentFaces(nfaces, edge2faces, variables, infer,
                                       faces_overlap),
            {M_PI_2, M_PI, M_PI_2 * 3});

        // edge_msda_each_face_set
        double edge_msda_each_face_set = 0.0;
        {
          std::vector<double> each_face_set(face_sets.size(), 0.0);
          for (int i = 0; i < face_sets.size(); i++) {
            each_face_set[i] =
                MeanSquaredDeviationOfContainer(AnglesBetweenAdjacentEdges(
                    vert2dir, config.anno.coplanar_points, variables, infer,
                    [&points2edge](int p1, int p2) -> bool {
                      return Contains(points2edge, MakeOrderedPair(p1, p2));
                    },
                    [i, &face_sets](int face) -> bool {
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
                AnglesBetweenAdjacentFaces(
                    nfaces, edge2faces, variables, infer, faces_overlap,
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
          std::vector<double> edge_proj_ratios(nedges);
          for (int edge = 0; edge < nedges; edge++) {
            double length2d = edge2line[edge].length();
            auto &edge_corners = config.anno.edges[edge];
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
      size_t disabled_vps_num = vps_table.size() - enabled_vps_in_use.size();
      ParallelRun(
          repeat_num, std::thread::hardware_concurrency() - 1,
          [&constraints, &point2dir, &rng, &energy_fun, disabled_vps_num,
           &results_table_under_cur_constraints](int repeat_id) {
            std::vector<int> fundamental_verts;
            std::vector<Point3> cur_corner2position;
            PerformReconstructionParam param;
            param.max_iters = 1000;
            double cur_energy =
                PerformReconstruction(constraints, point2dir,
                                      std::uniform_int_distribution<int>(
                                          0, point2dir.size() - 1)(rng),
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
        for (int e = 0; e < nedges; e++) {
          auto &p1 = best_result_under_cur_constraints
                         .component[config.anno.edges[e].first];
          auto &p2 = best_result_under_cur_constraints
                         .component[config.anno.edges[e].second];
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
        BeamSearch(vps_table.size(), beam_search_energy_fun, 1);
    config.saveCache("best_result", best_result);
  }

  core::Println("final energy after beam searching = ", best_result.weight());

  if (true) { // show final result
    gui::SceneBuilder sb;
    for (int e = 0; e < nedges; e++) {
      auto &p1 = best_result.component[config.anno.edges[e].first];
      auto &p2 = best_result.component[config.anno.edges[e].second];
      sb.add(Line3(p1, p2));
    }
    sb.show(true, false, gui::RenderOptions()
                             .camera(cam)
                             .renderMode(gui::Lines)
                             .fixUpDirectionInCameraMove(false)
                             .winName("final result"));
  }
}