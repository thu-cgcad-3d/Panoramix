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

struct Input {
  std::string model_name, cam_name;
  LineDrawingAnnotation anno;

  struct GroundTruth {
    std::vector<Point3> points;
    PerspectiveCamera camera;
  };
  std::shared_ptr<GroundTruth> gt;

  template <class... Ts>
  bool loadCache(const std::string &what, Ts &... ts) const {
    return misc::LoadCache(model_name + "_" + cam_name, what, ts...);
  }
  template <class... Ts>
  bool saveCache(const std::string &what, Ts &... ts) const {
    return misc::SaveCache(model_name + "_" + cam_name, what, ts...);
  }

  // FromObjFile
  static Input FromObjFile(const std::string &model_name,
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
    auto gt_ptr = std::make_shared<Input::GroundTruth>();
    gt_ptr->camera = std::move(cam);
    gt_ptr->points = std::move(point3ds);

    return Input{model_name, cam_name, std::move(anno), std::move(gt_ptr)};
  }

  // FromImageAnnotation
  static Input FromImageAnnotation(const std::string &image_name) {
    static const std::string folder =
        PANORAMIX_TEST_DATA_DIR_STR "\\linedrawing\\images\\";
    std::string anno_path = folder + image_name + ".ldanno.cereal";
    LineDrawingAnnotation anno;
    bool succ = LoadFromDisk(anno_path, anno);
    if (!succ) {
      core::Println("Failed to load ", image_name, "!");
      return Input();
    }
    return Input{"image_" + image_name, "no_cam", std::move(anno), nullptr};
  }
};

void RoutineReconstruct1() {

	int gui_level = 0;
  bool rerun_preprocess = false;
  bool rerun_find_vps = false;
  bool rerun_estimate_orientations = false;
  bool rerun_reconstruction = false;

  std::vector<Input> inputs = {Input::FromObjFile("plane", "cam1"),
                               Input::FromObjFile("towers", "cam1"),
                               Input::FromObjFile("tower", "cam1"),
                               Input::FromObjFile("towerx", "cam1"),
                               Input::FromObjFile("car", "cam1"),
                               Input::FromObjFile("gate", "cam1"),
                               Input::FromObjFile("smallgate", "cam1"),
                               Input::FromObjFile("castle", "cam1"),
                               Input::FromObjFile("desk", "cam1"),
                               Input::FromImageAnnotation("house_sketch.jpg")};

  for (auto &input : inputs) {
    core::Println("#### input.model_name = [", input.model_name,
                  "], input.cam_name = [", input.cam_name, "]");
    const size_t npoints = input.anno.points.size();
    const size_t nedges = input.anno.edges.size();
    const size_t nfaces = input.anno.coplanar_points.size();

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
        auto &e = input.anno.edges[edge];
        points2edge[e] = edge;
        point2edges[e.first].insert(edge);
        point2edges[e.second].insert(edge);
      }
      for (int face = 0; face < nfaces; face++) {
        auto &ps = input.anno.coplanar_points[face];
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
        auto &cs = input.anno.coplanar_points[face];
        auto &far_corners_for_edges = face2far_corners_for_edges[face];
        TriangulatePolygon(
            cs.begin(), cs.end(),
            [&input](int corner) -> const Point2 & {
              return input.anno.points.at(corner);
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
      // the we compute face overlapnesss
      for (int edge = 0; edge < nedges; edge++) {
        Line2 line(input.anno.points[input.anno.edges[edge].first],
                   input.anno.points[input.anno.edges[edge].second]);
        auto &fs = edge2faces.at(edge);
        if (fs.size() < 2) {
          continue;
        }
        for (auto iti = fs.begin(); iti != fs.end(); ++iti) {
          int far_corner_i = face2far_corners_for_edges.at(*iti).at(edge);
          bool on_left_side_of_edge_i = IsOnLeftSide(
              input.anno.points.at(far_corner_i), line.first, line.second);
          for (auto itj = std::next(iti); itj != fs.end(); ++itj) {
            int far_corner_j = face2far_corners_for_edges.at(*itj).at(edge);
            bool on_left_side_of_edge_j = IsOnLeftSide(
                input.anno.points.at(far_corner_j), line.first, line.second);
            faces_overlap[MakeOrderedPair(*iti, *itj)] =
                on_left_side_of_edge_i == on_left_side_of_edge_j;
          }
        }
      }

      // construct 2d lines for edges
      for (int edge = 0; edge < nedges; edge++) {
        edge2line[edge].first = input.anno.points[input.anno.edges[edge].first];
        edge2line[edge].second =
            input.anno.points[input.anno.edges[edge].second];
      }
      lines_box = BoundingBoxOfContainer(edge2line);
    }

    ////////////////////////////////////////////////
    // STEP 2: decompose into face sets and estimate camera parameters
    std::vector<std::set<int>> face_sets;
    std::vector<CameraParam> cam_params;
    if (rerun_preprocess ||
        !input.loadCache("preprocess", face_sets, cam_params)) {
      // decompose the drawing
      face_sets = DecomposeFaces(input.anno.coplanar_points, input.anno.points);
      // calibrate camera
      cam_params = CalibrateCamera(
          BoundingBoxOfContainer(input.anno.points), face_sets,
          [&input, &points2edge](int face) {
            std::vector<std::vector<int>> pids = {{}};
            auto &vs = input.anno.coplanar_points[face];
            for (int v : vs) {
              if (pids.back().empty() ||
                  Contains(points2edge,
                           MakeOrderedPair(pids.back().back(), v))) {
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
                chain.append(input.anno.points[p]);
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
                assert(face_loops.size() == 1);
              }
            }
            return face_loops;
          },
          20);
      input.saveCache("preprocess", face_sets, cam_params);
    }

    std::default_random_engine rng;

    // the best result
    Weighted<std::vector<Point3>> best_result;
    best_result.weight() = std::numeric_limits<double>::infinity();
    int best_cam_param_id = -1;

    for (int cam_param_id = 0; cam_param_id < cam_params.size();
         cam_param_id++) {
      core::Println("cam param id = ", cam_param_id);
      auto decorated_name_with_cam_param_id = [cam_param_id](
          const std::string &name) {
        return name + "_cam_param_" + std::to_string(cam_param_id);
      };

      ////////////////////////////////////////////////
      // STEP 3: recover camera and merge lines
      PerspectiveCamera cam;
      std::vector<Vec3> point2dir;
      std::vector<int> edge2merged_line;
      std::vector<Line2> merged_lines;
      std::vector<Line3> edge2line3;
      std::vector<Line3> merged_line3s;
      {
        assert(!cam_params.empty());
        auto &cur_cam_param = cam_params[cam_param_id];
        core::Println("current focal = ", cur_cam_param.focal, " pp = ",
                      cur_cam_param.pp);
        cam = PerspectiveCamera(input.anno.image.cols, input.anno.image.rows,
                                cur_cam_param.pp, cur_cam_param.focal);

        // calculate point directions
        point2dir = MakeRange(input.anno.points)
                        .transform([&cam](const Point2 &p) -> Vec3 {
                          return normalize(cam.direction(p));
                        })
                        .evalAsStdVector();

        // merge lines
        merged_lines = MergeColinearLines(
            edge2line, cur_cam_param, DegreesToRadians(0.1), &edge2merged_line);
        edge2line3 = MakeRange(edge2line)
                         .transform([&cam](const Line2 &line2) -> Line3 {
                           return Line3(normalize(cam.direction(line2.first)),
                                        normalize(cam.direction(line2.second)));
                         })
                         .evalAsStdVector();
        merged_line3s =
            MakeRange(merged_lines)
                .transform([&cam](const Line2 &line2) -> Line3 {
                  return Line3(normalize(cam.direction(line2.first)),
                               normalize(cam.direction(line2.second)));
                })
                .evalAsStdVector();
        // show merged lines
        if (gui_level >= 2) {
          auto canvas = gui::MakeCanvas(input.anno.image);
          canvas.color(gui::LightGray);
          canvas.colorTable(
              gui::CreateRandomColorTableWithSize(merged_lines.size()));
          canvas.thickness(2);
          for (int i = 0; i < merged_lines.size(); i++) {
            canvas.add(ClassifyAs(merged_lines[i], i));
          }
          // show noted corners
          canvas.color(gui::Black);
          canvas.paintingOptions().fontScale = 1.0;
          for (int i = 0; i < npoints; i++) {
            canvas.add(NoteAs(input.anno.points[i], std::to_string(i)));
          }
          // show noted edges
          canvas.color(gui::Blue);
          canvas.paintingOptions().fontScale = 0.9;
          for (int i = 0; i < nedges; i++) {
            canvas.add(NoteAs((input.anno.points[input.anno.edges[i].first] +
                               input.anno.points[input.anno.edges[i].second]) /
                                  2.0,
                              "[" + std::to_string(i) + "]"));
          }
          canvas.show(0, "merged lines");
        }
      }

      ////////////////////////////////////////////////
      // STEP 4: find potential vps
      std::vector<Weighted<Vec3>> vps;
      if (rerun_find_vps ||
          !input.loadCache(decorated_name_with_cam_param_id("vps"), vps)) {
        std::vector<std::pair<int, int>> intersection_merged_line_pairs;
        std::vector<Vec3> intersections = CollectLineIntersections(
            merged_line3s, &intersection_merged_line_pairs);
        assert(intersections.size() == intersection_merged_line_pairs.size());
        auto vote_map = MakeView<float>(PanoramicCamera(128));
        for (auto &inter : intersections) {
          auto p1 = ToPixel(vote_map.camera.toScreen(inter));
          auto p2 = ToPixel(vote_map.camera.toScreen(-inter));
          assert(p1 != p2);
          auto d1 = vote_map.camera.direction(p1);
          auto d2 = vote_map.camera.direction(p2);
          assert(Distance(normalize(d1), -normalize(d2)) < 0.1);
          vote_map(inter) += 1.0f;
          vote_map(-inter) += 1.0f;
        }

        // find strong directions peaky in votes
        Imagef suppressed_votes;
        std::vector<Pixel> peaky_pixels;
        NonMaximaSuppression(vote_map.image, suppressed_votes, 5,
                             &peaky_pixels);

        core::Println("all peaky pixel num = ", peaky_pixels.size());
        vps.reserve(peaky_pixels.size());
        for (auto &pixel : peaky_pixels) {
          if (pixel.y > vote_map.image.rows / 2) {
            continue;
          }
          Vec3 dir = vote_map.camera.direction(pixel);
          // remove pixels that are too close to existing points
          bool too_close_to_existing_point = false;
          for (auto &pdir : point2dir) {
            if (AngleBetweenUndirected(pdir, dir) < DegreesToRadians(5)) {
              too_close_to_existing_point = true;
              break;
            }
          }
          if (too_close_to_existing_point) {
            continue;
          }
          double votes = 0.0;
          ForEachPixelWithinViewCone(vote_map, dir, DegreesToRadians(2),
                                     [&votes, &vote_map](const Pixel &p) {
                                       votes += vote_map.image(p);
                                       vote_map.image(p) = 0.0f;
                                     });
          ForEachPixelWithinViewCone(vote_map, -dir, DegreesToRadians(2),
                                     [&votes, &vote_map](const Pixel &p) {
                                       votes += vote_map.image(p);
                                       vote_map.image(p) = 0.0f;
                                     });
          if (votes >= 3) {
            vps.push_back(
                WeightAs(normalize(vote_map.camera.direction(pixel)), votes));
          }
        }
        core::Println("initial vps num = ", vps.size());

        std::vector<int> initial_vps_ortho_count(vps.size(), 0);
        for (int i = 0; i < vps.size(); i++) {
          for (int j = i + 1; j < vps.size(); j++) {
            double angle =
                AngleBetweenUndirected(vps[i].component, vps[j].component);
            if (abs(angle - M_PI_2) < DegreesToRadians(5)) {
              initial_vps_ortho_count[i]++;
              initial_vps_ortho_count[j]++;
            }
          }
        }
        auto filtered_vps =
            MakeIotaRange<int>(vps.size())
                .filter([&initial_vps_ortho_count](int i) -> bool {
                  return initial_vps_ortho_count[i] > 0;
                })
                .transform([&vps](int i) { return vps[i]; })
                .evalAsStdVector();
        vps = std::move(filtered_vps);

        core::Println("num of final vps with orthogonal relations = ",
                      vps.size());

        input.saveCache(decorated_name_with_cam_param_id("vps"), vps);

        //{
        //  auto canvas = gui::MakeCanvas(Image3ub(vote_map.size()));
        //  for (auto &line3 : merged_line3s) {
        //    double angle = AngleBetweenDirected(line3.first, line3.second);
        //    double step_angle = DegreesToRadians(0.1);
        //    for (double a = angle; a <= M_PI * 2; a += step_angle) {
        //      auto d1 = RotateDirection(line3.first, line3.second, a -
        //      step_angle);
        //      auto d2 = RotateDirection(line3.first, line3.second, a);
        //      auto p1 = vote_map.camera.toScreen(d1);
        //      auto p2 = vote_map.camera.toScreen(d2);
        //      canvas.color(gui::Red);
        //      canvas.add(Line2(p1, p2));
        //    }
        //    for (double a = step_angle; a <= angle; a += step_angle) {
        //      auto d1 = RotateDirection(line3.first, line3.second, a -
        //      step_angle);
        //      auto d2 = RotateDirection(line3.first, line3.second, a);
        //      auto p1 = vote_map.camera.toScreen(d1);
        //      auto p2 = vote_map.camera.toScreen(d2);
        //	canvas.color(gui::Black);
        //      canvas.add(Line2(p1, p2));
        //    }
        //  }
        //  canvas.show();
        //}
      }

      ////////////////////////////////////////////////
      // STEP 5: orient some lines!
      std::vector<int> merged_line2vp(merged_lines.size(), -1);
      std::vector<int> edge2vp(nedges, -1);
      std::set<int> vps_with_edges;
      if (rerun_estimate_orientations ||
          !input.loadCache(decorated_name_with_cam_param_id("orientations"),
                           merged_line2vp, edge2vp, vps_with_edges)) {
        // build relation of adjacent merged lines
        std::vector<std::pair<int, int>> adjacent_merged_line_pairs;
        for (int p = 0; p < npoints; p++) {
          auto &edges = point2edges[p];
          for (auto i = edges.begin(); i != edges.end(); ++i) {
            int mlinei = edge2merged_line[*i];
            for (auto j = std::next(i); j != edges.end(); ++j) {
              int mlinej = edge2merged_line[*j];
              if (mlinei != mlinej) {
                adjacent_merged_line_pairs.emplace_back(mlinei, mlinej);
              }
            }
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

        EstimateEdgeOrientationsParam param;
        param.solve_max_iter = 10;
        param.angle_thres_allowed_vp_line_deviation = DegreesToRadians(5);
        merged_line2vp = EstimateEdgeOrientations(
            merged_line3s,
            MakeRange(vps)
                .transform([](auto &&scored) { return scored.component; })
                .evalAsStdVector(),
            adjacent_merged_line_pairs, face2merged_lines, param);

        std::map<int, std::set<int>> vp2edges;
        for (int edge = 0; edge < nedges; edge++) {
          int vp = merged_line2vp[edge2merged_line[edge]];
          edge2vp[edge] = vp;
          if (vp != -1) {
            vp2edges[vp].insert(edge);
            vps_with_edges.insert(vp);
          }
        }

        // show line orientation estimation
        if (gui_level >= 2) {
          for (auto &vp_edges : vp2edges) {
            auto canvas = gui::MakeCanvas(input.anno.image);
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
            canvas.show(-1, "vp_" + std::to_string(vp_edges.first));
          }
          cv::waitKey();
        }

        input.saveCache(decorated_name_with_cam_param_id("orientations"),
                        merged_line2vp, edge2vp, vps_with_edges);
      }

      ////////////////////////////////////////////////
      // STEP 6: collect basic constraints
      std::vector<PlaneConstraint> basic_constraints;
      {
        // add face coplanarities
        for (int face = 0; face < nfaces; face++) {
          basic_constraints.push_back(PlaneConstraint{
              input.anno.coplanar_points[face], MakePlaneMatrix()});
        }
        // add colinearity from assignment
        for (auto &ps : input.anno.colinear_points) {
          assert(ps.size() >= 3);
          auto dir1 = cam.direction(input.anno.points[ps.front()]);
          auto dir2 = cam.direction(input.anno.points[ps.back()]);
          Vec3 line_perp_dir = dir1.cross(dir2);
          assert(
              (!HasValue(line_perp_dir, [](auto e) { return IsInfOrNaN(e); })));
          basic_constraints.push_back(
              PlaneConstraint{std::vector<int>(ps.begin(), ps.end()),
                              MakePlaneMatrixAlongDirection(line_perp_dir)});
        }
        // add colinearity for edges that share same merged lines
        std::vector<std::set<int>> colinear_edge_groups(nedges);
        colinear_edge_groups.resize(ConnectedComponents(
            MakeIotaIterator<int>(0), MakeIotaIterator<int>(nedges),
            [&input, &edge2merged_line, &point2edges](int edge) {
              std::set<int> adj_colinear_edges;
              if (edge2merged_line[edge] == -1) {
                return adj_colinear_edges;
              }
              int corners[] = {input.anno.edges[edge].first,
                               input.anno.edges[edge].second};
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
                {input.anno.edges[edge].first, input.anno.edges[edge].second});
          }
          if (colinear_corners.size() < 3) {
            continue;
          }
          auto &merged_line =
              merged_lines[edge2merged_line[*colinear_edges.begin()]];
          auto dir1 = cam.direction(merged_line.first);
          auto dir2 = cam.direction(merged_line.second);
          Vec3 line_perp_dir = dir1.cross(dir2);
          basic_constraints.push_back(
              PlaneConstraint{std::vector<int>(colinear_corners.begin(),
                                               colinear_corners.end()),
                              MakePlaneMatrixAlongDirection(line_perp_dir)});
        }
      }

      ////////////////////////////////////////////////
      // STEP 7: build the energy function
      auto energy_fun = [&](const Inferencer &infer, const DenseMatd &vars,
                            const std::vector<Vec3> &vert2dir) -> double {
        DenseMatd variables = cv::abs(vars) / cv::norm(vars);
        // edge msda & edge ortho
        auto edge_angles = AnglesBetweenAdjacentEdges(
            vert2dir, input.anno.coplanar_points, variables, infer,
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
                    vert2dir, input.anno.coplanar_points, variables, infer,
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
            auto &edge_corners = input.anno.edges[edge];
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
        {
          e = 100.0 * edge_msda + 0 * face_msda + 0 * edge_ortho +
              100 * face_ortho + 100.0 * edge_msda_each_face_set +
              100 * face_ortho_each_face_set +
              std::max(max_edge_skew_ratio_score - 2.5, 0.0) * 100;
        }
        return e;
      };

      ////////////////////////////////////////////////
      // STEP 8: beam search
      // BeamSearch
      if (rerun_reconstruction ||
          !input.loadCache(decorated_name_with_cam_param_id("best_result"),
                           best_result, best_cam_param_id)) {

        std::vector<int> vp_ids(vps_with_edges.begin(), vps_with_edges.end());

        auto beam_search_fun = [&](
            const std::vector<bool> &vps_in_use_disable_flags) -> double {
          std::set<int> vps_used_in_this_node;
          for (int i = 0; i < vp_ids.size(); i++) {
            if (!vps_in_use_disable_flags[i]) {
              vps_used_in_this_node.insert(vp_ids[i]);
            }
          }

          std::string msg = "searching with enabled_vps_in_use = {";
          for (int vp : vps_used_in_this_node) {
            msg += std::to_string(vp) + " ";
          }
          msg += "}";
          misc::Clock clock(msg);

          // add more constraints provided by oriented edges
          std::vector<PlaneConstraint> constraints = basic_constraints;
          for (int edge = 0; edge < edge2vp.size(); edge++) {
            if (Contains(vps_used_in_this_node, edge2vp[edge])) { //
              auto &line = edge2line[edge];
              const Vec3 &vp = vps[edge2vp[edge]].component;
              Vec3 line_perp_dir =
                  cam.direction(line.first).cross(cam.direction(line.second));
              Vec3 normal = normalize(vp.cross(line_perp_dir));
              constraints.push_back(PlaneConstraint{
                  {input.anno.edges[edge].first, input.anno.edges[edge].second},
                  MakePlaneMatrixTowardDirection(normal)});
            }
          }

          // perform reconstruction in this node
          static const int repeat_num = 3;
          std::vector<Weighted<std::vector<Point3>>> results_in_this_node(
              repeat_num);
          size_t disabled_vps_num =
              vp_ids.size() - vps_used_in_this_node.size();
          ParallelRun(repeat_num, std::thread::hardware_concurrency() - 1,
                      [&constraints, &point2dir, &rng, &energy_fun,
                       disabled_vps_num, &results_in_this_node](int repeat_id) {
                        std::vector<int> fundamental_verts;
                        std::vector<Point3> cur_corner2position;
                        PerformReconstructionParam param;
                        param.max_iters = 500;
                        double cur_energy =
                            PerformReconstruction(
                                constraints, point2dir,
                                std::uniform_int_distribution<int>(
                                    0, point2dir.size() - 1)(rng),
                                energy_fun, rng, cur_corner2position,
                                &fundamental_verts, param) +
                            disabled_vps_num * 0.1; // compensation
                        // fill cur_energy and cur_corner2position to the tables
                        results_in_this_node[repeat_id].component =
                            std::move(cur_corner2position);
                        results_in_this_node[repeat_id].weight() = cur_energy;
                      });

          // get the best result in this node
          auto &best_result_in_this_node = *std::min_element(
              results_in_this_node.begin(), results_in_this_node.end());
          core::Println("best energy in this node = ",
                        best_result_in_this_node.weight());

          if (gui_level >= 1) { // show current node result
            gui::SceneBuilder sb;
            for (int e = 0; e < nedges; e++) {
              auto &p1 =
                  best_result_in_this_node.component[input.anno.edges[e].first];
              auto &p2 = best_result_in_this_node
                             .component[input.anno.edges[e].second];
              sb.add(Line3(p1, p2));
            }
            sb.show(true, true, gui::RenderOptions()
                                    .renderMode(gui::Lines)
                                    .fixUpDirectionInCameraMove(false)
                                    .winName("best result in this node"));
          }
          if (best_result_in_this_node < best_result) {
            best_result = best_result_in_this_node;
            best_cam_param_id = cam_param_id;
          }
          return best_result_in_this_node.weight();
        };

        auto best_vps_in_use_disable_flags =
            BeamSearch(vp_ids.size(), beam_search_fun, 1);

        input.saveCache(decorated_name_with_cam_param_id("best_result"),
                        best_result, best_cam_param_id);
      }
    }

    auto best_cam = PerspectiveCamera(
        input.anno.image.cols, input.anno.image.rows,
        cam_params[best_cam_param_id].pp, cam_params[best_cam_param_id].focal);
    
		input.saveCache("final_result", best_result, best_cam);

    if (gui_level >= 1) { // show final result
      gui::SceneBuilder sb;
      for (int e = 0; e < nedges; e++) {
        auto &p1 = best_result.component[input.anno.edges[e].first];
        auto &p2 = best_result.component[input.anno.edges[e].second];
        sb.add(Line3(p1, p2));
      }
      sb.show(true, false, gui::RenderOptions()
                               .camera(best_cam)
                               .renderMode(gui::Lines)
                               .fixUpDirectionInCameraMove(false)
                               .winName("final result"));
    }
  }
}