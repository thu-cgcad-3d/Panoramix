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
    auto gt_ptr = std::make_unique<Input::GroundTruth>();
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

  //auto input = Input::FromObjFile("hex", "cam1");
  auto input = Input::FromImageAnnotation("house_sketch.jpg");
  // auto input = Configuration::FromImageAnnotation("washington.jpg");

  const size_t npoints = input.anno.points.size();
  const size_t nedges = input.anno.edges.size();
  const size_t nfaces = input.anno.coplanar_points.size();

  bool rerun_preprocess = true;
  bool rerun_line_weights = true;
  bool rerun_reconstruction = true;



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

    // compute faces_overlap : whether two faces sharing a common edge overlap?
    // first for each edge of a face, we select a corner of the face which
    // satisfies:
    // 1. the corner is not on the edge
    // 2. the triangle spanned by the corner and the edge is covered by the face
    // polygon in 2d drawing
    std::vector<std::map<int, int>> face2far_corners_for_edges(nfaces);
    for (int face = 0; face < nfaces; face++) {
      auto &cs = input.anno.coplanar_points[face];
      auto &far_corners_for_edges = face2far_corners_for_edges[face];
      TriangulatePolygon(cs.begin(), cs.end(),
                         [&input](int corner) -> const Point2 & {
                           return input.anno.points.at(corner);
                         },
                         [&far_corners_for_edges, &input,
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
      edge2line[edge].second = input.anno.points[input.anno.edges[edge].second];
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

  for (int cam_param_id = 0; cam_param_id < cam_params.size(); cam_param_id++) {
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
      if (false) {
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
    // STEP 4: collect intersections and build the lines voting panel
    std::vector<Vec3> intersections;
    std::vector<std::pair<int, int>> intersection_merged_line_pairs;
    auto vote_map = MakeView<float>(PanoramicCamera(32));
    ElementBinaryRelations<float> merged_lines_parallelism_weights(
        merged_lines.size(), 0.0f);
    if (rerun_line_weights ||
        !input.loadCache(decorated_name_with_cam_param_id("line_weights"),
                         intersections, intersection_merged_line_pairs,
                         vote_map, merged_lines_parallelism_weights)) {
      intersections = CollectLineIntersections(merged_line3s,
                                               &intersection_merged_line_pairs);
      assert(intersections.size() == intersection_merged_line_pairs.size());
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

      // find directions peaky in votes
      Imagef suppressed_votes;
      std::vector<Pixel> peaky_pixels;
      NonMaximaSuppression(vote_map.image, suppressed_votes, 5, &peaky_pixels);

      int orthogonal_pairs_num = 0;
      MakeRange(peaky_pixels)
          .forEachTwo([&vote_map, &orthogonal_pairs_num](const Pixel &p1,
                                                         const Pixel &p2) {
            double angle = AngleBetweenUndirected(
                vote_map.camera.direction(p1), vote_map.camera.direction(p2));
            if (abs(angle - M_PI_2) < DegreesToRadians(5)) {
              orthogonal_pairs_num++;
            }
          });
			int all_pairs_num = peaky_pixels.size() * (peaky_pixels.size() - 1) / 2;
			core::Println("orthogonal pairs ratio = ", orthogonal_pairs_num / double(all_pairs_num));

      for (int i = 0; i < intersections.size(); i++) {
        merged_lines_parallelism_weights(
            intersection_merged_line_pairs[i].first,
            intersection_merged_line_pairs[i].second) =
            vote_map.at(intersections[i]);
      }
      input.saveCache(decorated_name_with_cam_param_id("line_weights"),
                      intersections, intersection_merged_line_pairs, vote_map,
                      merged_lines_parallelism_weights);

      /*  auto d1 = vote_map.camera.direction(Point2(47, 54));
        auto d2 = vote_map.camera.direction(Point2(58, 38));

        double angle12 = AngleBetweenUndirected(d1, d2);

        auto d3 = vote_map.camera.direction(Point2(148, 44));
        auto d4 = vote_map.camera.direction(Point2(158, 61));

        double angle34 = AngleBetweenUndirected(d3, d4);

        double angle13 = AngleBetweenUndirected(d1, d3);
        double angle24 = AngleBetweenUndirected(d2, d4);
        double angle14 = AngleBetweenUndirected(d1, d4);
        double angle23 = AngleBetweenUndirected(d2, d3);*/

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

      // enhance orthogonal peaks in vote_map

      /*Imagef new_votes(vote_map.size(), 0.0f);
      static const double kernal_angle = DegreesToRadians(10.0);
      Image_<Vec3> dirs_map(vote_map.size());
      for (auto it = dirs_map.begin(); it != dirs_map.end(); ++it) {
        *it = cam.direction(it.pos());
      }
      const double radius = kernal_angle * vote_map.camera.focal() / 2.0;
      for (int y = 0; y < new_votes.rows; y++) {
        int y2_lb = std::max(y - radius, 0.0);
        int y2_ub = std::min<int>(y + radius, new_votes.rows - 1);
        double pixel_weight =
            std::max(cos((y2_lb - new_votes.rows / 2.0) / new_votes.rows *
      M_PI),
                     cos((y2_ub - new_votes.rows / 2.0) / new_votes.rows *
      M_PI));
        double xspan = radius / pixel_weight;
        if (xspan > new_votes.cols) {
          xspan = new_votes.cols;
        }
        for (int x = 0; x < new_votes.cols; x++) {
          int count = 0;
          for (int y2 = y2_lb; y2 <= y2_ub; y2++) {
            for (int x2 = x - xspan / 2; x2 <= x + xspan / 2; x2++) {
              Pixel p2(WrapBetween(x2, 0, new_votes.cols), y2);
              double angle = AngleBetweenUndirected(dirs_map(y, x),
      dirs_map(p2));
              if (angle >= kernal_angle / 5.0) {
                continue;
              }
              count ++;
              double v =
                  (Gaussian(angle, kernal_angle / 10.0)) * (vote_map.image(p2));
              if (v > new_votes(y, x)) {
                new_votes(y, x) = v;
              }
            }
          }
        }
      }
      vote_map.image = new_votes;*/
     
    }

    ////////////////////////////////////////////////
    // STEP 5: collect constraints
    // basic_constraints
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
        basic_constraints.push_back(PlaneConstraint{
            std::vector<int>(colinear_corners.begin(), colinear_corners.end()),
            MakePlaneMatrixAlongDirection(line_perp_dir)});
      }
    }

    ////////////////////////////////////////////////
    // STEP 6: build the energy function
    auto energy_fun = [&](const Inferencer &infer,
                          const DenseMatd &vars) -> double {
      DenseMatd variables = cv::abs(vars) / cv::norm(vars);

    };

	}
}