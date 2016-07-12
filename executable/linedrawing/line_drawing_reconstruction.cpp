#include "line_drawing_reconstruction.hpp"
#include "line_drawing_evaluation.hpp"
#include "line_drawing_tools.hpp"

namespace pano {
namespace experimental {

// FromObjFile
LineDrawingInput LineDrawingInput::FromObjFile(const std::string &model_name,
                                               const std::string &cam_name) {
  std::string folder =
      PANORAMIX_TEST_DATA_DIR_STR "\\linedrawing\\objs\\" + model_name + "\\";

  // load 3d line drawing from file
  LineDrawingInput input;
  input.model_name = model_name;
  input.cam_name = cam_name;

  input.groundtruth = std::make_shared<LineDrawingGroundTruth>();
  input.groundtruth->line_drawing =
      LineDrawing3FromObjFile(folder + model_name + "_w_intf.obj");

  // setup a camera
  std::string cam_path = folder + model_name + ".obj." + cam_name + ".cereal";
  if (!LoadFromDisk(cam_path, input.groundtruth->camera)) {
    gui::SceneBuilder sb;
    for (auto &f : input.groundtruth->line_drawing.topo.coplanar_points) {
      assert(f.size() >= 3);
      for (int i = 0; i < f.size(); i++) {
        sb.add(Line3(
            input.groundtruth->line_drawing.points[f[i]],
            input.groundtruth->line_drawing.points[f[(i + 1) % f.size()]]));
      }
    }
    input.groundtruth->camera =
        sb.show(true, true, gui::RenderOptions()
                                .renderMode(gui::Lines)
                                .fixUpDirectionInCameraMove(false))
            .camera();
    SaveToDisk(cam_path, input.groundtruth->camera);
  }
  core::Println("gt focal = ", input.groundtruth->camera.focal(), " gt pp = ",
                input.groundtruth->camera.principlePoint());

  // build projection data
  // points
  input.projection.line_drawing.points.resize(
      input.groundtruth->line_drawing.points.size());
  for (int i = 0; i < input.projection.line_drawing.points.size(); i++) {
    input.projection.line_drawing.points[i] =
        input.groundtruth->camera.toScreen(
            input.groundtruth->line_drawing.points[i]);
  }
  // topo
  input.projection.line_drawing.topo = input.groundtruth->line_drawing.topo;

  // image (white board)
  input.projection.image =
      Image3ub(input.groundtruth->camera.screenSize(), Vec3ub(255, 255, 255));

  return input;
}

// FromImageAnnotation
LineDrawingInput
LineDrawingInput::FromImageAnnotation(const std::string &image_name) {
  static const std::string folder =
      PANORAMIX_TEST_DATA_DIR_STR "\\linedrawing\\images\\";
  std::string anno_path = folder + image_name + ".ldanno.cereal";
  LineDrawingImageProjection anno;
  bool succ = LoadFromDisk(anno_path, anno);
  if (!succ) {
    core::Println("Failed to load ", image_name, "!");
    return LineDrawingInput();
  }
  return LineDrawingInput{"image_" + image_name, "no_cam", std::move(anno),
                          nullptr};
}

// RunLineDrawingReconstruction
void RunLineDrawingReconstruction(const LineDrawingInput &input,
                                  LineDrawingReconstructionSteps &steps,
                                  GUILevel gui_level) {

  core::Println("#### input.model_name = [", input.model_name,
                "], input.cam_name = [", input.cam_name, "]");
  std::string id =
      "linedrawing[" + input.model_name + "][" + input.cam_name + "]";

  const size_t npoints = input.projection.line_drawing.points.size();
  const size_t nedges = input.projection.line_drawing.topo.edges.size();
  const size_t nfaces =
      input.projection.line_drawing.topo.coplanar_points.size();

  ////////////////////////////////////////////////
  // STEP 1: compute additional data
  AuxiliaryData<LineDrawing2> aux_data =
      MakeAuxiliary(input.projection.line_drawing);
  const auto &points2edge = aux_data.points2edge;
  const auto &point2edges = aux_data.point2edges;
  const auto &face2edges = aux_data.face2edges;
  const auto &edge2faces = aux_data.edge2faces;
  const auto &edge2line = aux_data.edge2line;
  Box2 lines_box = BoundingBoxOfContainer(edge2line);

  auto faces_overlap =
      ComputeFacesOverlap(input.projection.line_drawing, aux_data);

  ////////////////////////////////////////////////
  // STEP 2: decompose faces and estimate camera parameters
  std::vector<std::set<int>> face_sets;
  std::vector<CameraParam> cam_params;
  steps.preprocess.perform(id, face_sets, cam_params)([&]() {
    face_sets =
        DecomposeFaces(input.projection.line_drawing.topo.coplanar_points,
                       input.projection.line_drawing.points);

    // calibrate camera
    cam_params = CalibrateCamera(
        BoundingBoxOfContainer(input.projection.line_drawing.points), face_sets,
        [&input, &points2edge](int face) {
          std::vector<std::vector<int>> pids = {{}};
          auto &vs = input.projection.line_drawing.topo.coplanar_points[face];
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
              chain.append(input.projection.line_drawing.points[p]);
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
  });

  std::default_random_engine rng;

  // the best result
  Weighted<std::vector<Point3>> best_result;
  best_result.weight() = std::numeric_limits<double>::infinity();
  int best_cam_param_id = -1;

  for (int cam_param_id = 0; cam_param_id < cam_params.size(); cam_param_id++) {
    core::Println("cam param id = ", cam_param_id);
    auto decorated_name_with_cam_param_id =
        [cam_param_id](const std::string &name) {
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
      cam = PerspectiveCamera(input.projection.image.cols,
                              input.projection.image.rows, cur_cam_param.pp,
                              cur_cam_param.focal);

      // calculate point directions
      point2dir = MakeRange(input.projection.line_drawing.points)
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
      if (gui_level >= GUILevelAll) {
        auto canvas = gui::MakeCanvas(input.projection.image);
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
          canvas.add(NoteAs(input.projection.line_drawing.points[i],
                            std::to_string(i)));
        }
        // show noted edges
        canvas.color(gui::Blue);
        canvas.paintingOptions().fontScale = 0.9;
        for (int i = 0; i < nedges; i++) {
          canvas.add(NoteAs(
              (input.projection.line_drawing
                   .points[input.projection.line_drawing.topo.edges[i].first] +
               input.projection.line_drawing.points
                   [input.projection.line_drawing.topo.edges[i].second]) /
                  2.0,
              "[" + std::to_string(i) + "]"));
        }
        canvas.show(0, "merged lines");
      }
    }

    ////////////////////////////////////////////////
    // STEP 4: find potential vps
    std::vector<Weighted<Vec3>> vps;
    steps.find_vps.perform(decorated_name_with_cam_param_id(id), vps)([&]() {
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
      NonMaximaSuppression(vote_map.image, suppressed_votes, 5, &peaky_pixels);

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

      //// filter
      // double thres = std::max_element(vps.begin(), vps.end())->score / 2.0;
      // vps = MakeRange(vps)
      //          .filter([thres](auto &vp) { return vp.score > thres; })
      //          .evalAsStdVector();
      // core::Println("filtered vps num = ", vps.size());

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
      /*std::sort(vps.begin(), vps.end(), std::greater<>());
      vps.resize(std::min<int>(10, vps.size()));*/

      core::Println("num of final vps with orthogonal relations = ",
                    vps.size());

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
    });

    ////////////////////////////////////////////////
    // STEP 5: orient some lines!
    std::vector<int> merged_line2vp(merged_lines.size(), -1);
    std::vector<int> edge2vp(nedges, -1);
    std::set<int> vps_with_edges;
    steps.estimate_orientations.perform(decorated_name_with_cam_param_id(id),
                                        merged_line2vp, edge2vp,
                                        vps_with_edges)([&]() {
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
      if (gui_level >= GUILevelAll) {
        for (auto &vp_edges : vp2edges) {
          auto canvas = gui::MakeCanvas(input.projection.image);
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
    });

    ////////////////////////////////////////////////
    // STEP 6: collect basic constraints
    std::vector<PlaneConstraint> basic_constraints;
    {
      // add face coplanarities
      for (int face = 0; face < nfaces; face++) {
        basic_constraints.push_back(PlaneConstraint{
            input.projection.line_drawing.topo.coplanar_points[face],
            MakePlaneMatrix()});
      }
      // add colinearity from assignment
      for (auto &ps : input.projection.line_drawing.topo.colinear_points) {
        assert(ps.size() >= 3);
        auto dir1 =
            cam.direction(input.projection.line_drawing.points[ps.front()]);
        auto dir2 =
            cam.direction(input.projection.line_drawing.points[ps.back()]);
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
            int corners[] = {
                input.projection.line_drawing.topo.edges[edge].first,
                input.projection.line_drawing.topo.edges[edge].second};
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
              {input.projection.line_drawing.topo.edges[edge].first,
               input.projection.line_drawing.topo.edges[edge].second});
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
    // STEP 7: build the energy function
    auto energy_fun = [&](const Inferencer &infer, const DenseMatd &vars,
                          const std::vector<Vec3> &vert2dir) -> double {
      DenseMatd variables = cv::abs(vars) / cv::norm(vars);
      // edge msda & edge ortho
      auto edge_angles = AnglesBetweenAdjacentEdges(
          vert2dir, input.projection.line_drawing.topo.coplanar_points,
          variables, infer, [&points2edge](int p1, int p2) -> bool {
            return Contains(points2edge, MakeOrderedPair(p1, p2));
          }); //// the angles are in [0, M_PI]
      double edge_msda = MeanSquaredDeviationOfContainer(edge_angles);
      double edge_ortho =
          MeanSquaredDeviationOfContainer(edge_angles, {M_PI, M_PI_2}); ////

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
                  vert2dir, input.projection.line_drawing.topo.coplanar_points,
                  variables, infer,
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
              AnglesBetweenAdjacentFaces(nfaces, edge2faces, variables, infer,
                                         faces_overlap,
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
          auto &edge_corners = input.projection.line_drawing.topo.edges[edge];
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
      /* {
         e = 100.0 * edge_msda + 0 * face_msda + 0 * edge_ortho +
             100 * face_ortho + 100.0 * edge_msda_each_face_set +
             100 * face_ortho_each_face_set +
             std::max(max_edge_skew_ratio_score - 2.5, 0.0) * 100;
       }*/
      {
        e = 100.0 * edge_msda + 0 * face_msda + 30 * edge_ortho +
            100 * face_ortho + 100.0 * edge_msda_each_face_set +
            100 * face_ortho_each_face_set +
            std::max(max_edge_skew_ratio_score - 2.5, 0.0) * 100;
      }
      return e;
    };

    ////////////////////////////////////////////////
    // STEP 8: beam search
    // BeamSearch
    steps.reconstruction.perform(decorated_name_with_cam_param_id(id),
                                 best_result, best_cam_param_id)([&] {

      std::vector<int> vp_ids(vps_with_edges.begin(), vps_with_edges.end());

      auto beam_search_fun =
          [&](const std::vector<bool> &vps_in_use_disable_flags) -> double {
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
                {input.projection.line_drawing.topo.edges[edge].first,
                 input.projection.line_drawing.topo.edges[edge].second},
                MakePlaneMatrixTowardDirection(normal)});
          }
        }

        // perform reconstruction in this node
        static const int repeat_num = 3;
        std::vector<Weighted<std::vector<Point3>>> results_in_this_node(
            repeat_num);
        size_t disabled_vps_num = vp_ids.size() - vps_used_in_this_node.size();
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

        if (gui_level >= GUILevelNecessary) { // show current node result
          gui::SceneBuilder sb;
          for (int e = 0; e < nedges; e++) {
            auto &p1 = best_result_in_this_node.component
                           [input.projection.line_drawing.topo.edges[e].first];
            auto &p2 = best_result_in_this_node.component
                           [input.projection.line_drawing.topo.edges[e].second];
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
    });

    core::Println("best energy till cam param id = ", cam_param_id, " is ",
                  best_result.weight());
    if (gui_level >= GUILevelAll) { // best result current
      gui::SceneBuilder sb;
      for (int e = 0; e < nedges; e++) {
        auto &p1 =
            best_result
                .component[input.projection.line_drawing.topo.edges[e].first];
        auto &p2 =
            best_result
                .component[input.projection.line_drawing.topo.edges[e].second];
        sb.add(Line3(p1, p2));
      }
      sb.show(true, false,
              gui::RenderOptions()
                  .camera(PerspectiveCamera(
                      input.projection.image.cols, input.projection.image.rows,
                      cam_params[best_cam_param_id].pp,
                      cam_params[best_cam_param_id].focal))
                  .renderMode(gui::Lines)
                  .fixUpDirectionInCameraMove(false)
                  .winName("best result till cam param id = " +
                           std::to_string(cam_param_id)));
    }
  }

  auto best_cam = PerspectiveCamera(
      input.projection.image.cols, input.projection.image.rows,
      cam_params[best_cam_param_id].pp, cam_params[best_cam_param_id].focal);

  if (gui_level >= GUILevelNecessary) { // show final result
    gui::SceneBuilder sb;
    for (int e = 0; e < nedges; e++) {
      auto &p1 =
          best_result
              .component[input.projection.line_drawing.topo.edges[e].first];
      auto &p2 =
          best_result
              .component[input.projection.line_drawing.topo.edges[e].second];
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
}