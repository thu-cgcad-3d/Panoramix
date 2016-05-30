#include "tools.hpp"

int main(int argc, char **argv, char **env) {
  gui::Singleton::SetCmdArgs(argc, argv, env);
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");
  misc::Matlab matlab;

  // static const std::string name = "gate";
   
  //static const std::string name = "towerx";
  //static const std::string cam_name = "cam2";

  // static const std::string name = "tower";
  // static const std::string name = "hex";
  // static const std::string name = "triangle";
  // static const std::string name = "twotriangles";
  
  static const std::string name = "bridge";
  static const std::string cam_name = "cam1";

  // static const std::string name = "gundam";

  static constexpr bool reset_cam = false;

  std::string obj_file_name =
      "F:\\LineDrawings\\manifold\\" + name + "\\" + name + ".obj";
  std::string cam_file = "F:\\LineDrawings\\manifold\\" + name + "\\" + name +
                         ".obj." + cam_name + ".cereal";

  //// [Load Mesh]
  auto mesh = LoadMeshFromObjFile(obj_file_name);




  if (true) { // show original mesh
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XLines;
    sb.installingOptions().lineWidth = 3;
    for (auto &h : mesh.halfedges()) {
      auto line = Line3(mesh.data(h.topo.from()), mesh.data(h.topo.to()));
      if (h.topo.face.valid() && mesh.topo(h.topo.opposite).face.valid()) {
        sb.installingOptions().discretizeOptions.color(gui::Red);
        sb.add(line);
      } else {
        sb.installingOptions().discretizeOptions.color(gui::Black);
        sb.add(line);
      }
    }

    sb.show(true, true, gui::RenderOptions()
                            .backgroundColor(gui::White)
                            .renderMode(gui::All)
                            .bwTexColor(0.0)
                            .bwColor(1.0)
                            .fixUpDirectionInCameraMove(false)
                            .cullBackFace(false)
                            .cullFrontFace(false));
  }

  auto mesh_proxy = MakeMeshProxy(mesh);

  //// [Decompose]
  auto cut_face_pairs =
      DecomposeAll(mesh_proxy, [](HalfHandle hh1, HalfHandle hh2) -> bool {
        // TODO: fix this
        return false;
      });
  {
    // check validity
    for (auto &h : mesh_proxy.halfedges()) {
      HalfHandle hh = h.data;
      if (hh.invalid()) {
        hh = mesh_proxy.data(h.topo.opposite);
      }
      assert(hh.valid());
    }
  }

  auto submeshes = ExtractSubMeshes(mesh_proxy);
  Println("found ", submeshes.size(), " submeshes");

  //// [Load Camera]
  PerspectiveCamera cam;
  if (!LoadFromDisk(cam_file, cam) || reset_cam) {
    auto sphere = BoundingBoxOfContainer(mesh.vertices()).outerSphere();
    PerspectiveCamera projCam(500, 500, Point2(250, 250), 200,
                              sphere.center + Vec3(1, 2, 3) * sphere.radius * 2,
                              sphere.center);

    const gui::ColorTable ctable =
        gui::CreateRandomColorTableWithSize(submeshes.size());

    HandledTable<HalfHandle, int> hh_proxy2submesh_id(
        mesh_proxy.internalHalfEdges().size(), -1);
    for (auto &h : mesh_proxy.halfedges()) {
      int &id = hh_proxy2submesh_id[h.topo.hd];
      for (int i = 0; i < submeshes.size(); i++) {
        if (Contains(submeshes[i], h.topo.hd)) {
          id = i;
          break;
        }
      }
    }
    HandledTable<FaceHandle, int> fh_proxy2submesh_id(
        mesh_proxy.internalFaces().size(), -1);
    for (auto &f : mesh_proxy.faces()) {
      int &id = fh_proxy2submesh_id[f.topo.hd];
      for (int i = 0; i < submeshes.size(); i++) {
        if (Contains(submeshes[i], f.topo.hd)) {
          id = i;
          break;
        }
      }
    }

    // show each subMesh
    if (true) {
      for (int i = 0; i < submeshes.size(); i++) {
        Println("subMesh - ", i);
        gui::SceneBuilder sb;
        sb.installingOptions().defaultShaderSource =
            gui::OpenGLShaderSourceDescriptor::XTriangles;
        sb.installingOptions().discretizeOptions.color(gui::Black);
        sb.installingOptions().lineWidth = 0.03;
        AddToScene(
            sb, mesh_proxy,
            [&submeshes, i](auto h) { return Contains(submeshes[i], h); },
            [&mesh](VertHandle vh) { return mesh.data(vh); },
            [&hh_proxy2submesh_id, &ctable](HalfHandle hh) {
              return gui::Black;
            },
            [&fh_proxy2submesh_id, &ctable](FaceHandle fh) {
              return ctable[fh_proxy2submesh_id[fh]];
            });
        sb.show(true, false, gui::RenderOptions()
                                 .camera(projCam)
                                 .backgroundColor(gui::White)
                                 .renderMode(gui::All)
                                 .bwTexColor(0.0)
                                 .bwColor(1.0)
                                 .fixUpDirectionInCameraMove(false)
                                 .cullBackFace(false)
                                 .cullFrontFace(false));
      }
    }

    { // show together
      gui::SceneBuilder sb;
      sb.installingOptions().defaultShaderSource =
          gui::OpenGLShaderSourceDescriptor::XTriangles;
      sb.installingOptions().discretizeOptions.color(gui::Black);
      sb.installingOptions().lineWidth = 0.03;
      AddToScene(
          sb, mesh_proxy, [](auto) { return true; },
          [&mesh](VertHandle vh) { return mesh.data(vh); },
          [&hh_proxy2submesh_id, &ctable](HalfHandle hh) { return gui::Black; },
          [&fh_proxy2submesh_id, &ctable](FaceHandle fh) {
            return ctable[fh_proxy2submesh_id[fh]];
          });
      cam = sb.show(true, false, gui::RenderOptions()
                                     .camera(projCam)
                                     .backgroundColor(gui::White)
                                     .renderMode(gui::All)
                                     .bwTexColor(0.0)
                                     .bwColor(1.0)
                                     .fixUpDirectionInCameraMove(false)
                                     .cullBackFace(false)
                                     .cullFrontFace(false))
                .camera();
    }
    SaveToDisk(cam_file, cam);
  }

  Println("gt camera: focal = ", cam.focal(), " pp = ", cam.principlePoint());

  //// [Make 2D Mesh]
  // convert to 2d
  auto mesh2d = Transform(
      mesh, [&cam](const Point3 &p) -> Point2 { return cam.toScreen(p); });

  // add offset noise
  Vec2 offset_noise = Vec2(20, -20);
  for (auto &v : mesh2d.vertices()) {
    v.data += offset_noise;
  }

  if (true) {
    Image3ub im(cam.screenSize(), Vec3ub(255, 255, 255));
    auto canvas = gui::MakeCanvas(im);
    canvas.color(gui::Black);
    canvas.thickness(2);
    for (auto &h : mesh2d.halfedges()) {
      auto &p1 = mesh2d.data(h.topo.from());
      auto &p2 = mesh2d.data(h.topo.to());
      canvas.add(Line2(p1, p2));
    }
    canvas.show(0, "mesh2d");
  }

  // find all vh proxies lying on an edge
  std::unordered_set<VertHandle> non_corner_vh_proxies;
  for (auto &v : mesh_proxy.vertices()) {
    VertHandle vh = v.data;
    const Point2 &phere = mesh2d.data(vh);
    auto &hh_proxies = v.topo.halfedges;

    bool parallel_edge_found = false;
    for (int i = 0; !parallel_edge_found && i < hh_proxies.size(); i++) {
      VertHandle another_vhi =
          mesh_proxy.data(mesh_proxy.topo(hh_proxies[i]).to());
      if (another_vhi == vh) {
        another_vhi = mesh_proxy.data(mesh_proxy.topo(hh_proxies[i]).from());
      }
      assert(another_vhi != vh);
      const Point2 &ptherei = mesh2d.data(another_vhi);

      for (int j = i + 1; !parallel_edge_found && j < hh_proxies.size(); j++) {
        VertHandle another_vhj =
            mesh_proxy.data(mesh_proxy.topo(hh_proxies[j]).to());
        if (another_vhj == vh) {
          another_vhj = mesh_proxy.data(mesh_proxy.topo(hh_proxies[j]).from());
        }
        assert(another_vhj != vh);
        const Point2 &ptherej = mesh2d.data(another_vhj);

        if (IsFuzzyZero(AngleBetweenDirected(ptherei - phere, ptherej - phere) -
                            M_PI,
                        1e-8)) {
          non_corner_vh_proxies.insert(v.topo.hd);
          parallel_edge_found = true;
          break;
        }
      }
    }
  }

  auto pp_focal_groups =
      EstimateCameraParameters(mesh2d, mesh_proxy, submeshes);

  //// estimate vanishing points
  HandledTable<HalfHandle, int> hh2d2edge;
  std::vector<std::map<int, double>> vp2edge_with_angles;
  std::vector<std::vector<Scored<int>>> edge2ordered_vp_and_angles;
  std::vector<Line2> edge2line;
  auto vp_positions = EstimateVanishingPoints(
      mesh2d, cam.screenSize(), &hh2d2edge, &vp2edge_with_angles,
      &edge2ordered_vp_and_angles, &edge2line);

  int nvps = vp_positions.size();
  int nedges = edge2ordered_vp_and_angles.size();

  // collect edge intersections and
  // get vp_positions from the intersections
  std::vector<int> edge2vp;
  std::vector<std::vector<int>> vp2edges;

  // factor graph
  if (false ||
      !misc::LoadCache(obj_file_name + "-" + cam_name, "edge2vp_vp2edges",
                       edge2vp, vp2edges)) {
    EstimateEdgeOrientations(mesh2d, vp_positions, hh2d2edge,
                             vp2edge_with_angles, edge2ordered_vp_and_angles,
                             edge2vp, vp2edges);
    misc::SaveCache(obj_file_name + "-" + cam_name, "edge2vp_vp2edges", edge2vp,
                    vp2edges);
  }

  if (false) { // show line classification results
    for (int i = 0; i < nvps; i++) {
      if (vp2edges[i].empty()) {
        continue;
      }
      Image3ub im(cam.screenSize(), Vec3ub(255, 255, 255));
      auto canvas = gui::MakeCanvas(im);
      canvas.color(gui::LightGray);
      canvas.thickness(2);
      for (auto &line : edge2line) {
        canvas.add(line);
      }
      canvas.color(gui::Gray);
      canvas.thickness(2);
      for (int edge : vp2edges[i]) {
        canvas.add(edge2line[edge].ray());
      }
      canvas.color(gui::Black);
      for (int edge : vp2edges[i]) {
        canvas.add(edge2line[edge]);
      }
      canvas.show(0, "optimized vp_" + std::to_string(i));
    }
  }

  assert(edge2vp.size() == nedges);
  assert(vp2edges.size() == nvps);

  // analyze the fh/vh dependencies in each subMesh
  FaceVertDependency dependency;
  {
    std::vector<FaceHandle> fhs;
    fhs.reserve(mesh2d.internalFaces().size());
    for (auto &f : mesh2d.faces()) {
      fhs.push_back(f.topo.hd);
    }
    auto fundamental_vhs_set = SortFacesWithPlanarityDependency(
        mesh2d, fhs.begin(), fhs.end(),
        [&mesh2d](auto vhs_begin, auto vhs_end) -> bool {
          auto trans_fun = [&mesh2d](VertHandle vh) -> const Point2 & {
            return mesh2d.data(vh);
          };
          return IsFuzzyColinear(MakeTransformIterator(vhs_begin, trans_fun),
                                 MakeTransformIterator(vhs_end, trans_fun),
                                 1e-6);
        });
    dependency.ordered_fhs = std::move(fhs);
    dependency.fundamental_vhs = std::vector<VertHandle>(
        fundamental_vhs_set.begin(), fundamental_vhs_set.end());
  }

  matlab << "clear;";
  // prepare matlab data
  size_t nvars = dependency.fundamental_vhs.size();
  matlab.setVar("nvars", (double)nvars);

  // collect adjacent fh pairs and put into matlab
  DenseMati adj_fpairs;
  {
    std::set<std::pair<int, int>> fpairs;
    for (auto &h : mesh2d.halfedges()) {
      FaceHandle fh1 = h.topo.face;
      FaceHandle fh2 = mesh2d.topo(h.topo.opposite).face;
      fpairs.insert(std::make_pair(fh1.id, fh2.id));
    }
    adj_fpairs = DenseMati::zeros(fpairs.size(), 2);
    int row = 0;
    for (auto &p : fpairs) {
      adj_fpairs(row, 0) = p.first;
      adj_fpairs(row, 1) = p.second;
      row++;
    }
  }
  matlab.setVar("adj_fpairs", adj_fpairs);
  matlab << "adj_fpairs = adj_fpairs + 1;";

  //// collect vh triplets forming an angle in a face
  // DenseMati adj_vtriplets;
  //{
  // std::vector<std::tuple<int, int, int>> vtriplets;

  //}
  // matlab.setVar("adj_vtriplets", adj_vtriplets);
  // matlab << "adj_vtriplets = adj_vtriplets + 1;";

  // for each pp focal candidate
  std::vector<Mesh3> reconstructions;
  for (int config_id = 0;
       config_id < std::min(5ull, pp_focal_groups.size()) &&
       pp_focal_groups[config_id].first.size() * 30 >= pp_focal_groups.size();
       config_id++) {

    double focal = pp_focal_groups[config_id].second.focal;
    auto &pp = pp_focal_groups[config_id].second.pp;
    Println("current camera: focal = ", focal, " pp = ", pp);

    PerspectiveCamera cur_cam(cam.screenWidth(), cam.screenHeight(), pp, focal);
    std::vector<Vec3> vp2dir(nvps);
    for (int i = 0; i < nvps; i++) {
      vp2dir[i] = cur_cam.direction(vp_positions[i]);
    }

    // collect directions and matrices
    HandledTable<VertHandle, Vec3> vh2dir = mesh2d.createVertexTable(Vec3());

    // construct matrices
    std::map<VertHandle, DenseMatd> vh2matrix;
    std::map<FaceHandle, DenseMatd> fh2matrix;
    {
      ConstructSingleViewTransformMatrices(
          dependency.fundamental_vhs,
          dependency.ordered_fhs, // face2related_verts
          [&mesh2d, &cur_cam](FaceHandle fh) -> std::vector<VertHandle> {
            std::vector<VertHandle> related_vhs;
            for (HalfHandle hh : mesh2d.topo(fh).halfedges) {
              related_vhs.push_back(mesh2d.topo(hh).from());
            }
            return related_vhs;
          }, // vert2direction
          [&mesh2d, &cur_cam](VertHandle vh) -> Vec3 {
            return normalize(cur_cam.direction(mesh2d.data(vh)));
          },
          vh2matrix, fh2matrix);
    }

    // set VM FM using vh2matrix and fh2matrix
    DenseMatd VM, FM;
    {
      VM = DenseMatd::zeros(vh2matrix.size(), nvars);
      for (auto &vm : vh2matrix) {
        assert(vm.second.rows == 1 && vm.second.cols == nvars);
        for (int i = 0; i < nvars; i++) {
          VM(vm.first.id, i) = vm.second(0, i);
        }
      }
      matlab.setVar("VM", VM);
      FM = DenseMatd::zeros(fh2matrix.size() * 3, nvars);
      for (auto &fm : fh2matrix) {
        assert(fm.second.rows == 3 && fm.second.cols == nvars);
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < nvars; j++) {
            FM(fm.first.id * 3 + i, j) = fm.second(i, j);
          }
        }
      }
      matlab.setVar("FM", FM);
    }

    // define the reconstruction functor
    auto reconstruct_fun = [&](const DenseMatd &X) -> Mesh3 {
      DenseMatd inversed_depths = VM * X;
      Mesh3 reconstructed =
          Transform(mesh2d, [&cur_cam, &mesh2d](const Point2 &p2) {
            return normalize(cur_cam.direction(p2));
          });
      for (auto &v : reconstructed.vertices()) {
        v.data /= inversed_depths(v.topo.hd.id, 0);
      }
      return reconstructed;
    };

    // define the final energy functor
    auto final_energy_fun = [&](const DenseMatd &X) -> DenseMatd {
      auto reconstructed = reconstruct_fun(X);

    };

    // collect orientated vh pairs and put into matlab
    DenseMati orient_vpairs;
    DenseMatd orient_ratios; // second depth / first depth,
                             // or first inv-depth / second inv-depth
    {
      std::map<std::pair<int, int>, double> vpair_ratios;
      for (auto &h : mesh2d.halfedges()) {
        int edge = hh2d2edge[h.topo.hd];
        if (edge2vp[edge] == -1) {
          continue;
        }
        Vec3 vp = normalize(vp2dir[edge2vp[edge]]);
        VertHandle vh1 = h.topo.from();
        VertHandle vh2 = h.topo.to();
        auto vpair = std::make_pair(vh1.id, vh2.id);
        if (Contains(vpair_ratios, vpair) ||
            Contains(vpair_ratios, ReversePair(vpair))) {
          continue;
        }

        Vec3 dir1 = normalize(cur_cam.direction(mesh2d.data(vh1)));
        Vec3 dir2 = normalize(cur_cam.direction(mesh2d.data(vh2)));
        Ray3 ray1(dir1, vp);
        Ray3 ray2(Origin(), dir2);
        Point3 interp = DistanceBetweenTwoLines(ray1, ray2).second.second;
        double ratio = norm(interp);

        vpair_ratios[vpair] = ratio;
      }
      orient_vpairs = DenseMati::zeros(vpair_ratios.size(), 2);
      orient_ratios = DenseMatd::zeros(vpair_ratios.size(), 1);
      int row = 0;
      for (auto &p : vpair_ratios) {
        orient_vpairs(row, 0) = p.first.first;
        orient_vpairs(row, 1) = p.first.second;
        orient_ratios(row, 0) = p.second;
        row++;
      }
    }
    matlab.setVar("orient_vpairs", orient_vpairs);
    matlab.setVar("orient_ratios", orient_ratios);
    matlab << "orient_vpairs = orient_vpairs + 1;";

    // solve vh depths
    // simulated anealing first
    DenseMatd X = DenseMatd::ones(nvars, 1);

    // SimulatedAnnealing(X, )

    matlab << "X_prev = ones(nvars, 1);";

    while (false) {
      matlab << "cvx_begin";
      matlab << "variable X(nvars);";

      auto construct_energy = [&matlab](const std::string &objective_name) {
        // vert inv depths
        matlab << "vert_inv_depths = VM * X;";
        matlab << "vert_inv_depths_first = "
                  "vert_inv_depths(orient_vpairs(:, 1));";
        matlab << "vert_inv_depths_second = "
                  "vert_inv_depths(orient_vpairs(:, 2));";
        matlab << "diff1 = vert_inv_depths_first - vert_inv_depths_second .* "
                  "orient_ratios;";
        matlab << "diff2 = vert_inv_depths_first ./ orient_ratios - "
                  "vert_inv_depths_second;";

        // face equations
        matlab << "face_plane_eqs = reshape(FM * X, [3, size(FM, 1) / 3])';";
        matlab << "face_plane_eqs_first = face_plane_eqs(adj_fpairs(:, 1), "
                  ":);";
        matlab << "face_plane_eqs_prev = reshape(FM * X_prev, [3, size(FM, "
                  "1) / 3])';";
        matlab << "face_plane_eqs_norm_prev = sqrt(sum(face_plane_eqs_prev.^2, "
                  "2));";
        matlab << "face_plane_eqs_second_prev = "
                  "face_plane_eqs_prev(adj_fpairs(:, 2), :);";
        matlab << "diff3 = dot(face_plane_eqs_first, "
                  "face_plane_eqs_second_prev, 2) ./ "
                  "face_plane_eqs_norm_prev(adj_fpairs(:, 1)) ./ "
                  "face_plane_eqs_norm_prev(adj_fpairs(:, 2));";

        matlab << (objective_name + " = sum_square(diff1) + "
                                    "sum_square(diff2) + "
                                    "sum_square(diff3) * 0;");
      };

      construct_energy("energy");
      matlab << "minimize energy * 1e5";
      matlab << "subject to";
      matlab << "	ones(size(vert_inv_depths)) <= vert_inv_depths;";
      matlab << "cvx_end";
      matlab << "X = X ./ norm(X);";

      X = matlab.var("X").toCVMat();
      assert(X.rows == nvars && X.cols == 1);
      construct_energy("energy");
      Println("energy = ", matlab.var("energy").scalar());

      { // show current X
        gui::SceneBuilder sb;
        sb.installingOptions().defaultShaderSource =
            gui::OpenGLShaderSourceDescriptor::XLines;
        sb.add(reconstruct_fun(X));
        sb.show(true, true, gui::RenderOptions()
                                .renderMode(gui::Lines)
                                .fixUpDirectionInCameraMove(false));
      }
      matlab << "X_prev = X;";
    }

    while (true) {
      matlab << "clear;";
      Mesh3 mesh_reconstructed =
          ReconstructWithOrientations(edge2line, edge2vp, vp2dir, hh2d2edge,
                                      cur_cam, mesh_proxy, mesh2d, matlab, 1e3);
      OptimizeWithoutOrientations(mesh_reconstructed, cur_cam, mesh_proxy,
                                  submeshes, non_corner_vh_proxies);

      {
        // udpate edge2vp
        for (auto &h : mesh2d.halfedges()) {
          int edge = hh2d2edge[h.topo.hd];
          int vp = edge2vp[edge];
          if (vp == -1) {
            continue;
          }
          const Vec3 &vpdir = vp2dir[vp];
          VertHandle vh1 = h.topo.from();
          VertHandle vh2 = h.topo.to();
          Line3 line(mesh_reconstructed.data(vh1),
                     mesh_reconstructed.data(vh2));
          // check whether the reconstruced line satisfy the vp direction
          double angle = AngleBetweenUndirected(line.direction(), vpdir);
          //Println("angle = ", angle);
          if (angle > DegreesToRadians(10)) {
            edge2vp[edge] = -1;
          }
        }
      }
      if (true) {
        gui::SceneBuilder sb;
        sb.installingOptions().defaultShaderSource =
            gui::OpenGLShaderSourceDescriptor::XLines;
        sb.installingOptions().discretizeOptions.color(gui::Black);
        sb.installingOptions().lineWidth = 10;
        sb.add(mesh_reconstructed);
        sb.show(true, true, gui::RenderOptions()
                                .camera(cur_cam)
                                .renderMode(gui::Lines)
                                .winName("After Holistic Optimization")
                                .fixUpDirectionInCameraMove(false));
      }
    }
  }

  return 0;
}
