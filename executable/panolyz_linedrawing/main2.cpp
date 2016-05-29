#include "tools.hpp"

// EnergyWeights
struct EnergyWeights {
  double vert_msda_weight;
  double face_msda_weight;
  double face_angle_weight;
  double all_msda_weight;
};

// ComputeEnergy
template <class VT, class HT, class FT, class FaceHandleIterT,
          class VertHandleToPoint3FunT, class FaceHandleToVec3EquationFunT,
          class VertHandle2BeIgnoredFunT>
std::vector<double> ComputeEnergy(
    const Mesh<VT, HT, FT> &mesh, FaceHandleIterT fhs_begin,
    FaceHandleIterT fhs_end, const EnergyWeights &weights,
    VertHandleToPoint3FunT vh2point3, FaceHandleToVec3EquationFunT fh2vec3_eq,
    VertHandle2BeIgnoredFunT vh2ignored, const Point3 &eye = Point3()) {
  using namespace Eigen;
  std::vector<double> energy_terms;

  std::vector<double> all_angles;
  std::map<VertHandle, std::vector<double>> vh2angles;
  for (FaceHandle fh : MakeRange(fhs_begin, fhs_end)) {
    auto &loophhs = mesh.topo(fh).halfedges;
    std::vector<double> face_angles(loophhs.size(), 0.0);
    for (int k = 0; k < loophhs.size(); k++) {
      int knext = (k + 1) % loophhs.size();
      HalfHandle hh1 = loophhs[k];
      HalfHandle hh2 = loophhs[knext];
      assert(mesh.topo(hh1).to() == mesh.topo(hh2).from());
      VertHandle v1 = mesh.topo(hh1).from();
      VertHandle v2 = mesh.topo(hh1).to();
      VertHandle v3 = mesh.topo(hh2).to();

      if (vh2ignored(v2)) {
        continue;
      }

      auto &p1 = vh2point3(v1);
      auto &p2 = vh2point3(v2);
      auto &p3 = vh2point3(v3);
      double angle = AngleBetweenDirected(p1 - p2, p3 - p2);
      face_angles[k] = angle;
      vh2angles[v2].push_back(angle);
      all_angles.push_back(angle);
    }
    double face_mean_angle =
        std::accumulate(face_angles.begin(), face_angles.end(), 0.0) /
        face_angles.size();
    for (double a : face_angles) {
      energy_terms.push_back((a - face_mean_angle) * weights.face_msda_weight /
                             face_angles.size());
    }
  }

  for (auto &vhangles : vh2angles) {
    VertHandle vh = vhangles.first;
    auto &angles = vhangles.second;
    double vert_mean_angle =
        std::accumulate(angles.begin(), angles.end(), 0.0) / angles.size();
    for (double a : angles) {
      energy_terms.push_back((a - vert_mean_angle) * weights.vert_msda_weight /
                             angles.size());
    }
  }

  double all_mean_angle =
      std::accumulate(all_angles.begin(), all_angles.end(), 0.0) /
      all_angles.size();
  for (double a : all_angles) {
    energy_terms.push_back((a - all_mean_angle) * weights.all_msda_weight /
                           all_angles.size());
  }

  std::set<HalfHandle> valid_hhs;
  //for (auto &h : mesh.halfedges()) {
  //  auto hh = h.topo.hd;
  //  auto fh1 = mesh.topo(hh).face;
  //  auto fh2 = mesh.topo(mesh.topo(hh).opposite).face;
  //  if (std::find(fhs_begin, fhs_end, fh1) != fhs_end &&
  //      std::find(fhs_begin, fhs_end, fh2) != fhs_end) {
  //    valid_hhs.push_back(hh);
  //  }
  //}
  for (FaceHandle fh1 : MakeRange(fhs_begin, fhs_end)) {
    for (HalfHandle hh : mesh.topo(fh1).halfedges) {
      HalfHandle oppohh = mesh.topo(hh).opposite;
      if (Contains(MakeRange(fhs_begin, fhs_end), mesh.topo(oppohh).face)) {
        valid_hhs.insert(hh);
        valid_hhs.insert(oppohh);
      }
    }
  }

  for (auto &hh : valid_hhs) {
    auto fh1 = mesh.topo(hh).face;
    auto planeEq1 = fh2vec3_eq(fh1);
    auto fh2 = mesh.topo(mesh.topo(hh).opposite).face;
    auto planeEq2 = fh2vec3_eq(fh2);
    double cos_angle =
        normalize(Vec3(planeEq1[0], planeEq1[1], planeEq1[2]))
            .dot(normalize(Vec3(planeEq2[0], planeEq2[1], planeEq2[2])));
    assert(!IsInfOrNaN(cos_angle));
    energy_terms.push_back((cos_angle)*weights.face_angle_weight /
                           valid_hhs.size());
  }

  return energy_terms;
}



int main(int argc, char **argv, char **env) {
  gui::Singleton::SetCmdArgs(argc, argv, env);
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");
  misc::Matlab matlab;

  // static const std::string name = "gate";
  // static const std::string name = "towerx";
  // static const std::string name = "tower";
  // static const std::string name = "hex";
  // static const std::string name = "triangle";
  // static const std::string name = "twotriangles";
  static const std::string name = "bridge";
  //static const std::string name = "gundam";

  static const std::string cam_name = "cam1";
  static constexpr bool reset_cam = false;

  std::string obj_file_name = "F:\\LineDrawings\\manifold\\" + name +
                        "\\" + name + ".obj";
  std::string cam_file = "F:\\LineDrawings\\manifold\\" + name +
                        "\\" + name + ".obj." + cam_name + ".cereal";

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
            [&hh_proxy2submesh_id, &ctable](HalfHandle hh) { return gui::Black; },
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
    size_t nvars = dependency.fundamental_vhs.size();
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

    // put vh2matrix into matlab
    matlab << "clear;";
    matlab.setVar("nvars", (double)nvars);
    DenseMatd VM = DenseMatd::zeros(vh2matrix.size(), nvars);
    for (auto &vm : vh2matrix) {
      assert(vm.second.rows == 1 && vm.second.cols == nvars);
      for (int i = 0; i < nvars; i++) {
        VM(vm.first.id, i) = vm.second(0, i);
      }
    }
    matlab.setVar("VM", VM);
    DenseMatd FM = DenseMatd::zeros(fh2matrix.size() * 3, nvars);
    for (auto &fm : fh2matrix) {
      assert(fm.second.rows == 3 && fm.second.cols == nvars);
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < nvars; j++) {
          FM(fm.first.id * 3 + i, j) = fm.second(i, j);
        }
      }
    }
    matlab.setVar("FM", FM);

    // collect orientated vh pairs and put into matlab
    DenseMati orient_vpairs;
    DenseMatd orient_ratios;
    {
      std::vector<std::pair<std::pair<int, int>, double>> vpair_ratios;
      for (auto &h : mesh2d.halfedges()) {
        int edge = hh2d2edge[h.topo.hd];
        if (edge2vp[edge] == -1) {
          continue;
        }
        Vec3 vp = normalize(vp2dir[edge2vp[edge]]);
        VertHandle vh1 = h.topo.from();
        VertHandle vh2 = h.topo.to();
        Vec3 dir1 = normalize(cur_cam.direction(mesh2d.data(vh1)));
        Vec3 dir2 = normalize(cur_cam.direction(mesh2d.data(vh2)));
        Ray3 ray1(dir1, vp);
        Ray3 ray2(Origin(), dir2);
        Point3 interp = DistanceBetweenTwoLines(ray1, ray2).second.second;
        double ratio = norm(interp);

        vpair_ratios.emplace_back(std::make_pair(vh1.id, vh2.id), ratio);
      }
      orient_vpairs = DenseMati::zeros(vpair_ratios.size(), 2);
      orient_ratios = DenseMatd::zeros(vpair_ratios.size(), 1);
      for (int i = 0; i < vpair_ratios.size(); i++) {
        orient_vpairs(i, 0) = vpair_ratios[i].first.first;
        orient_vpairs(i, 1) = vpair_ratios[i].first.second;
        orient_ratios(i, 0) = vpair_ratios[i].second;
      }
    }
    matlab.setVar("orient_vpairs", orient_vpairs);
    matlab.setVar("orient_ratios", orient_ratios);

    // solve vh depths
    // DenseMatd fundamental_vh_depths;

    while (true) {

    }

    //Mesh3 mesh_reconstructed;
    //static const bool no_orientation_estimation = false;
    //if (!no_orientation_estimation) {
    //  mesh_reconstructed =
    //      ReconstructWithOrientations(edge2line, edge2vp, vp2dir, hh2d2edge,
    //                                  cur_cam, mesh_proxy, mesh2d, matlab);
    //} else {
    //  mesh_reconstructed =
    //      Transform(mesh2d, [&cur_cam](const Point2 &p) -> Point3 {
    //        return normalize(cur_cam.direction(p));
    //      });
    //}
    //OptimizeWithoutOrientations(mesh_reconstructed, cur_cam, mesh_proxy,
    //                            submeshes, non_corner_vh_proxies);

    //reconstructions.push_back(std::move(mesh_reconstructed));
  }

  return 0;
}
