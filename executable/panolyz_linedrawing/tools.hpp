#pragma once

#include "cache.hpp"
#include "clock.hpp"
#include "factor_graph.hpp"
#include "parallel.hpp"

#include "canvas.hpp"
#include "gui_util.hpp"
#include "qttools.hpp"
#include "scene.hpp"
#include "singleton.hpp"

#include "line_drawing.hpp"
#include "mesh_advanced_util.hpp"
#include "pi_graph_annotation.hpp"
#include "pi_graph_cg.hpp"
#include "pi_graph_control.hpp"
#include "pi_graph_occlusion.hpp"
#include "pi_graph_optimize.hpp"
#include "pi_graph_vis.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

// FaceVertDependency
struct FaceVertDependency {
  std::vector<FaceHandle> ordered_fhs;
  std::vector<VertHandle> fundamental_vhs;
};

// PPFocalCandidate
struct PPFocalCandidate {
  Point2 pp;
  double focal;
};

// VHFSelectedFunT: (VertHandle/HalfHandle/FaceHandle) -> bool
// VertPositionFunT: (VertHandle)->Point3
// HalfColorerFunT: (HalfHandle)->gui::Color
// FaceColorerFunT: (FaceHandle)->gui::Color
template <class VertDataT, class HalfDataT, class FaceDataT,
          class VHFSelectedFunT, class VertPositionFunT, class HalfColorerFunT,
          class FaceColorerFunT>
void AddToScene(gui::SceneBuilder &sb,
                const Mesh<VertDataT, HalfDataT, FaceDataT> &m,
                VHFSelectedFunT selected, VertPositionFunT vertPosFun,
                HalfColorerFunT colorHalf, FaceColorerFunT colorFace);

// PossibleKeyVanishingPoints
std::vector<Point2> PossibleKeyVanishingPoints(const Chain2 &chain);

// EstimateCameraParameters
std::vector<std::pair<std::set<int>, PPFocalCandidate>>
EstimateCameraParameters(
    const Mesh2 &mesh2d,
    const Mesh<VertHandle, HalfHandle, FaceHandle> &mesh_proxy,
    const std::vector<SubMesh> &submeshes);

// EstimateVanishingPoints
std::vector<Point2> EstimateVanishingPoints(
    const Mesh2 &mesh2d, const Sizei &drawing_size,
    HandledTable<HalfHandle, int> *hh2d2edge_ptr = nullptr,
    std::vector<std::map<int, double>> *vp2edge_with_angles_ptr = nullptr,
    std::vector<std::vector<Scored<int>>> *edge2ordered_vp_and_angles_ptr =
        nullptr,
    std::vector<Line2> *edge2line_ptr = nullptr);

// EstimateEdgeOrientations
void EstimateEdgeOrientations(
    const Mesh2 &mesh2d, const std::vector<Point2> &vp_positions,
    const HandledTable<HalfHandle, int> &hh2d2edge,
    const std::vector<std::map<int, double>> &vp2edge_with_angles,
    const std::vector<std::vector<Scored<int>>> &edge2ordered_vp_and_angles,
    std::vector<int> &edge2vp, std::vector<std::vector<int>> &vp2edges);

// ReconstructWithOrientations
Mesh3 ReconstructWithOrientations(
    const std::vector<Line2> &edge2line, const std::vector<int> &edge2vp,
    const std::vector<Vec3> &vp2dir,
    const HandledTable<HalfHandle, int> &hh2d2edge,
    const PerspectiveCamera &cur_cam,
    const Mesh<VertHandle, HalfHandle, FaceHandle> &mesh_proxy,
    const Mesh2 &mesh2d, misc::Matlab &matlab,
    double energy_ratio_con_over_ortho = 1e3);

// OptimizeWithoutOrientations
void OptimizeWithoutOrientations(
    Mesh3 &cur_reconstruction, const PerspectiveCamera &cur_cam,
    const Mesh<VertHandle, HalfHandle, FaceHandle> &mesh_proxy,
    const std::vector<SubMesh> &submeshes,
    const std::unordered_set<VertHandle> &non_corner_vh_proxies);



// VHFSelectedFunT: (VertHandle/HalfHandle/FaceHandle) -> bool
// VertPositionFunT: (VertHandle)->Point3
// HalfColorerFunT: (HalfHandle)->gui::Color
// FaceColorerFunT: (FaceHandle)->gui::Color
template <class VertDataT, class HalfDataT, class FaceDataT,
          class VHFSelectedFunT, class VertPositionFunT, class HalfColorerFunT,
          class FaceColorerFunT>
void AddToScene(gui::SceneBuilder &sb,
                const Mesh<VertDataT, HalfDataT, FaceDataT> &m,
                VHFSelectedFunT selected, VertPositionFunT vertPosFun,
                HalfColorerFunT colorHalf, FaceColorerFunT colorFace) {
  sb.installingOptions().lineWidth = 1.0;
  sb.installingOptions().defaultShaderSource =
      gui::OpenGLShaderSourceDescriptor::XLines;
  HandledTable<HalfHandle, int> added(m.internalHalfEdges().size(), false);
  for (auto &h : m.halfedges()) {
    if (!selected(h.topo.hd) || !selected(h.topo.from()) ||
        !selected(h.topo.to())) {
      continue;
    }
    Line3 line(vertPosFun(m.data(h.topo.from())),
               vertPosFun(m.data(h.topo.to())));
    auto hh = h.topo.hd;
    auto fh = h.topo.face;
    auto oppohh = h.topo.opposite;
    auto oppofh = oppohh.valid() ? m.topo(oppohh).face : FaceHandle();
    bool hasFace = fh.valid();
    bool hasOppo = oppohh.valid();
    gui::Color color = colorHalf(hh);
    if (color.isTransparent()) {
      continue;
    }
    if (!hasFace && hasOppo) {
      color = gui::Red;
    } else if (hasFace && !hasOppo) {
      color = gui::Blue;
    } else if (!hasFace && !hasOppo) {
      color = gui::Yellow;
    }
    if (added[oppohh]) {
      continue;
    }
    sb.add(gui::ColorAs(line, color), [hh, fh, oppohh, oppofh](auto &...) {
      std::cout << "halfedge id: " << hh.id
                << ", opposite halfedge id: " << oppohh.id
                << ", face id: " << fh.id << ", opposite face id: " << oppofh.id
                << '\n';
    });
    added[hh] = true;
  }
  sb.installingOptions().defaultShaderSource =
      gui::OpenGLShaderSourceDescriptor::XTriangles;
  for (auto &f : m.faces()) {
    if (!selected(f.topo.hd)) {
      continue;
    }
    Polygon3 poly;
    for (auto h : f.topo.halfedges) {
      auto v = m.topo(h).to();
      poly.corners.push_back(vertPosFun(m.data(v)));
    }
    assert(poly.corners.size() > 2);
    poly.normal = (poly.corners[0] - poly.corners[1])
                      .cross(poly.corners[0] - poly.corners[2]);
    auto fh = f.topo.hd;
    gui::Color color = colorFace(f.topo.hd);
    if (color.isTransparent()) {
      continue;
    }
    sb.add(gui::ColorAs(poly, color),
           [fh](auto &...) { std::cout << "face id: " << fh.id << '\n'; });
  }
}


// PerspectiveReconstruction
template <class VertT, class FaceT, class Face2RelatedVertsFunT>
std::map<VertT, double> PerspectiveReconstruction(
	const std::map<VertT, Vec3> & v2dir,
    const std::map<VertT, DenseMatd> &v2matrix,
    const std::map<FaceT, DenseMatd> &f2matrix, const std::vector<Vec3> &vps,
    const std::map<std::pair<VertT, VertT>, int> &vpair2vp,
    misc::Matlab &matlab) {

  assert(!v2matrix.empty());
  int nvars = v2matrix->begin()->second.cols;
  assert(nvars > 0);
  for (auto &vm : v2matrix) {
    assert(vm.second.cols == nvars);
  }
  for (auto &fm : f2matrix) {
    assert(fm.second.cols == nvars);
  }

  int nverts = v2matrix.size();
  int nfaces = f2matrix.size();

  std::map<VertT, int> vert2position;
  DenseMatd verts_mat(v2matrix.size(), 3);
  for (auto & vm : v2matrix) {
	verts_mat.row(vert2position.size()).setTo(vm.second);
    vert2position[vm.first] = vert2position.size();
  }
  std::map<FaceT, int> face2positon;
  DenseMatd faces_mat(f2matrix.size() * 3, 3);
  for (auto & fm : f2matrix) {
    faces_mat
        .rowRange(
            cv::Range(face2positon.size() * 3, face2positon.size() * 3 + 3))
        .setTo(fm.second);
    face2positon[fm.first] = face2positon.size();
  }

  matlab.setVar("nverts", nverts);
  matlab.setVar("nfaces", nfaces);
  matlab.setVar("verts_mat", verts_mat);
  matlab.setVar("faces_mat", faces_mat);

  int nvpairs = vpair2vp.size();
  DenseMati vpair_lefts(npairs), vpair_rights(npairs);
  DenseMati vpair_vpids(npairs);
  DenseMatd vpair_ratios(npairs);
  int vpair_i = 0;
  for (auto &vpair_vpid : vpair2vp) {
    vpair_lefts[vpair_i] = vert2position.at(vpair_vpid.first.first);
    vpair_rights[vpair_i] = vert2position.at(vpair_vpid.first.second);
    vpair_vpids[vpair_i] = vpair_vpid.second;

	// compute ratio
	Vec3 vp = normalize(vps[vpair_vpid.second]);
	Vec3 dir1 = normalize(v2dir.at(vpair_vpid.first.first));
	Vec3 dir2 = normalize(v2dir.at(vpair_vpid.first.second));
	Ray3 ray1(dir1, vp);
	Ray3 ray2(Origin(), dir2);
	Point3 interp = Intersection(ray1, ray2);
	double ratio = norm(interp);
	vpair_ratios[vpair_i] = ratio;

	vpair_i ++;
  }

  matlab.setVar("nvpairs", nvpairs);
  matlab.setVar("vpair_lefts", vpair_lefts);
  matlab.setVar("vpair_rights", vpair_rights);
  matlab.setVar("vpair_vpids", vpair_vpids);
  matlab.setVar("vpair_ratios", vpair_ratios);

  matlab.setVar("vps", vps);
  // TODO:
  
  std::map<VertT, double> vert2depth;
  // TODO:
  return vert2depth;
}


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
  // for (auto &h : mesh.halfedges()) {
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
