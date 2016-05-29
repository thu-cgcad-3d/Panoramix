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
    const Mesh2 &mesh2d, misc::Matlab &matlab);

//// OptimizeWithoutOrientations
//void OptimizeWithoutOrientations(
//    Mesh3 &cur_reconstruction, const PerspectiveCamera &cur_cam,
//    const Mesh<VertHandle, HalfHandle, FaceHandle> &mesh_proxy,
//    const std::vector<SubMesh> &submeshes,
//    const std::unordered_set<VertHandle> &non_corner_vh_proxies);
//


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
