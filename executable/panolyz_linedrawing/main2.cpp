#include <QtWidgets>

#include "factor_graph.hpp"
#include "parallel.hpp"
#include "cache.hpp"
#include "clock.hpp"

#include "canvas.hpp"
#include "qttools.hpp"
#include "scene.hpp"
#include "singleton.hpp"
#include "gui_util.hpp"

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

#define DISABLED_main MACRO_CONCAT(main_, __COUNTER__)

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

// FaceVertDependency
struct FaceVertDependency {
  std::vector<FaceHandle> ordered_fhs;
  std::unordered_set<VertHandle> fundamental_vhs;
};

// PossibleKeyVanishingPoints
std::vector<Point2> PossibleKeyVanishingPoints(const Chain2 &chain) {
  assert(chain.size() > 2);
  if (chain.size() == 3) {
    return {};
  }
  if (chain.size() == 4) {
    return {Intersection(chain.edge(0).ray(), chain.edge(2).ray()),
            Intersection(chain.edge(1).ray(), chain.edge(3).ray())};
  }
  if (chain.size() % 2 == 0) {
    std::vector<Point2> vp_positions;
    vp_positions.reserve(chain.size() * chain.size() / 4);
	for (int i = 0; i < chain.size(); i++) {
      int j = (i + chain.size() / 2) % chain.size();
      vp_positions.push_back(
          Intersection(chain.edge(i).ray(), chain.edge(j).ray()));
    }
    return vp_positions;
  } else {
    std::vector<Point2> vp_positions;
    for (int i = 0; i < chain.size(); i++) {
      int j = (i + chain.size() / 2) % chain.size();
      vp_positions.push_back(
          Intersection(chain.edge(i).ray(), chain.edge(j).ray()));
      vp_positions.push_back(
          Intersection(chain.edge(i).ray(), chain.edge(j + 1).ray()));
    }
    return vp_positions;
  }
}

// PPFocalCandidate
struct PPFocalCandidate {
  Point2 pp;
  double focal;
};

// EstimateCameraParameters
std::vector<std::pair<std::set<int>, PPFocalCandidate>> EstimateCameraParameters(
    const Mesh2 &mesh2d,
    const Mesh<VertHandle, HalfHandle, FaceHandle> &mesh_proxy,
    const std::vector<SubMesh> &submeshes) {
  //// [Estimate PP & Focal Candidates from 2D Mesh]
  auto point2dAtProxy = [&mesh2d, &mesh_proxy](VertHandle vhInProxy) -> Point2 {
    return mesh2d.data(mesh_proxy.data(vhInProxy));
  };
  auto line2dAtProxy = [&mesh2d, &mesh_proxy,
                        point2dAtProxy](HalfHandle hhInProxy) -> Line2 {
    return Line2(point2dAtProxy(mesh_proxy.topo(hhInProxy).from()),
                 point2dAtProxy(mesh_proxy.topo(hhInProxy).to()));
  };
  Box2 box = BoundingBoxOfContainer(mesh2d.vertices());
  double scale = box.outerSphere().radius;

  // collect pp focal candidates
  std::vector<PPFocalCandidate> pp_focal_candidates;
  pp_focal_candidates.reserve(submeshes.size() * 3);

  for (int submesh_id = 0; submesh_id < submeshes.size(); submesh_id++) {
    // collect edge intersections in each face
    std::vector<Point2> interps;
    for (auto fh_proxy : submeshes[submesh_id].fhs) {
      auto &hh_proxies = mesh_proxy.topo(fh_proxy).halfedges;
      Chain2 corners;
      for (auto hh_proxy : hh_proxies) {
        corners.append(point2dAtProxy(mesh_proxy.topo(hh_proxy).to()));
      }
      auto key_vps = PossibleKeyVanishingPoints(corners);
      // remove all vps that lie on certain lines
      for (auto &key_vp : key_vps) {
        bool overlaped = false;
        for (int i = 0; i < corners.size(); i++) {
          Line2 edge = corners.edge(i);
          if (Distance(key_vp, edge) < 1e-3) {
            overlaped = true;
            break;
          }
        }
        if (!overlaped) {
          interps.push_back(key_vp);
        }
      }
    }

    for (int i = 0; i < interps.size(); i++) {
      const Point2 &p1 = interps[i];
      for (int j = i + 1; j < interps.size(); j++) {
        const Point2 &p2 = interps[j];
        for (int k = j + 1; k < interps.size(); k++) {
          const Point2 &p3 = interps[k];
          // compute pp and focal
          Point2 pp;
          double focal = 0.0;
          std::tie(pp, focal) = ComputePrinciplePointAndFocalLength(p1, p2, p3);
          if (HasValue(pp, IsInfOrNaN<double>) || IsInfOrNaN(focal)) {
            continue;
          }
          if (!IsBetween(focal, scale / 5.0, scale * 5.0) ||
              Distance(pp, box.center()) > scale) {
            continue;
          }
          pp_focal_candidates.push_back(PPFocalCandidate{pp, focal});
        }
      }
    }
  }

  std::sort(pp_focal_candidates.begin(), pp_focal_candidates.end(),
            [](auto &a, auto &b) { return a.focal < b.focal; });
  // naive clustering
  std::vector<std::pair<std::set<int>, PPFocalCandidate>> pp_focal_groups;
  {
    std::vector<int> ppFocalId2group(pp_focal_candidates.size(), -1);
    int ngroups = 0;
    RTreeMap<Vec3, int> ppFocalIdTree;
    for (int i = 0; i < pp_focal_candidates.size(); i++) {
      Vec3 coordinate =
          cat(pp_focal_candidates[i].pp, pp_focal_candidates[i].focal);
      const double thres = scale / 50.0;
      // find the nearest ppFocal sample point
      int nearestPPFocalCandId = -1;
      double min_dist = thres;
      ppFocalIdTree.search(BoundingBox(coordinate).expand(thres * 2),
                           [&nearestPPFocalCandId, &min_dist,
                            &coordinate](const std::pair<Vec3, int> &cand) {
                             double dist = Distance(cand.first, coordinate);
                             if (dist < min_dist) {
                               min_dist = dist;
                               nearestPPFocalCandId = cand.second;
                             }
                             return true;
                           });
      if (nearestPPFocalCandId != -1) { // if found, assign to the same group
        ppFocalId2group[i] = ppFocalId2group[nearestPPFocalCandId];
      } else { // otherwise, create a new group
        ppFocalId2group[i] = ngroups++;
      }
      ppFocalIdTree.emplace(coordinate, i);
    }

    pp_focal_groups.resize(ngroups);
    for (auto &g : pp_focal_groups) {
      g.second.focal = 0.0;
      g.second.pp = Point2();
    }
    for (int i = 0; i < ppFocalId2group.size(); i++) {
      auto &g = pp_focal_groups[ppFocalId2group[i]];
      g.first.insert(i);
      g.second.focal += pp_focal_candidates[i].focal;
      g.second.pp += pp_focal_candidates[i].pp;
    }
    for (auto &g : pp_focal_groups) {
      g.second.focal /= g.first.size();
      g.second.pp /= double(g.first.size());
    }

    std::sort(
        pp_focal_groups.begin(), pp_focal_groups.end(),
        [](auto &g1, auto &g2) { return g1.first.size() > g2.first.size(); });
  }

  return pp_focal_groups;
}

// EstimateVanishingPoints
std::vector<Point2> EstimateVanishingPoints(
    const Mesh2 &mesh2d, const Sizei &drawing_size,
    HandledTable<HalfHandle, int> *hh2d2edge_ptr = nullptr,
    std::vector<std::map<int, double>> *vp2edge_with_angles_ptr = nullptr,
    std::vector<std::vector<Scored<int>>> *edge2ordered_vp_and_angles_ptr =
        nullptr,
    std::vector<Line2> *edge2line_ptr = nullptr) {

  Box2 box = BoundingBoxOfContainer(mesh2d.vertices());
  double scale = box.outerSphere().radius;

  std::vector<std::pair<HalfHandle, HalfHandle>> edge2hh2ds;
  std::vector<Line2> edge2line;
  HandledTable<HalfHandle, int> hh2d2edge(mesh2d.internalHalfEdges().size(),
                                          -1);

  int nedges = 0;
  {
    for (auto &h : mesh2d.halfedges()) {
      auto hh = h.topo.hd;
      auto oppohh = h.topo.opposite;
      if (hh2d2edge[hh] == -1 && hh2d2edge[oppohh] == -1) {
        hh2d2edge[hh] = hh2d2edge[oppohh] = nedges;
        nedges++;
        edge2hh2ds.push_back(MakeOrderedPair(hh, oppohh));
        edge2line.push_back(Line2(mesh2d.data(mesh2d.topo(hh).from()),
                                  mesh2d.data(mesh2d.topo(hh).to())));
      }
    }
    assert(edge2hh2ds.size() == nedges && edge2line.size() == nedges);
  }

  std::vector<Point2> vp_positions;
  int nvps = 0;
  std::vector<std::vector<Scored<int>>> edge2ordered_vp_and_angles(nedges);

  std::vector<Point2> intersections;
  std::vector<std::pair<int, int>> intersection2edges;
  intersections.reserve(nedges * (nedges - 1) / 2);
  intersection2edges.reserve(nedges * (nedges - 1) / 2);
  for (int i = 0; i < nedges; i++) {
    const Line2 &linei = edge2line[i];
    for (int j = i + 1; j < nedges; j++) {
      const Line2 &linej = edge2line[j];
      Point2 interp = Intersection(linei.ray(), linej.ray());
      if (IsFuzzyParallel(linei.direction(), linej.direction(), 1e-8)) {
        interp = normalize(linei.direction()) * 1e8;
      }
      if (std::min(Distance(interp, linei), Distance(interp, linej)) <=
          scale / 10.0) {
        continue;
      }
      intersections.push_back(interp);
      intersection2edges.emplace_back(i, j);
    }
  }
  assert(intersections.size() == intersection2edges.size());

  std::vector<int> intersection2raw_vp(intersections.size(), -1);
  RTreeMap<Point2, int> intersection_tree;
  for (int i = 0; i < intersections.size(); i++) {
    const double thres = scale / 30.0;
    auto &p = intersections[i];
    int nearest_intersection_id = -1;
    double min_dist = thres;
    intersection_tree.search(
        BoundingBox(p).expand(thres * 2),
        [&nearest_intersection_id, &min_dist,
         &p](const std::pair<Point2, int> &location_and_intersection_id) {
          double dist = Distance(location_and_intersection_id.first, p);
          if (dist < min_dist) {
            min_dist = dist;
            nearest_intersection_id = location_and_intersection_id.second;
          }
          return true;
        });
    if (nearest_intersection_id != -1) {
      intersection2raw_vp[i] = intersection2raw_vp[nearest_intersection_id];
    } else {
      intersection2raw_vp[i] = nvps++;
    }
    intersection_tree.emplace(p, i);
  }

  // further merge vps
  // if any two vps share two or more edges, merge them
  std::vector<std::set<int>> raw_vp2edges(nvps);
  for (int i = 0; i < intersection2raw_vp.size(); i++) {
    int vpid = intersection2raw_vp[i];
    raw_vp2edges[vpid].insert(i);
  }
  std::vector<std::set<int>> raw_vp2raw_vp_should_merge(nvps);
  for (int vp1 = 0; vp1 < nvps; vp1++) {
    auto &edges1 = raw_vp2edges[vp1];
    for (int vp2 = vp1 + 1; vp2 < nvps; vp2++) {
      auto &edges2 = raw_vp2edges[vp2];
      std::set<int> common_edges;
      std::set_intersection(edges1.begin(), edges1.end(), edges2.begin(),
                            edges2.end(),
                            std::inserter(common_edges, common_edges.begin()));
      if (common_edges.size() >= 2) {
        raw_vp2raw_vp_should_merge[vp1].insert(vp2);
        raw_vp2raw_vp_should_merge[vp2].insert(vp1);
      }
    }
  }
  std::map<int, std::set<int>> new_vp2raw_vps;
  std::vector<int> raw_vp2new_vp(nvps, -1);
  std::vector<int> raw_vp_ids(nvps);
  std::iota(raw_vp_ids.begin(), raw_vp_ids.end(), 0);
  nvps = ConnectedComponents(
      raw_vp_ids.begin(), raw_vp_ids.end(),
      [&raw_vp2raw_vp_should_merge](int vp) -> const std::set<int> & {
        return raw_vp2raw_vp_should_merge.at(vp);
      },
      [&new_vp2raw_vps, &raw_vp2new_vp](int raw_vp, int newVP) {
        new_vp2raw_vps[newVP].insert(raw_vp);
        raw_vp2new_vp[raw_vp] = newVP;
      });

  vp_positions.resize(nvps, Origin<2>());
  std::vector<std::set<int>> vp2intersections(nvps);
  for (int i = 0; i < intersection2raw_vp.size(); i++) {
    int raw_vp = intersection2raw_vp[i];
    int newVP = raw_vp2new_vp[raw_vp];
    vp_positions[newVP] += intersections[i]; // TODO: what if some
                                             // intersections are oppsite far
                                             // points?
    vp2intersections[newVP].insert(i);
  }
  for (int i = 0; i < nvps; i++) {
    vp_positions[i] /= double(vp2intersections[i].size());
    assert(!HasValue(vp_positions[i], [](auto &&e) { return IsInfOrNaN(e); }));
  }
  // erase invalid vp_positions
  vp_positions.erase(std::remove_if(vp_positions.begin(), vp_positions.end(),
                                    [](const Point2 &pos) {
                                      return IsInfOrNaN(pos[0]) ||
                                             IsInfOrNaN(pos[1]);
                                    }),
                     vp_positions.end());
  nvps = vp_positions.size();

  // initial edge vp bindings
  std::vector<std::map<int, double>> vp2edge_with_angles(nvps);
  std::vector<bool> vp_is_good(nvps, true);
  for (int vp = 0; vp < nvps; vp++) {
    auto &vp_pos = vp_positions[vp];
    for (int edge = 0; edge < nedges; edge++) {
      auto &line = edge2line[edge];
      double lambda = ProjectionOfPointOnLine(vp_pos, line).ratio;
      static const double thres = 0.1;
      if (lambda >= -thres && lambda <= 1.0 + thres) {
        continue;
      }
      double angle =
          AngleBetweenUndirected(line.direction(), vp_pos - line.center());
      assert(!IsInfOrNaN(angle));
      static const double theta = DegreesToRadians(5); ////// TODO
      if (angle >= theta) {
        continue;
      }
      vp2edge_with_angles[vp][edge] = angle;
    }
    vp_is_good[vp] = vp2edge_with_angles[vp].size() >= 3;
  }

  if (false) {
    for (int i = 0; i < std::min(10, nvps); i++) {
      Image3ub im(drawing_size, Vec3ub(255, 255, 255));
      auto canvas = gui::MakeCanvas(im);
      canvas.color(gui::LightGray);
      canvas.thickness(2);
      for (auto &line : edge2line) {
        canvas.add(line);
      }
      canvas.color(gui::Gray);
      canvas.thickness(1);
      for (auto &edge_with_angle : vp2edge_with_angles[i]) {
        canvas.add(edge2line[edge_with_angle.first].ray());
      }
      canvas.color(gui::Black);
      for (auto &edge_with_angle : vp2edge_with_angles[i]) {
        canvas.add(edge2line[edge_with_angle.first]);
      }
      canvas.show(0, "before removing bad vps: raw vp_" + std::to_string(i));
    }
  }

  // remove bad vps
  int newvp = 0;
  for (int oldvp = 0; oldvp < nvps; oldvp++) {
    if (vp_is_good[oldvp]) {
      // update vp_positions
      // update vp2edge_with_angles
      vp_positions[newvp] = vp_positions[oldvp];
      vp2edge_with_angles[newvp] = std::move(vp2edge_with_angles[oldvp]);
      newvp++;
    }
  }
  nvps = newvp;
  vp_positions.resize(nvps);
  vp2edge_with_angles.resize(nvps);

  // construct edge2OrderedVPs
  for (int vp = 0; vp < nvps; vp++) {
    for (auto &edgeAndAngle : vp2edge_with_angles[vp]) {
      int edge = edgeAndAngle.first;
      double angle = edgeAndAngle.second;
      edge2ordered_vp_and_angles[edge].push_back(ScoreAs(vp, angle));
    }
  }
  for (auto &vp_and_angles : edge2ordered_vp_and_angles) {
    std::sort(vp_and_angles.begin(), vp_and_angles.end());
  }

  if (false) {
    for (int i = 0; i < std::min(10, nvps); i++) {
      Image3ub im(drawing_size, Vec3ub(255, 255, 255));
      auto canvas = gui::MakeCanvas(im);
      canvas.color(gui::LightGray);
      canvas.thickness(2);
      for (auto &line : edge2line) {
        canvas.add(line);
      }
      canvas.color(gui::Gray);
      canvas.thickness(2);
      for (auto &edge_with_angle : vp2edge_with_angles[i]) {
        canvas.add(edge2line[edge_with_angle.first].ray());
      }
      canvas.color(gui::Black);
      for (auto &edge_with_angle : vp2edge_with_angles[i]) {
        canvas.add(edge2line[edge_with_angle.first]);
      }
      canvas.show(0, "raw vp_" + std::to_string(i));
    }
  }

  if (hh2d2edge_ptr) {
    *hh2d2edge_ptr = std::move(hh2d2edge);
  }
  if (vp2edge_with_angles_ptr) {
    *vp2edge_with_angles_ptr = std::move(vp2edge_with_angles);
  }
  if (edge2ordered_vp_and_angles_ptr) {
    *edge2ordered_vp_and_angles_ptr = std::move(edge2ordered_vp_and_angles);
  }
  if (edge2line_ptr) {
    *edge2line_ptr = std::move(edge2line);
  }
  return vp_positions;
}

// EstimateEdgeOrientations
void EstimateEdgeOrientations(
    const Mesh2 &mesh2d, const std::vector<Point2> &vp_positions,
    const HandledTable<HalfHandle, int> &hh2d2edge,
    const std::vector<std::map<int, double>> &vp2edge_with_angles,
    const std::vector<std::vector<Scored<int>>> &edge2ordered_vp_and_angles,
    std::vector<int> &edge2vp, std::vector<std::vector<int>> &vp2edges) {

  Box2 box = BoundingBoxOfContainer(mesh2d.vertices());
  double scale = box.outerSphere().radius;

  int nedges = edge2ordered_vp_and_angles.size();
  int nvps = vp2edge_with_angles.size();
  edge2vp.resize(nedges, -1);
  vp2edges.resize(nvps);

  // construct a factor graph to optimize edge-vp bindings
  FactorGraph fg;
  std::vector<int> edge2vh(nedges, -1);
  {
    for (int edge = 0; edge < nedges; edge++) {
      auto vc =
          fg.addVarCategory(edge2ordered_vp_and_angles[edge].size() + 1, 1.0);
      edge2vh[edge] = fg.addVar(vc);
    }

    // potential 1: the edge should bind to some vp, should prefer better
    // scored
    for (int edge = 0; edge < nedges; edge++) {
      auto vh = edge2vh[edge];
      auto &relatedVPAndAngles = edge2ordered_vp_and_angles[edge];
      auto fc = fg.addFactorCategory(
          [&relatedVPAndAngles, nedges](const std::vector<int> & varlabels,
                                        void *givenData) -> double {
            assert(varlabels.size() == 1);
            int label = varlabels[0];
            assert(label <= relatedVPAndAngles.size());
            const double K = 50.0 / nedges;
            if (label == relatedVPAndAngles.size()) { // not bind to any vp
              return K;
            }
            double angle = relatedVPAndAngles[label].score;
            assert(!IsInfOrNaN(angle));
            return (1.0 - Gaussian(angle, DegreesToRadians(3))) * K;
          },
          1.0);
	  fg.addFactor(fc, {vh});
    }

    // potential 2: two adjacent edges should not bind to a near vp
    int ncorners = 0;
    for (auto &f : mesh2d.faces()) {
      ncorners += f.topo.halfedges.size();
    }
    for (auto &f : mesh2d.faces()) {
      auto &hhs = f.topo.halfedges;
      for (int i = 0; i < hhs.size(); i++) {
        auto hh1 = hhs[i];
        auto hh2 = hhs[(i + 1) % hhs.size()];
        int edge1 = hh2d2edge[hh1];
        int edge2 = hh2d2edge[hh2];
        auto &relatedVPAndAngles1 = edge2ordered_vp_and_angles[edge1];
        auto &relatedVPAndAngles2 = edge2ordered_vp_and_angles[edge2];
        auto vh1 = edge2vh[edge1];
        auto vh2 = edge2vh[edge2];
        auto fc = fg.addFactorCategory(
            [edge1, edge2, &relatedVPAndAngles1, &relatedVPAndAngles2,
             &vp_positions, ncorners, scale](const std::vector<int> &varlabels,
                                             void *givenData) -> double {
              assert(varlabels.size() == 2);
              int bindedVP1 =
                  varlabels[0] == relatedVPAndAngles1.size()
                      ? -1
                      : (relatedVPAndAngles1[varlabels[0]].component);
              int bindedVP2 =
                  varlabels[1] == relatedVPAndAngles2.size()
                      ? -1
                      : (relatedVPAndAngles2[varlabels[1]].component);
              if (bindedVP1 == -1 || bindedVP2 == -1) {
                return 0;
              }
              auto &vpPos1 = vp_positions[bindedVP1];
              auto &vpPos2 = vp_positions[bindedVP2];
              const double thres = scale / 10.0;
              const double K = 10.0 / ncorners;
              if (Distance(vpPos1, vpPos2) < thres) { // todo
                return K;
              }
              return 0.0;
            },
            1.0);
        fg.addFactor(fc, {vh1, vh2});
      }
    }

    // potential 3: the vp_positions of edges sharing a same face should lie
    // on
    // the same
    // line (the vanishing line of the face)
    int ntris = 0;
    for (auto &f : mesh2d.faces()) {
      auto &hhs = f.topo.halfedges;
      if (hhs.size() <= 3) {
        continue;
      }
      ntris = hhs.size() > 4 ? (2 * hhs.size()) : hhs.size();
    }
    for (auto &f : mesh2d.faces()) {
      auto &hhs = f.topo.halfedges;
      if (hhs.size() <= 3) {
        continue;
      }
      for (int i = 0; i < hhs.size(); i++) {
        int maxGap = hhs.size() > 4 ? 2 : 1;
        for (int gap = 1; gap <= maxGap; gap++) {
          int prevEdge = hh2d2edge[hhs[(i + hhs.size() - gap) % hhs.size()]];
          int edge = hh2d2edge[hhs[i]];
          int nextEdge = hh2d2edge[hhs[(i + gap) % hhs.size()]];
          auto fc = fg.addFactorCategory(
              [&edge2ordered_vp_and_angles, &vp_positions, prevEdge, edge,
               nextEdge, ntris](const std::vector<int> &varlabels,
                                void *givenData) -> double {
                assert(varlabels.size() == 3);
                int vp1 =
                    varlabels[0] == edge2ordered_vp_and_angles[prevEdge].size()
                        ? -1
                        : (edge2ordered_vp_and_angles[prevEdge][varlabels[0]]
                               .component);
                if (vp1 == -1) {
                  return 0.0;
                }
                int vp2 =
                    varlabels[1] == edge2ordered_vp_and_angles[edge].size()
                        ? -1
                        : (edge2ordered_vp_and_angles[edge][varlabels[1]]
                               .component);
                if (vp2 == -1) {
                  return 0.0;
                }
                int vp3 =
                    varlabels[2] == edge2ordered_vp_and_angles[nextEdge].size()
                        ? -1
                        : (edge2ordered_vp_and_angles[nextEdge][varlabels[2]]
                               .component);
                if (vp3 == -1) {
                  return 0.0;
                }
                if (vp1 == vp2 || vp2 == vp3 || vp1 == vp3) {
                  return 0.0;
                }
                double angle = AngleBetweenUndirected(
                    vp_positions[vp1] - vp_positions[vp2],
                    vp_positions[vp3] - vp_positions[vp2]);
                assert(!IsInfOrNaN(angle));
                const double K = 30.0 / ntris;
                return (1.0 - Gaussian(angle, DegreesToRadians(10))) * K;
              },
              1.0);
          fg.addFactor(fc,
                       {edge2vh[prevEdge], edge2vh[edge], edge2vh[nextEdge]});
        }
      }
    }

    //// potential 4: a vp is good only when edges are bound to it!
    //// a vp should be with either no edges or >= 3 edges
    // for (int vpid = 0; vpid < nvps; vpid++) {
    //  std::vector<FactorGraph::VarHandle> relatedVhs;
    //  relatedVhs.reserve(vp2edges[vpid].size());
    //  for (int edge : vp2edges[vpid]) {
    //    relatedVhs.push_back(edge2vh[edge]);
    //  }
    //  auto fc = fg.addFactorCategory(
    //      [&edge2ordered_vp_and_angles, vpid, &vp2edges](
    //          const int *varlabels, size_t nvar,
    //          FactorGraph::FactorCategoryId fcid, void *givenData) ->
    //          double
    //          {
    //        int bindedEdges = 0;
    //        auto &relatedEdges = vp2edges[vpid];
    //        assert(nvar == relatedEdges.size());
    //        for (int i = 0; i < relatedEdges.size(); i++) {
    //          int edge = relatedEdges[i];
    //          int edgeLabel = varlabels[i];
    //          if (edgeLabel == edge2ordered_vp_and_angles[edge].size()) {
    //            continue; // not bound to any vp
    //          }
    //          int edgeBindedVPId =
    //          edge2ordered_vp_and_angles[edge][edgeLabel].component;
    //          if (edgeBindedVPId == vpid) {
    //            bindedEdges++;
    //          }
    //        }
    //        // todo
    //        if (bindedEdges == 1 || bindedEdges == 2) {
    //          return 10.0;
    //        }
    //        return 0.0;
    //      },
    //      1.0);
    //  fg.addFactor(relatedVhs.begin(), relatedVhs.end(), fc);
    //}
  }

  // solve the factor graph
  auto result = fg.solve(5, 1, [](int epoch, double energy, double denergy,
                                  const std::vector<int> &results) -> bool {
    Println("epoch: ", epoch, "  energy: ", energy);
    return true;
  });

  for (int edge = 0; edge < nedges; edge++) {
    int id = result[edge2vh[edge]];
    if (id == edge2ordered_vp_and_angles[edge].size()) {
      continue;
    }
    edge2vp[edge] = edge2ordered_vp_and_angles[edge][id].component;
    vp2edges[edge2vp[edge]].push_back(edge);
  }

  // invalidate the edge bindings for vps who have <= 3 edges
  for (int vp = 0; vp < nvps; vp++) {
    if (vp2edges[vp].size() <= 3) { //
      for (int edge : vp2edges[vp]) {
        edge2vp[edge] = -1;
      }
      vp2edges[vp].clear();
    }
  }
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

// ReconstructWithOrientations
Mesh3 ReconstructWithOrientations(
    const std::vector<Line2> &edge2line, const std::vector<int> &edge2vp,
    const std::vector<Vec3> &vp2dir,
    const HandledTable<HalfHandle, int> &hh2d2edge,
    const PerspectiveCamera &cur_cam,
    const Mesh<VertHandle, HalfHandle, FaceHandle> &mesh_proxy,
    const Mesh2 &mesh2d, misc::Matlab &matlab) {

  int nedges = edge2line.size();
  assert(edge2vp.size() == nedges);

  // get the vh2d2reconstructedPosition data
  // entity
  using SupportingPlane = PIConstraintGraph::Entity::SupportingPlane;
  struct EntityBase {
    SupportingPlane supportingPlane;
    EntityBase(const SupportingPlane &sp) : supportingPlane(sp) {}
    virtual ~EntityBase() {}
    virtual FaceHandle fh_proxy() const { return FaceHandle(); }
    virtual int edge() const { return -1; }
    bool isFace() const { return fh_proxy().valid(); }
    bool isEdge() const { return edge() != -1; }
  };
  struct FaceEntity : EntityBase {
    FaceHandle fh;
    FaceEntity(FaceHandle fh, const Vec3 &center)
        : EntityBase(SupportingPlane(SegControl{-1, -1}, center, {})), fh(fh) {}
    virtual ~FaceEntity() {}
    virtual FaceHandle fh_proxy() const override { return fh; }
  };
  struct EdgeEntity : EntityBase {
    int e;
    EdgeEntity(int e, const Classified<Line3> &line,
               const std::vector<Vec3> &vp2dir)
        : EntityBase(SupportingPlane(line, vp2dir)), e(e) {}
    virtual ~EdgeEntity() {}
    virtual int edge() const override { return e; }
  };

  std::vector<std::unique_ptr<EntityBase>> entities;

  std::vector<int> edge2ent(nedges, -1);
  HandledTable<FaceHandle, int> fhProxy2ent(mesh_proxy.internalFaces().size(),
                                            -1);

  // install supporting planes from edges and faces
  // from edges
  for (int edge = 0; edge < nedges; edge++) {
    auto &line2 = edge2line[edge];
    Line3 line3(normalize(cur_cam.direction(line2.first)),
                normalize(cur_cam.direction(line2.second)));
    int vp = edge2vp[edge];
    entities.push_back(
        std::make_unique<EdgeEntity>(edge, ClassifyAs(line3, vp), vp2dir));
    int ent = entities.size() - 1;
    edge2ent[edge] = ent;
  }
  // from faces
  HandledTable<FaceHandle, Point2> fhProxy2center2d(
      mesh_proxy.internalFaces().size());
  for (auto &f : mesh_proxy.faces()) {
    Point2 center2 = Origin<2>();
    for (auto hh : f.topo.halfedges) {
      center2 += mesh2d.data(mesh_proxy.data(mesh_proxy.topo(hh).to()));
    }
    center2 /= double(f.topo.halfedges.size());
    fhProxy2center2d[f.topo.hd] = center2;
    entities.push_back(std::make_unique<FaceEntity>(
        f.topo.hd, normalize(cur_cam.direction(center2))));
    int ent = entities.size() - 1;
    fhProxy2ent[f.topo.hd] = ent;
  }

  //// install connections
  // std::map<std::pair<int, int>, std::vector<Vec3>> ents2anchors;
  // for (auto &f : mesh_proxy.faces()) {
  //  for (auto hh_proxy : f.topo.halfedges) {
  //    HalfHandle hh2d = mesh_proxy.data(hh_proxy);
  //    if (hh2d.invalid()) {
  //      hh2d = mesh_proxy.data(mesh_proxy.topo(hh_proxy).opposite);
  //    }
  //    assert(hh2d.valid());
  //    int edge = hh2d2edge[hh2d];

  //    auto &line2 = edge2line[edge];
  //    ents2anchors[std::make_pair(edge2ent[edge], fhProxy2ent[f.topo.hd])]
  //    =
  //        {normalize(cur_cam.direction(line2.first)),
  //         normalize(cur_cam.direction(line2.second))};
  //  }
  //}

  //// install face angles constraints
  // std::vector<std::pair<int, int>> adjFaceEnts;
  // std::vector<int> adjFace2SubMeshId;
  // std::vector<bool> adjFaceOverlapInView;
  // for (auto &half : mesh_proxy.halfedges()) {
  //  FaceHandle fh1 = half.topo.face;
  //  FaceHandle fh2 = mesh_proxy.topo(half.topo.opposite).face;

  //  bool onSameSide = IsOnLeftSide(faceCorners[0], line.first,
  //  line.second) ==
  //                    IsOnLeftSide(faceCorners[1], line.first,
  //                    line.second);

  //  adjFaceOverlapInView.push_back(onSameSide);
  //}
  // int nadjFaces = adjFaceEnts.size();

  // Start Building Matrices
  // vert start position in variable vector
  std::vector<int> ent2varPosition(entities.size(), -1);
  std::vector<int> ent2nvar(entities.size(), -1);
  std::vector<DenseMatd> ent2matFromVarToPlaneCoeffs(entities.size());
  std::vector<int> var2ent;
  int nvars = 0;
  for (int ent = 0; ent < entities.size(); ent++) {
    auto &e = *entities[ent];
    int nvar = e.supportingPlane.dof;
    ent2nvar[ent] = nvar;
    ent2varPosition[ent] = nvars;
    ent2matFromVarToPlaneCoeffs[ent] =
        e.supportingPlane.matFromVarsToPlaneCoeffs();
    var2ent.insert(var2ent.end(), (size_t)nvar, ent);
    nvars += nvar;
  }

  // P : [(3 * nents) x nvars]
  // P * X -> plane coefficients
  int nents = entities.size();
  std::vector<SparseMatElementd> Ptriplets;
  for (int ent = 0; ent < nents; ent++) {
    int varpos = ent2varPosition[ent];
    // [3 x k]
    auto &matFromVarsToPlaneCoeffs = ent2matFromVarToPlaneCoeffs[ent];
    assert(matFromVarsToPlaneCoeffs.rows == 3 &&
           matFromVarsToPlaneCoeffs.cols == ent2nvar[ent]);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < ent2nvar[ent]; j++) {
        assert(varpos + j < nvars);
        Ptriplets.emplace_back(i + ent * 3, varpos + j,
                               matFromVarsToPlaneCoeffs(i, j));
      }
    }
  }

  // anchorPositions: [nanchors x 3]
  int nanchors = 0;
  std::vector<SparseMatElementd> anchorPositions;
  HandledTable<VertHandle, int> vhProxy2anchor(
      mesh_proxy.internalVertices().size(), -1);
  for (auto &vert : mesh_proxy.vertices()) {
    auto vh2d = vert.data;
    Vec3 anchor = normalize(cur_cam.direction(mesh2d.data(vh2d)));
    for (int i = 0; i < 3; i++) {
      anchorPositions.emplace_back(nanchors, i, anchor[i]);
    }
    vhProxy2anchor[vert.topo.hd] = nanchors;
    nanchors++;
  }

  // connection2anchor: [nconnection -> nanchors],
  // connection2leftEnt, connection2rightEnt:
  // [nconnection -> nents]
  std::vector<double> connection2anchor, connection2leftEnt,
      connection2rightEnt;
  // angleHead2anchor, angleTailA2anchor, angleTailB2anchor
  std::vector<double> angleHead2anchor, angleTailA2anchor, angleTailB2anchor;
  // angle2faceEnt
  std::vector<double> angle2faceEnt;

  int nconnections = 0;
  int nangles = 0;
  for (auto &f : mesh_proxy.faces()) {
    auto &hhProxies = f.topo.halfedges;
    for (int i = 0; i < hhProxies.size(); i++) {
      auto hh_proxy = hhProxies[i];
      HalfHandle hh2d = mesh_proxy.data(hh_proxy);
      if (hh2d.invalid()) {
        hh2d = mesh_proxy.data(mesh_proxy.topo(hh_proxy).opposite);
      }
      assert(hh2d.valid());
      int edge = hh2d2edge[hh2d];
      int edgeEnt = edge2ent[edge];
      int faceEnt = fhProxy2ent[f.topo.hd];

      VertHandle curVhProxy = mesh_proxy.topo(hh_proxy).to();
      int curAnchorId = vhProxy2anchor[curVhProxy];
      int prevAnchorId = vhProxy2anchor[mesh_proxy.topo(hh_proxy).from()];
      assert(mesh_proxy
                 .topo(hhProxies[(i + hhProxies.size() - 1) % hhProxies.size()])
                 .to() == mesh_proxy.topo(hh_proxy).from());
      int nextAnchorId =
          vhProxy2anchor[mesh_proxy.topo(hhProxies[(i + 1) % hhProxies.size()])
                             .to()];

      connection2leftEnt.push_back(edgeEnt);
      connection2rightEnt.push_back(faceEnt);
      connection2anchor.push_back(curAnchorId);
      connection2leftEnt.push_back(edgeEnt);
      connection2rightEnt.push_back(faceEnt);
      connection2anchor.push_back(prevAnchorId);

      angleHead2anchor.push_back(curAnchorId);
      angleTailA2anchor.push_back(prevAnchorId);
      angleTailB2anchor.push_back(nextAnchorId);
      angle2faceEnt.push_back(faceEnt);

      nconnections += 2;
      nangles++;
    }
  }
  assert(angle2faceEnt.size() == nangles);

  // solve!
  matlab << "clear()";

  // matrices
  matlab.setVar("P", MakeSparseMatFromElements(
                         3 * nents, nvars, Ptriplets.begin(), Ptriplets.end()));
  matlab.setVar("AnchorPositions",
                MakeSparseMatFromElements(nanchors, 3, anchorPositions.begin(),
                                          anchorPositions.end()));
  // index mappings
  matlab.setVar("connection2anchor", DenseMatd(connection2anchor));
  matlab.setVar("connection2leftEnt", DenseMatd(connection2leftEnt));
  matlab.setVar("connection2rightEnt", DenseMatd(connection2rightEnt));

  matlab.setVar("angleHead2anchor", DenseMatd(angleHead2anchor));
  matlab.setVar("angleTailA2anchor", DenseMatd(angleTailA2anchor));
  matlab.setVar("angleTailB2anchor", DenseMatd(angleTailB2anchor));
  matlab.setVar("angle2faceEnt", DenseMatd(angle2faceEnt));

  // sizes
  matlab.setVar("nvars", nvars);
  matlab.setVar("nents", nents);
  matlab.setVar("nanchors", nanchors);
  matlab.setVar("nconnections", nconnections);
  matlab.setVar("nangles", nangles);

  double minE = std::numeric_limits<double>::infinity();
  matlab << "FinalX = zeros(nvars, 1);";
  static const int maxIter = 10;

  matlab << "Prev_InverseDepthsOnLeftOfConnections = ones(nconnections, 1);";
  matlab << "Prev_InverseDepthsOnRightOfConnections = ones(nconnections, 1);";

  matlab << "Prev_Scalar1 = ones(nangles, 1);";
  matlab << "Prev_Scalar2 = ones(nangles, 1);";
  matlab << "Prev_Scalar3 = ones(nangles, 1);";
  // matlab << "Prev_AvgAngleCos = zeros(nangles, 1);";

  matlab.setPrintMessage(false);
  for (int t = 0; t < maxIter; t++) {
    matlab << "cvx_begin quiet";
    matlab << "variable X(nvars);";

    auto constructProblem = [&matlab](const std::string &objectiveVarName) {
      matlab << "PlaneEqs = reshape(P * X, [3 nents])';"; // [nents x 3]

      matlab << "PlaneEqsOnLeftOfConnections = "
                "PlaneEqs(connection2leftEnt + 1, :);";
      matlab << "PlaneEqsOnRightOfConnections = "
                "PlaneEqs(connection2rightEnt + 1, :);";

      matlab << "ConnectionPositions = "
                "AnchorPositions(connection2anchor + 1, :);";
      matlab << "InverseDepthsOnLeftOfConnections = "
                "dot(PlaneEqsOnLeftOfConnections, ConnectionPositions, 2);";
      matlab << "InverseDepthsOnRightOfConnections = "
                "dot(PlaneEqsOnRightOfConnections, ConnectionPositions, 2);";

      // energy of connections
      matlab << "EnergyOfConnections = "
                "sum_square((InverseDepthsOnLeftOfConnections - "
                "InverseDepthsOnRightOfConnections) ./ "
                "Prev_InverseDepthsOnLeftOfConnections ./ "
                "Prev_InverseDepthsOnRightOfConnections);";
      const std::string trueEnergyOfConnectionsExpr =
          "sum_square((InverseDepthsOnLeftOfConnections - "
          "InverseDepthsOnRightOfConnections) ./ "
          "InverseDepthsOnLeftOfConnections ./ "
          "InverseDepthsOnRightOfConnections)";

      matlab << "PlaneEqsOfAngles = PlaneEqs(angle2faceEnt + 1, :);";

      matlab << "AngleHeadPositions = "
                "AnchorPositions(angleHead2anchor + 1, :);";
      matlab << "AngleTailAPositions = "
                "AnchorPositions(angleTailA2anchor + 1, :);";
      matlab << "AngleTailBPositions = "
                "AnchorPositions(angleTailB2anchor + 1, :);";

      // m : AngleTailAPositions
      // n : AngleTailBPositions
      matlab << "Scalar1 = dot(AngleHeadPositions, PlaneEqsOfAngles, 2);";
      matlab << "Scalar2 = dot(AngleTailBPositions, PlaneEqsOfAngles, 2);";
      matlab << "Scalar3 = dot(AngleTailAPositions, PlaneEqsOfAngles, 2);";

      // px : AngleHeadPositions(:, 1)
      // py : AngleHeadPositions(:, 2)
      // pz : AngleHeadPositions(:, 3)

      // mx : AngleTailAPositions(:, 1)
      // my : AngleTailAPositions(:, 2)
      // mz : AngleTailAPositions(:, 3)

      // nx : AngleTailBPositions(:, 1)
      // ny : AngleTailBPositions(:, 2)
      // nz : AngleTailBPositions(:, 3)

      // / mx   px \ / nx   px \   / my   py \ / ny   py \   / mz   pz \ / nz   pz \
        // | -- - -- | | -- - -- | + | -- - -- | | -- - -- | + | -- - -- | | -- - -- |
      // \ #3   #1 / \ #2   #1 /   \ #3   #1 / \ #2   #1 /   \ #3   #1 / \
          // #2
      // #1 /
      matlab << "AngleCosineNumerator1 = "
                "(Scalar1 .* AngleTailAPositions(:, 1) - "
                " Scalar3 .* AngleHeadPositions(:, 1)) .*"
                "(Prev_Scalar1 .* AngleTailBPositions(:, 1) - "
                "Prev_Scalar2 .* AngleHeadPositions(:, 1));";
      matlab << "AngleCosineNumerator2 = "
                "(Scalar1 .* AngleTailAPositions(:, 2) - "
                " Scalar3 .* AngleHeadPositions(:, 2)) .*"
                "(Prev_Scalar1 .* AngleTailBPositions(:, 2) - "
                "Prev_Scalar2 .* AngleHeadPositions(:, 2));";
      matlab << "AngleCosineNumerator3 = "
                "(Scalar1 .* AngleTailAPositions(:, 3) - "
                " Scalar3 .* AngleHeadPositions(:, 3)) .*"
                "(Prev_Scalar1 .* AngleTailBPositions(:, 3) - "
                "Prev_Scalar2 .* AngleHeadPositions(:, 3));";

      matlab << "AngleCosines = "
                "(AngleCosineNumerator1 + AngleCosineNumerator2 + "
                "AngleCosineNumerator3) "
                "./ Prev_Scalar1 ./ Prev_Scalar1 "
                "./ Prev_Scalar2 ./ Prev_Scalar3;";

      // restrict face angles to be orhtogonal, aka, cosine values to be
      // zero
      matlab << "EnergyOfOrthogonalFaceAngles = sum_square(AngleCosines);";

      // energy
      matlab << (objectiveVarName + " = 1e3 * EnergyOfConnections + "
                                    "EnergyOfOrthogonalFaceAngles;");
    };

    constructProblem("ObjectiveEnergy");
    matlab << "minimize 1e3 * ObjectiveEnergy";
    matlab << "subject to";
    matlab << "   ones(nconnections, 1) <= "
              "InverseDepthsOnLeftOfConnections;";
    matlab << "   ones(nconnections, 1) <= "
              "InverseDepthsOnRightOfConnections;";
    matlab << "cvx_end";

    matlab << "Prev_InverseDepthsOnLeftOfConnections = "
              "InverseDepthsOnLeftOfConnections ./ "
              "norm(InverseDepthsOnLeftOfConnections);";
    matlab << "Prev_InverseDepthsOnRightOfConnections = "
              "InverseDepthsOnRightOfConnections ./ "
              "norm(InverseDepthsOnRightOfConnections);";
    matlab << "Prev_Scalar1 = Scalar1;";
    matlab << "Prev_Scalar2 = Scalar2;";
    matlab << "Prev_Scalar3 = Scalar3;";

    constructProblem("TrueObjectiveEnergy");

    double curEnergy = matlab.var("TrueObjectiveEnergy");
    Println("cur energy - ", curEnergy);
    if (IsInfOrNaN(curEnergy)) {
      break;
    }
    if (curEnergy < minE) {
      minE = curEnergy;
      matlab << "FinalX = 2 * X ./ median(InverseDepthsOnLeftOfConnections + "
                "InverseDepthsOnRightOfConnections);";
      matlab << "norm(FinalX)";
    }
  }
  matlab.setPrintMessage(true);

  // use the solved X to recover supporting planes in each ent
  matlab << "PlaneEqsVector = P * FinalX;"; // [nents*3 x 1]
  DenseMatd planeEqsVector = matlab.var("PlaneEqsVector");
  for (int ent = 0; ent < nents; ent++) {
    auto &entity = *entities[ent];
    entity.supportingPlane.reconstructed = Plane3FromEquation(
        planeEqsVector(ent * 3 + 0), planeEqsVector(ent * 3 + 1),
        planeEqsVector(ent * 3 + 2));
  }

  if (true) {
    // show reconstructed edge lines directly
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XLines;
    sb.installingOptions().discretizeOptions.color(gui::Black);
    sb.installingOptions().lineWidth = 10;
    for (int edge = 0; edge < nedges; edge++) {
      auto &line2 = edge2line[edge];
      auto &plane = entities[edge2ent[edge]]->supportingPlane.reconstructed;
      Line3 line3(
          Intersection(Ray3(cur_cam.eye(), cur_cam.direction(line2.first)),
                       plane),
          Intersection(Ray3(cur_cam.eye(), cur_cam.direction(line2.second)),
                       plane));
      sb.add(line3);
    }
    sb.show(true, true, gui::RenderOptions()
                            .winName("After CLLS")
                            .fixUpDirectionInCameraMove(false));
  }

  // reconstruct vertex positions
  Mesh3 curReconstruction =
      Transform(mesh2d, [](const Point2 &) { return Origin(); });
  HandledTable<VertHandle, int> vh2d2faceProxyCount(
      mesh2d.internalVertices().size(), 0);
  for (auto &faceProxy : mesh_proxy.faces()) {
    auto fh_proxy = faceProxy.topo.hd;
    int ent = fhProxy2ent[fh_proxy];
    const Plane3 &plane = entities[ent]->supportingPlane.reconstructed;
    for (auto hh_proxy : faceProxy.topo.halfedges) {
      VertHandle vhProxy = mesh_proxy.topo(hh_proxy).to();
      VertHandle vh2d = mesh_proxy.data(vhProxy);
      Vec3 direction = normalize(cur_cam.direction(mesh2d.data(vh2d)));
      Point3 position = Intersection(Ray3(cur_cam.eye(), direction), plane);
      curReconstruction.data(vh2d) += position;
      vh2d2faceProxyCount[vh2d]++;
    }
  }
  for (auto &vert : curReconstruction.vertices()) {
    vert.data /= double(vh2d2faceProxyCount[vert.topo.hd]);
  }

  if (false) { // show current reconstruction
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XTriangles;
    sb.installingOptions().discretizeOptions.color(gui::Black);
    sb.installingOptions().lineWidth = 10;
    AddToScene(sb, curReconstruction, [](auto) { return true; },
               [](const Point3 &pos) { return pos; },
               [](HalfHandle hh) { return gui::Black; },
               [](FaceHandle fh) {
                 return gui::Color(
                     Vec3i(rand() % 256, rand() % 256, rand() % 256));
               });
    sb.show(true, true, gui::RenderOptions()
                            .backgroundColor(gui::White)
                            .renderMode(gui::All)
                            .bwTexColor(0.0)
                            .bwColor(1.0)
                            .fixUpDirectionInCameraMove(false)
                            .cullBackFace(false)
                            .cullFrontFace(false));
  }
  return curReconstruction;
}

// OptimizeWithoutOrientations
void OptimizeWithoutOrientations(
    Mesh3 &cur_reconstruction, const PerspectiveCamera &cur_cam,
    const Mesh<VertHandle, HalfHandle, FaceHandle> &mesh_proxy,
    const std::vector<SubMesh> &submeshes,
    const std::unordered_set<VertHandle> &non_corner_vh_proxies) {
  using namespace Eigen;

  if(false){ // show together
    auto color_table = gui::CreateRandomColorTableWithSize(submeshes.size());
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XLines;
    sb.installingOptions().lineWidth = 10;
    for (int i = 0; i < submeshes.size(); i++) {
      sb.installingOptions().discretizeOptions.color(color_table[i]);
      for (HalfHandle hh : submeshes[i].hhs) {
        VertHandle vh1 = mesh_proxy.topo(hh).from();
        VertHandle vh2 = mesh_proxy.topo(hh).to();
        Point3 p1 = cur_reconstruction.data(mesh_proxy.data(vh1));
        Point3 p2 = cur_reconstruction.data(mesh_proxy.data(vh2));
        sb.add(Line3(p1, p2));
      }
    }
    sb.show(true, true, gui::RenderOptions()
                            .winName("Together Before Optimization")
                            .fixUpDirectionInCameraMove(false));
  }

  EnergyWeights weights;
  weights.face_angle_weight = 1e2;
  weights.face_msda_weight = 1e1;
  weights.vert_msda_weight = 1e3;
  weights.all_msda_weight = 1e5;

  // holistic optimization
  {
    std::vector<FaceHandle> fhs;
    for (auto &f : cur_reconstruction.faces()) {
      fhs.push_back(f.topo.hd);
    }
    // get fundamental vhs and sort fhs
    auto fundamental_vhs_set = SortFacesWithPlanarityDependency(
        cur_reconstruction, fhs.begin(), fhs.end(),
        [&cur_reconstruction](auto vhs_begin, auto vhs_end) -> bool {
          auto trans_fun =
              [&cur_reconstruction](VertHandle vh) -> const Point3 & {
            return cur_reconstruction.data(vh);
          };
          return IsFuzzyColinear(MakeTransformIterator(vhs_begin, trans_fun),
                                 MakeTransformIterator(vhs_end, trans_fun),
                                 1e-6);
        });

    std::vector<VertHandle> fundamental_vhs(fundamental_vhs_set.begin(),
                                            fundamental_vhs_set.end());

	// get v2matrix and f2matrix
    std::map<VertHandle, DenseMatd> v2matrix;
    std::map<FaceHandle, DenseMatd> f2matrix;
    ConstructSingleViewTransformMatrices(
        fundamental_vhs, fhs, // face2related_verts
        [&cur_reconstruction](FaceHandle fh) -> std::vector<VertHandle> {
          std::vector<VertHandle> related_vhs;
          for (HalfHandle hh : cur_reconstruction.topo(fh).halfedges) {
            related_vhs.push_back(cur_reconstruction.topo(hh).from());
          }
          return related_vhs;
        }, // vert2direction
        [&cur_reconstruction, &cur_cam](VertHandle vh) -> Vec3 {
          Vec3 dir = cur_reconstruction.data(vh);
          return normalize(dir);
        },
        v2matrix, f2matrix);

    DenseMatd inversed_depths_vector;
    {
      inversed_depths_vector = EstimateSingleViewInversedDepths(
          v2matrix, f2matrix,
          // vert2initial_depth
          [&cur_reconstruction, &cur_cam](VertHandle vh) -> double {
            // return 1.0 / vh2inversed_depth[vh];
            return Distance(cur_reconstruction.data(vh), cur_cam.eye());
          },
          // energy_fun
          [&cur_reconstruction, &mesh_proxy, &cur_cam, &submeshes, &weights,
           &non_corner_vh_proxies](auto &&vh2depth,
                                   auto &&fh2eq) -> std::vector<double> {
            std::vector<double> energy_terms;
            for (auto &submesh : submeshes) {
              auto energy_terms_submesh = ComputeEnergy(
                  mesh_proxy, submesh.fhs.begin(), submesh.fhs.end(), weights,
                  [&cur_reconstruction, &mesh_proxy,
                   &vh2depth](VertHandle vh_proxy) -> Point3 {
                    VertHandle vh = mesh_proxy.data(vh_proxy);
                    return normalize(cur_reconstruction.data(vh)) *
                           vh2depth(vh);
                  },
                  [&cur_reconstruction, &mesh_proxy, &vh2depth, &fh2eq,
                   &non_corner_vh_proxies](FaceHandle fh_proxy) -> Vec3 {
                    FaceHandle fh = mesh_proxy.data(fh_proxy);
                    if (fh.valid()) {
                      return fh2eq(fh);
                    } else {
                      auto &hh_proxies = mesh_proxy.topo(fh_proxy).halfedges;
                      assert(hh_proxies.size() >= 3);
                      for (int i1 = 0; i1 < hh_proxies.size(); i1++) {
                        VertHandle vh1 = mesh_proxy.data(
                            mesh_proxy.topo(hh_proxies[i1]).from());
                        Point3 point1 =
                            normalize(cur_reconstruction.data(vh1)) *
                            vh2depth(vh1);
                        for (int i2 = i1 + 1; i2 < hh_proxies.size(); i2++) {
                          VertHandle vh2 = mesh_proxy.data(
                              mesh_proxy.topo(hh_proxies[i2]).from());
                          Point3 point2 =
                              normalize(cur_reconstruction.data(vh2)) *
                              vh2depth(vh2);
                          for (int i3 = i2 + 1; i3 < hh_proxies.size(); i3++) {
                            VertHandle vh3 = mesh_proxy.data(
                                mesh_proxy.topo(hh_proxies[i3]).from());
                            Point3 point3 =
                                normalize(cur_reconstruction.data(vh3)) *
                                vh2depth(vh3);
                            Vec3 eq = Plane3ToEquation(
                                Plane3From3Points(point1, point2, point3));
                            if (norm(eq) > 1e-5 && !HasValue(eq, [](auto &v) {
                                  return IsInfOrNaN(v);
                                })) {
                              return eq;
                            }
                          }
                        }
                      }
                      ASSERT_OR_PANIC(false && "degenerated face!");
                      return Vec3();
                    }
                  },
                  [&non_corner_vh_proxies](VertHandle vh_proxy) -> bool {
                    return Contains(non_corner_vh_proxies, vh_proxy);
                  },
                  cur_cam.eye());
              energy_terms.insert(energy_terms.end(),
                                  energy_terms_submesh.begin(),
                                  energy_terms_submesh.end());
            }
            return energy_terms;
          });

      assert(std::all_of(inversed_depths_vector.begin(),
                         inversed_depths_vector.end(),
                         [](double e) { return e > 0; }));
    }

    if (true) {
      gui::SceneBuilder sb;
      sb.installingOptions().defaultShaderSource =
          gui::OpenGLShaderSourceDescriptor::XTriangles;
      sb.installingOptions().discretizeOptions.color(gui::Black);
      sb.installingOptions().lineWidth = 10;

      for (auto &h : cur_reconstruction.halfedges()) {
        HalfHandle hh = h.topo.hd;
        VertHandle vh1 = cur_reconstruction.topo(hh).from();
        VertHandle vh2 = cur_reconstruction.topo(hh).to();
        Point3 p1 = normalize(cur_reconstruction.data(vh1)) /
                    DenseMatd(v2matrix.at(vh1) * inversed_depths_vector)(0, 0);
        Point3 p2 = normalize(cur_reconstruction.data(vh2)) /
                    DenseMatd(v2matrix.at(vh2) * inversed_depths_vector)(0, 0);
        sb.add(Line3(p1, p2));
      }
      sb.show(true, true, gui::RenderOptions()
                               .camera(cur_cam)
                               .winName("After Holistic Optimization")
                               .fixUpDirectionInCameraMove(false));
    }
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
  auto mesh = LoadFromObjFile(obj_file_name);

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

  // check validity
  for (auto &h : mesh_proxy.halfedges()) {
    HalfHandle hh = h.data;
    if (hh.invalid()) {
      hh = mesh_proxy.data(h.topo.opposite);
    }
    assert(hh.valid());
  }
  std::unordered_map<FaceHandle, FaceHandle> cut_face2another;
  for (auto &cut_face_pair : cut_face_pairs) {
    cut_face2another[cut_face_pair.first] = cut_face_pair.second;
    cut_face2another[cut_face_pair.second] = cut_face_pair.first;
  }
  auto submeshes = ExtractSubMeshes(mesh_proxy);
  Println("found ", submeshes.size(), " submeshes");
  HandledTable<FaceHandle, int> fh_proxy2submesh_id(
      mesh_proxy.internalFaces().size(), -1);
  for (int i = 0; i < submeshes.size(); i++) {
    for (auto fh : submeshes[i].fhs) {
      fh_proxy2submesh_id[fh] = i;
    }
  }

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

  if (true) { // show line classification results
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
  std::vector<FaceVertDependency> dependencies(submeshes.size());
  for (int submesh_id = 0; submesh_id < submeshes.size(); submesh_id++) {
    // find the vhs that should serve as anchoring variables
    // and sort the order of the fhs to be determined
    dependencies[submesh_id].ordered_fhs = std::vector<FaceHandle>(
        submeshes[submesh_id].fhs.begin(), submeshes[submesh_id].fhs.end());
    dependencies[submesh_id].fundamental_vhs = SortFacesWithPlanarityDependency(
        mesh_proxy, dependencies[submesh_id].ordered_fhs.begin(),
        dependencies[submesh_id].ordered_fhs.end(),
        [&mesh_proxy, &mesh2d](auto vhs_begin, auto vhs_end) -> bool {
          // whether these vhs are colinear?
          auto trans_fun = [&mesh2d,
                            &mesh_proxy](VertHandle vh) -> const Point2 & {
            return mesh2d.data(mesh_proxy.data(vh));
          };
          return IsFuzzyColinear(MakeTransformIterator(vhs_begin, trans_fun),
                                 MakeTransformIterator(vhs_end, trans_fun),
                                 1e-6);
        });
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
    std::map<VertHandle, DenseMatd> vh_proxy2matrix;
    std::map<FaceHandle, DenseMatd> fh_proxy2matrix;
    // TODO: construct matrices

    // TODO: call PerspectiveReconstruction
    // 

    Mesh3 mesh_reconstructed;
    static const bool no_orientation_estimation = false;
    if (!no_orientation_estimation) {
      mesh_reconstructed =
          ReconstructWithOrientations(edge2line, edge2vp, vp2dir, hh2d2edge,
                                      cur_cam, mesh_proxy, mesh2d, matlab);
    } else {
      mesh_reconstructed =
          Transform(mesh2d, [&cur_cam](const Point2 &p) -> Point3 {
            return normalize(cur_cam.direction(p));
          });
    }
    OptimizeWithoutOrientations(mesh_reconstructed, cur_cam, mesh_proxy,
                                submeshes, non_corner_vh_proxies);

    reconstructions.push_back(std::move(mesh_reconstructed));
  }

  return 0;
}
