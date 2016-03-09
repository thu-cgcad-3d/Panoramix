#include <QtWidgets>

#include "../../src/core/factor_graph.hpp"
#include "../../src/core/parallel.hpp"
#include "../../src/misc/cache.hpp"
#include "../../src/misc/clock.hpp"

#include "../../src/gui/canvas.hpp"
#include "../../src/gui/qttools.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/singleton.hpp"
#include "../../src/gui/utility.hpp"

#include "../../src/experimental/line_drawing.hpp"
#include "../../src/experimental/pi_graph_annotation.hpp"
#include "../../src/experimental/pi_graph_cg.hpp"
#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
#include "../../src/experimental/pi_graph_optimize.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"

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

inline std::vector<Point2> PossibleKeyVanishingPoints(const Chain2 &chain) {
  assert(chain.size() > 2);
  if (chain.size() == 3) {
    return {};
  }
  if (chain.size() == 4) {
    return {Intersection(chain.edge(0).ray(), chain.edge(2).ray()),
            Intersection(chain.edge(1).ray(), chain.edge(3).ray())};
  }
  if (chain.size() % 2 == 0) {
    std::vector<Point2> vpPositions;
    vpPositions.reserve(chain.size() * chain.size() / 4);
    for (int i = 0; i < chain.size() / 2; i++) {
      for (int j = i + 1; j < i + chain.size() / 2; j++) {
        vpPositions.push_back(
            Intersection(Line2(chain.at(i), chain.at(j)).ray(),
                         Line2(chain.at(i + chain.size() / 2),
                               chain.at(j + chain.size() / 2))
                             .ray()));
      }
    }
    return vpPositions;
  } else {
    std::vector<Point2> vpPositions;
    for (int i = 0; i < chain.size(); i++) {
      for (int j = i + 1; j < chain.size(); j++) {
        vpPositions.push_back(
            Intersection(chain.edge(i).ray(), chain.edge(j).ray()));
      }
    }
    return vpPositions;
  }
}

int main(int argc, char **argv) {
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");

  std::string name = "hex";
  std::string camName = "cam1";
  bool resetCam = false;

  std::string objFile = "H:\\GitHub\\Panoramix\\data\\linedrawing\\" + name +
                        "\\" + name + ".obj";
  std::string camFile = "H:\\GitHub\\Panoramix\\data\\linedrawing\\" + name +
                        "\\" + name + ".obj." + camName + ".cereal";

  auto mesh = LoadFromObjFile(objFile);

  auto meshProxy = MakeMeshProxy(mesh);

  // decompose
  auto cutFacePairs = DecomposeAll(
      meshProxy, [](HalfHandle hh1, HalfHandle hh2) -> bool { return false; });
  std::unordered_map<FaceHandle, FaceHandle> cutFace2Another;
  for (auto &cutFacePair : cutFacePairs) {
    cutFace2Another[cutFacePair.first] = cutFacePair.second;
    cutFace2Another[cutFacePair.second] = cutFacePair.first;
  }

  auto subMeshes = ExtractSubMeshes(meshProxy,
                                    [](auto hhbegin, auto hhend) -> bool {
                                      return std::distance(hhbegin, hhend) <= 1;
                                    },
                                    10);
  Println("found ", subMeshes.size(), " subMeshes");

  PerspectiveCamera cam;
  if (!LoadFromDisk(camFile, cam) || resetCam) {
    auto sphere = BoundingBoxOfContainer(mesh.vertices()).outerSphere();
    PerspectiveCamera projCam(500, 500, Point2(250, 250), 200,
                              sphere.center + Vec3(1, 2, 3) * sphere.radius * 2,
                              sphere.center);

    const gui::ColorTable ctable =
        gui::CreateRandomColorTableWithSize(subMeshes.size());

    HandledTable<HalfHandle, int> hh2subMeshId(
        meshProxy.internalHalfEdges().size(), -1);
    for (auto &h : meshProxy.halfedges()) {
      int &id = hh2subMeshId[h.topo.hd];
      for (int i = 0; i < subMeshes.size(); i++) {
        if (subMeshes[i].contains(h.topo.hd)) {
          id = i;
          break;
        }
      }
    }
    HandledTable<FaceHandle, int> fh2subMeshId(meshProxy.internalFaces().size(),
                                               -1);
    for (auto &f : meshProxy.faces()) {
      int &id = fh2subMeshId[f.topo.hd];
      for (int i = 0; i < subMeshes.size(); i++) {
        if (subMeshes[i].contains(f.topo.hd)) {
          id = i;
          break;
        }
      }
    }

    // show each subMesh
    if (true) {
      for (int i = 0; i < subMeshes.size(); i++) {
        Println("subMesh - ", i);
        gui::SceneBuilder sb;
        sb.installingOptions().defaultShaderSource =
            gui::OpenGLShaderSourceDescriptor::XTriangles;
        sb.installingOptions().discretizeOptions.color(gui::Black);
        sb.installingOptions().lineWidth = 0.03;
        AddToScene(
            sb, meshProxy,
            [&subMeshes, i](auto h) { return subMeshes[i].contains(h); },
            [&mesh](VertHandle vh) { return mesh.data(vh); },
            [&hh2subMeshId, &ctable](HalfHandle hh) { return gui::Black; },
            [&fh2subMeshId, &ctable](FaceHandle fh) {
              return ctable[fh2subMeshId[fh]];
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
      AddToScene(sb, meshProxy, [](auto) { return true; },
                 [&mesh](VertHandle vh) { return mesh.data(vh); },
                 [&hh2subMeshId, &ctable](HalfHandle hh) { return gui::Black; },
                 [&fh2subMeshId, &ctable](FaceHandle fh) {
                   return ctable[fh2subMeshId[fh]];
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
    SaveToDisk(camFile, cam);
  }

  // convert to 2d
  auto mesh2d = Transform(
      mesh, [&cam](const Point3 &p) -> Point2 { return cam.toScreen(p); });

  // add offset noise
  Vec2 offsetNoise = Vec2(20, -20);
  for (auto &v : mesh2d.vertices()) {
    v.data += offsetNoise;
  }

  {
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

  auto point2dAt = [&mesh2d, &meshProxy](VertHandle vhInProxy) -> Point2 {
    return mesh2d.data(meshProxy.data(vhInProxy));
  };
  auto line2dAt = [&mesh2d, &meshProxy,
                   point2dAt](HalfHandle hhInProxy) -> Line2 {
    return Line2(point2dAt(meshProxy.topo(hhInProxy).from()),
                 point2dAt(meshProxy.topo(hhInProxy).to()));
  };

  Box2 box = BoundingBoxOfContainer(mesh2d.vertices());
  double scale = box.outerSphere().radius;

  struct PPFocalCandidate {
    Point2 pp;
    double focal;
  };

  // collect pp focal candidates
  std::vector<PPFocalCandidate> ppFocalCandidates;
  ppFocalCandidates.reserve(subMeshes.size() * 3);

  for (int subMeshId = 0; subMeshId < subMeshes.size(); subMeshId++) {
    // collect edge intersections in each face
    std::vector<Point2> interps;
    for (auto fh : subMeshes[subMeshId].fhs) {
      auto &hhs = meshProxy.topo(fh).halfedges;
      Chain2 corners;
      for (auto hh : hhs) {
        corners.append(point2dAt(meshProxy.topo(hh).to()));
      }
      auto keyVPs = PossibleKeyVanishingPoints(corners);
      interps.insert(interps.end(), keyVPs.begin(), keyVPs.end());
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
              Distance(pp, box.center()) > scale * 2.0) {
            continue;
          }
          ppFocalCandidates.push_back(PPFocalCandidate{pp, focal});
        }
      }
    }
  }

  std::sort(ppFocalCandidates.begin(), ppFocalCandidates.end(),
            [](auto &a, auto &b) { return a.focal < b.focal; });

  // naive clustering
  std::vector<std::pair<std::set<int>, PPFocalCandidate>> ppFocalGroups;
  {
    std::vector<int> ppFocalId2group(ppFocalCandidates.size(), -1);
    int ngroups = 0;
    RTreeMap<Vec3, int> ppFocalIdTree;
    for (int i = 0; i < ppFocalCandidates.size(); i++) {
      Vec3 coordinate =
          cat(ppFocalCandidates[i].pp, ppFocalCandidates[i].focal);
      const double thres = scale / 50.0;
      // find the nearest ppFocal sample point
      int nearestPPFocalCandId = -1;
      double minDist = thres;
      ppFocalIdTree.search(BoundingBox(coordinate).expand(thres * 2),
                           [&nearestPPFocalCandId, &minDist,
                            &coordinate](const std::pair<Vec3, int> &cand) {
                             double dist = Distance(cand.first, coordinate);
                             if (dist < minDist) {
                               minDist = dist;
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

    ppFocalGroups.resize(ngroups);
    for (auto &g : ppFocalGroups) {
      g.second.focal = 0.0;
      g.second.pp = Point2();
    }
    for (int i = 0; i < ppFocalId2group.size(); i++) {
      auto &g = ppFocalGroups[ppFocalId2group[i]];
      g.first.insert(i);
      g.second.focal += ppFocalCandidates[i].focal;
      g.second.pp += ppFocalCandidates[i].pp;
    }
    for (auto &g : ppFocalGroups) {
      g.second.focal /= g.first.size();
      g.second.pp /= double(g.first.size());
    }

    std::sort(
        ppFocalGroups.begin(), ppFocalGroups.end(),
        [](auto &g1, auto &g2) { return g1.first.size() > g2.first.size(); });
  }

  // record edges
  std::vector<std::pair<HalfHandle, HalfHandle>> edge2hhs;
  std::vector<Line2> edge2line;
  HandledTable<HalfHandle, int> hh2edge(mesh2d.internalHalfEdges().size(), -1);
  int nedges = 0;
  {
    for (auto &h : mesh2d.halfedges()) {
      auto hh = h.topo.hd;
      auto oppohh = h.topo.opposite;
      if (hh2edge[hh] == -1 && hh2edge[oppohh] == -1) {
        hh2edge[hh] = hh2edge[oppohh] = nedges;
        nedges++;
        edge2hhs.push_back(MakeOrderedPair(hh, oppohh));
        edge2line.push_back(Line2(mesh2d.data(mesh2d.topo(hh).from()),
                                  mesh2d.data(mesh2d.topo(hh).to())));
      }
    }
    assert(edge2hhs.size() == nedges && edge2line.size() == nedges);
  }

  // collect edge intersections
  std::vector<Point2> intersections;
  std::vector<std::pair<int, int>> intersection2edges;
  {
    intersections.reserve(nedges * (nedges - 1) / 2);
    intersection2edges.reserve(nedges * (nedges - 1) / 2);
    for (int i = 0; i < nedges; i++) {
      const Line2 &linei = edge2line[i];
      for (int j = i + 1; j < nedges; j++) {
        const Line2 &linej = edge2line[j];
        Point2 interp = Intersection(linei.ray(), linej.ray());
        if (std::min(Distance(interp, linei), Distance(interp, linej)) <=
            scale / 10.0) {
          continue;
        }
        intersections.push_back(interp);
        intersection2edges.emplace_back(i, j);
      }
    }
    assert(intersections.size() == intersection2edges.size());
  }

  // get vpPositions from the intersections
  std::vector<Point2> vpPositions;
  int nvps = 0;
  {
    std::vector<int> intersection2vp(intersections.size(), -1);
    RTreeMap<Point2, int> intersectionTree;
    for (int i = 0; i < intersections.size(); i++) {
      const double thres = scale / 50.0;
      auto &p = intersections[i];
      int nearestIntersectionId = -1;
      double minDist = thres;
      intersectionTree.search(
          BoundingBox(p).expand(thres * 2),
          [&nearestIntersectionId, &minDist,
           &p](const std::pair<Point2, int> &locationAndIntersectionId) {
            double dist = Distance(locationAndIntersectionId.first, p);
            if (dist < minDist) {
              minDist = dist;
              nearestIntersectionId = locationAndIntersectionId.second;
            }
            return true;
          });
      if (nearestIntersectionId != -1) {
        intersection2vp[i] = intersection2vp[nearestIntersectionId];
      } else {
        intersection2vp[i] = nvps++;
      }
      intersectionTree.emplace(p, i);
    }

    vpPositions.resize(nvps, Origin<2>());
    std::vector<std::set<int>> vp2intersections(nvps);
    for (int i = 0; i < intersection2vp.size(); i++) {
      int vpid = intersection2vp[i];
      vpPositions[vpid] += intersections[i];
      vp2intersections[vpid].insert(i);
    }
    for (int i = 0; i < nvps; i++) {
      vpPositions[i] /= double(vp2intersections[i].size());
    }
  }

  // regroup edges to vpPositions
  struct EVPBinding {
    int edgeId;
    int vpId;
    double angle;
  };
  std::vector<EVPBinding> evps;
  std::vector<std::set<int>> vp2evps(nvps);
  std::vector<std::vector<int>> edge2evps(nedges);

  {
    evps.reserve(nvps * 2);
    for (int vp = 0; vp < nvps; vp++) {
      auto &vpPos = vpPositions[vp];
      for (int edge = 0; edge < nedges; edge++) {
        auto &line = edge2line[edge];
        double angle =
            AngleBetweenUndirected(line.direction(), vpPos - line.center());
        static const double theta = DegreesToRadians(3);
        if (angle >= theta) {
          continue;
        }
        evps.push_back(EVPBinding{edge, vp, angle});
        int evp = evps.size() - 1;
        vp2evps[vp].insert(evp);
        edge2evps[edge].push_back(evp);
      }
    }
    std::vector<int> orderedVPIds(nvps);
    std::iota(orderedVPIds.begin(), orderedVPIds.end(), 0);
    std::sort(orderedVPIds.begin(), orderedVPIds.end(),
              [&vp2evps](int vpid1, int vpid2) {
                return vp2evps[vpid1].size() > vp2evps[vpid2].size();
              });
    std::vector<Point2> orderedVPs(nvps);
    std::vector<std::set<int>> orderedVP2bindings(nvps);
    for (int i = 0; i < nvps; i++) {
      orderedVPs[i] = vpPositions[orderedVPIds[i]];
      orderedVP2bindings[i] = std::move(vp2evps[orderedVPIds[i]]);
    }
    for (auto &b : evps) {
      b.vpId = std::find(orderedVPIds.begin(), orderedVPIds.end(), b.vpId) -
               orderedVPIds.begin();
    }
    vpPositions = std::move(orderedVPs);
    vp2evps = std::move(orderedVP2bindings);
  }

  //std::vector<std::set<int>> vp2edges(nvps);
  //for (int i = 0; i < nvps; i++) {
  //  for (int b : vp2evps[i]) {
  //    vp2edges[i].insert(evps[b].edgeId);
  //  }
  //}
  //std::vector<double> vp2weights(nvps, 0.0);
  //for (int i = 0; i < nvps; i++) {
  //  for (int edge : vp2edges[i]) {
  //    vp2weights[i] += edge2line[edge].length();
  //  }
  //}

  // construct a factor graph to optimize edge-vp bindings
  FactorGraph fg;
  {
    std::vector<FactorGraph::VarHandle> edge2vh(nedges);
    for (int edge = 0; edge < nedges; edge++) {
      auto vc = fg.addVarCategory(edge2evps[edge].size() + 1, 1.0);
      edge2vh[edge] = fg.addVar(vc);
    }

    // potential 1: the edge should bind to some vp (?)
    for (int edge = 0; edge < nedges; edge++) {
      auto vh = edge2vh[edge];
      auto &relatedEVPs = edge2evps[edge];
      auto fc = fg.addFactorCategory(
          [&evps, &relatedEVPs](const int *varlabels, size_t nvar,
                                FactorGraph::FactorCategoryId fcid,
                                void *givenData) -> double {
            assert(nvar == 1);
            int label = varlabels[0];
            assert(label <= relatedEVPs.size());
            if (label == relatedEVPs.size()) { // not bind to any vp
              return 1.0;
            }
            auto &bd = evps[relatedEVPs[label]];
            // todo
            return 1.0 - Gaussian(bd.angle, DegreesToRadians(3));
          },
          1.0);
      fg.addFactor({vh}, fc);
    }

    // potential 2: two adjacent edges should not bind to a same vp
    for (auto &f : mesh2d.faces()) {
      auto &hhs = f.topo.halfedges;
      for (int i = 0; i < hhs.size(); i++) {
        auto hh1 = hhs[i];
        auto hh2 = hhs[(i + 1) % hhs.size()];
        int edge1 = hh2edge[hh1];
        int edge2 = hh2edge[hh2];
        auto vh1 = edge2vh[edge1];
        auto vh2 = edge2vh[edge2];
        auto fc = fg.addFactorCategory(
            [&evps, &edge2evps, edge1, edge2](
                const int *varlabels, size_t nvar,
                FactorGraph::FactorCategoryId fcid, void *givenData) -> double {
              assert(nvar == 2);
              int bindedVP1 = varlabels[0] == edge2evps[edge1].size()
                                  ? -1
                                  : (evps[edge2evps[edge1][varlabels[0]]].vpId);
              int bindedVP2 = varlabels[1] == edge2evps[edge2].size()
                                  ? -1
                                  : (evps[edge2evps[edge2][varlabels[1]]].vpId);
              // todo
              if (bindedVP1 == -1 || bindedVP2 == -1) {
                return 0;
              }
              if (bindedVP1 == bindedVP2) {
                return 10.0;
              }
              return 0.0;
            },
            1.0);
        fg.addFactor({vh1, vh2}, fc);
      }
    }

    // potential 3: the vpPositions of edges sharing a same face should lie on
    // the same
    // line (the vanishing line of the face)
    for (auto &f : mesh2d.faces()) {
      auto &hhs = f.topo.halfedges;
      std::vector<FactorGraph::VarHandle> vhs;
      vhs.reserve(hhs.size());
      for (auto hh : hhs) {
        vhs.push_back(edge2vh[hh2edge[hh]]);
      }
      auto fc = fg.addFactorCategory(
          [&evps, &edge2evps, &hhs, &vpPositions](
              const int *varlabels, size_t nvar,
              FactorGraph::FactorCategoryId fcid, void *givenData) -> double {
            // todo

          },
          1.0);
      fg.addFactor(vhs.begin(), vhs.end(), fc);
    }

    // potential 4: a vp is good only when edges are bound to it!
    // a vp should be with either no edges or >= 3 edges
    for (int vpid = 0; vpid < nvps; vpid++) {
      std::vector<FactorGraph::VarHandle> relatedVhs;
      relatedVhs.reserve(vp2evps[vpid].size());
      for (int bind : vp2evps[vpid]) {
        int edge = evps[bind].edgeId;
        relatedVhs.push_back(edge2vh[edge]);
      }
      auto fc = fg.addFactorCategory(
          [&edge2evps, &evps, vpid, &vp2evps](
              const int *varlabels, size_t nvar,
              FactorGraph::FactorCategoryId fcid, void *givenData) -> double {
            int bindedEdges = 0;
            auto &bindings = vp2evps[vpid];
            int i = 0;
            assert(nvar == bindings.size());
            for (auto it = bindings.begin(); it != bindings.end(); ++it, ++i) {
              int bind = *it;
              int edge = evps[bind].edgeId;
              int edgeLabel = varlabels[i];
              if (edgeLabel == edge2evps[edge].size()) {
                continue; // not bound to any vp
              }
              int edgeBindedVPId = evps[edge2evps[edge][edgeLabel]].vpId;
              if (edgeBindedVPId == vpid) {
                bindedEdges++;
              }
            }
            // todo
            if (bindedEdges == 1 || bindedEdges == 2) {
              return 10.0;
            }
            return 0.0;
          },
          1.0);
      fg.addFactor(relatedVhs.begin(), relatedVhs.end(), fc);
    }
  }

  if (true) {
    for (int i = 0; i < nvps; i++) {
      Image3ub im(cam.screenSize(), Vec3ub(255, 255, 255));
      auto canvas = gui::MakeCanvas(im);
      canvas.color(gui::LightGray);
      canvas.thickness(2);
      for (auto &line : edge2line) {
        canvas.add(line);
      }
      canvas.color(gui::Gray);
      canvas.thickness(1);
      for (int binding : vp2evps[i]) {
        int edge = evps[binding].edgeId;
        canvas.add(edge2line[edge].ray());
      }
      canvas.color(gui::Black);
      for (int binding : vp2evps[i]) {
        int edge = evps[binding].edgeId;
        canvas.add(edge2line[edge]);
      }
      canvas.show(0, "vp_" + std::to_string(i));
    }
  }

  for (int configId = 0;
       configId < std::min(5ull, ppFocalGroups.size()) &&
       ppFocalGroups[configId].first.size() * 10 >= ppFocalCandidates.size();
       configId++) {

    double focal = ppFocalGroups[configId].second.focal;
    auto &pp = ppFocalGroups[configId].second.pp;

    PerspectiveCamera curCam(cam.screenWidth(), cam.screenHeight(), pp, focal);
    std::vector<Vec3> vp2dir(nvps);
    for (int i = 0; i < nvps; i++) {
      vp2dir[i] = curCam.direction(vpPositions[i]);
    }
  }
}