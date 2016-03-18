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
  misc::Matlab matlab;

  std::string name = "hex";
  std::string camName = "cam1";
  bool resetCam = false;

  std::string objFile = "H:\\GitHub\\Panoramix\\data\\linedrawing\\" + name +
                        "\\" + name + ".obj";
  std::string camFile = "H:\\GitHub\\Panoramix\\data\\linedrawing\\" + name +
                        "\\" + name + ".obj." + camName + ".cereal";

  //// [Load Mesh]
  auto mesh = LoadFromObjFile(objFile);
  auto meshProxy = MakeMeshProxy(mesh);

  //// [Decompose]
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
  HandledTable<FaceHandle, int> fhProxy2subMeshId(
      meshProxy.internalFaces().size(), -1);
  for (int i = 0; i < subMeshes.size(); i++) {
    for (auto fh : subMeshes[i].fhs) {
      fhProxy2subMeshId[fh] = i;
    }
  }

  //// [Load Camera]
  PerspectiveCamera cam;
  if (!LoadFromDisk(camFile, cam) || resetCam) {
    auto sphere = BoundingBoxOfContainer(mesh.vertices()).outerSphere();
    PerspectiveCamera projCam(500, 500, Point2(250, 250), 200,
                              sphere.center + Vec3(1, 2, 3) * sphere.radius * 2,
                              sphere.center);

    const gui::ColorTable ctable =
        gui::CreateRandomColorTableWithSize(subMeshes.size());

    HandledTable<HalfHandle, int> hhProxy2subMeshId(
        meshProxy.internalHalfEdges().size(), -1);
    for (auto &h : meshProxy.halfedges()) {
      int &id = hhProxy2subMeshId[h.topo.hd];
      for (int i = 0; i < subMeshes.size(); i++) {
        if (subMeshes[i].contains(h.topo.hd)) {
          id = i;
          break;
        }
      }
    }
    HandledTable<FaceHandle, int> fhProxy2subMeshId(
        meshProxy.internalFaces().size(), -1);
    for (auto &f : meshProxy.faces()) {
      int &id = fhProxy2subMeshId[f.topo.hd];
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
            [&hhProxy2subMeshId, &ctable](HalfHandle hh) { return gui::Black; },
            [&fhProxy2subMeshId, &ctable](FaceHandle fh) {
              return ctable[fhProxy2subMeshId[fh]];
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
          sb, meshProxy, [](auto) { return true; },
          [&mesh](VertHandle vh) { return mesh.data(vh); },
          [&hhProxy2subMeshId, &ctable](HalfHandle hh) { return gui::Black; },
          [&fhProxy2subMeshId, &ctable](FaceHandle fh) {
            return ctable[fhProxy2subMeshId[fh]];
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

  //// [Make 2D Mesh]
  // convert to 2d
  auto mesh2d = Transform(
      mesh, [&cam](const Point3 &p) -> Point2 { return cam.toScreen(p); });

  // add offset noise
  Vec2 offsetNoise = Vec2(20, -20);
  for (auto &v : mesh2d.vertices()) {
    v.data += offsetNoise;
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

  //// [Estimate PP & Focal Candidates from 2D Mesh]
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

  //// [Orient Edges]
  // record edges
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

  // collect edge intersections and
  // get vpPositions from the intersections
  std::vector<Point2> vpPositions;
  int nvps = 0;
  std::vector<std::vector<Scored<int>>> edge2OrderedVPAndAngles(nedges);
  {
    std::vector<Point2> intersections;
    std::vector<std::pair<int, int>> intersection2edges;
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

    std::vector<int> intersection2Rawvp(intersections.size(), -1);
    RTreeMap<Point2, int> intersectionTree;
    for (int i = 0; i < intersections.size(); i++) {
      const double thres = scale / 30.0;
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
        intersection2Rawvp[i] = intersection2Rawvp[nearestIntersectionId];
      } else {
        intersection2Rawvp[i] = nvps++;
      }
      intersectionTree.emplace(p, i);
    }

    // further merge vps
    // if any two vps share two or more edges, merge them
    std::vector<std::set<int>> rawVP2edges(nvps);
    for (int i = 0; i < intersection2Rawvp.size(); i++) {
      int vpid = intersection2Rawvp[i];
      rawVP2edges[vpid].insert(i);
    }
    std::vector<std::set<int>> rawVp2rawVpShouldMerge(nvps);
    for (int vp1 = 0; vp1 < nvps; vp1++) {
      auto &edges1 = rawVP2edges[vp1];
      for (int vp2 = vp1 + 1; vp2 < nvps; vp2++) {
        auto &edges2 = rawVP2edges[vp2];
        std::set<int> commonEdges;
        std::set_intersection(edges1.begin(), edges1.end(), edges2.begin(),
                              edges2.end(),
                              std::inserter(commonEdges, commonEdges.begin()));
        if (commonEdges.size() >= 2) {
          rawVp2rawVpShouldMerge[vp1].insert(vp2);
          rawVp2rawVpShouldMerge[vp2].insert(vp1);
        }
      }
    }
    std::map<int, std::set<int>> newVP2rawVPs;
    std::vector<int> rawVP2newVP(nvps, -1);
    std::vector<int> rawVPIds(nvps);
    std::iota(rawVPIds.begin(), rawVPIds.end(), 0);
    nvps = ConnectedComponents(
        rawVPIds.begin(), rawVPIds.end(),
        [&rawVp2rawVpShouldMerge](int vp) -> const std::set<int> & {
          return rawVp2rawVpShouldMerge.at(vp);
        },
        [&newVP2rawVPs, &rawVP2newVP](int rawVP, int newVP) {
          newVP2rawVPs[newVP].insert(rawVP);
          rawVP2newVP[rawVP] = newVP;
        });

    vpPositions.resize(nvps, Origin<2>());
    std::vector<std::set<int>> vp2intersections(nvps);
    for (int i = 0; i < intersection2Rawvp.size(); i++) {
      int rawVP = intersection2Rawvp[i];
      int newVP = rawVP2newVP[rawVP];
      vpPositions[newVP] += intersections[i]; // TODO: what if some
                                              // intersections are oppsite far
                                              // points?
      vp2intersections[newVP].insert(i);
    }
    for (int i = 0; i < nvps; i++) {
      vpPositions[i] /= double(vp2intersections[i].size());
    }

    // initial edge vp bindings
    std::vector<std::map<int, double>> vp2edgeWithAngles(nvps);
    std::vector<bool> vpIsGood(nvps, true);
    for (int vp = 0; vp < nvps; vp++) {
      auto &vpPos = vpPositions[vp];
      for (int edge = 0; edge < nedges; edge++) {
        auto &line = edge2line[edge];
        double lambda = ProjectionOfPointOnLine(vpPos, line).ratio;
        static const double thres = 0.1;
        if (lambda >= -thres && lambda <= 1.0 + thres) {
          continue;
        }
        double angle =
            AngleBetweenUndirected(line.direction(), vpPos - line.center());
        static const double theta = DegreesToRadians(5); ////// TODO
        if (angle >= theta) {
          continue;
        }
        vp2edgeWithAngles[vp][edge] = angle;
      }
      vpIsGood[vp] = vp2edgeWithAngles[vp].size() >= 3;
    }

    if (false) {
      for (int i = 0; i < std::min(10, nvps); i++) {
        Image3ub im(cam.screenSize(), Vec3ub(255, 255, 255));
        auto canvas = gui::MakeCanvas(im);
        canvas.color(gui::LightGray);
        canvas.thickness(2);
        for (auto &line : edge2line) {
          canvas.add(line);
        }
        canvas.color(gui::Gray);
        canvas.thickness(1);
        for (auto &edgeWithAngle : vp2edgeWithAngles[i]) {
          canvas.add(edge2line[edgeWithAngle.first].ray());
        }
        canvas.color(gui::Black);
        for (auto &edgeWithAngle : vp2edgeWithAngles[i]) {
          canvas.add(edge2line[edgeWithAngle.first]);
        }
        canvas.show(0, "before removing bad vps: raw vp_" + std::to_string(i));
      }
    }

    // remove bad vps
    int newvp = 0;
    for (int oldvp = 0; oldvp < nvps; oldvp++) {
      if (vpIsGood[oldvp]) {
        // update vpPositions
        // update vp2edgeWithAngles
        vpPositions[newvp] = vpPositions[oldvp];
        vp2edgeWithAngles[newvp] = std::move(vp2edgeWithAngles[oldvp]);
        newvp++;
      }
    }
    nvps = newvp;
    vpPositions.resize(nvps);
    vp2edgeWithAngles.resize(nvps);

    // construct edge2OrderedVPs
    for (int vp = 0; vp < nvps; vp++) {
      for (auto &edgeAndAngle : vp2edgeWithAngles[vp]) {
        int edge = edgeAndAngle.first;
        double angle = edgeAndAngle.second;
        edge2OrderedVPAndAngles[edge].push_back(ScoreAs(vp, angle));
      }
    }
    for (auto &vpAndAngles : edge2OrderedVPAndAngles) {
      std::sort(vpAndAngles.begin(), vpAndAngles.end());
    }

    if (false) {
      for (int i = 0; i < std::min(10, nvps); i++) {
        Image3ub im(cam.screenSize(), Vec3ub(255, 255, 255));
        auto canvas = gui::MakeCanvas(im);
        canvas.color(gui::LightGray);
        canvas.thickness(2);
        for (auto &line : edge2line) {
          canvas.add(line);
        }
        canvas.color(gui::Gray);
        canvas.thickness(2);
        for (auto &edgeWithAngle : vp2edgeWithAngles[i]) {
          canvas.add(edge2line[edgeWithAngle.first].ray());
        }
        canvas.color(gui::Black);
        for (auto &edgeWithAngle : vp2edgeWithAngles[i]) {
          canvas.add(edge2line[edgeWithAngle.first]);
        }
        canvas.show(0, "raw vp_" + std::to_string(i));
      }
    }
  }

  std::vector<int> edge2vp;
  std::vector<std::vector<int>> vp2edges;

  if (!misc::LoadCache(objFile + "-" + camName, "edge2vp_vp2edges", edge2vp,
                       vp2edges)) {
    edge2vp.resize(nedges, -1);
    vp2edges.resize(nvps);

    // construct a factor graph to optimize edge-vp bindings
    FactorGraph fg;
    std::vector<FactorGraph::VarHandle> edge2vh(nedges);
    {
      for (int edge = 0; edge < nedges; edge++) {
        auto vc =
            fg.addVarCategory(edge2OrderedVPAndAngles[edge].size() + 1, 1.0);
        edge2vh[edge] = fg.addVar(vc);
      }

      // potential 1: the edge should bind to some vp, should prefer better
      // scored
      for (int edge = 0; edge < nedges; edge++) {
        auto vh = edge2vh[edge];
        auto &relatedVPAndAngles = edge2OrderedVPAndAngles[edge];
        auto fc = fg.addFactorCategory(
            [&relatedVPAndAngles, nedges](const int *varlabels, size_t nvar,
                                          FactorGraph::FactorCategoryId fcid,
                                          void *givenData) -> double {
              assert(nvar == 1);
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
        fg.addFactor({vh}, fc);
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
          auto &relatedVPAndAngles1 = edge2OrderedVPAndAngles[edge1];
          auto &relatedVPAndAngles2 = edge2OrderedVPAndAngles[edge2];
          auto vh1 = edge2vh[edge1];
          auto vh2 = edge2vh[edge2];
          auto fc = fg.addFactorCategory(
              [edge1, edge2, &relatedVPAndAngles1, &relatedVPAndAngles2,
               &vpPositions, ncorners,
               scale](const int *varlabels, size_t nvar,
                      FactorGraph::FactorCategoryId fcid,
                      void *givenData) -> double {
                assert(nvar == 2);
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
                auto &vpPos1 = vpPositions[bindedVP1];
                auto &vpPos2 = vpPositions[bindedVP2];
                const double thres = scale / 10.0;
                const double K = 10.0 / ncorners;
                if (Distance(vpPos1, vpPos2) < thres) { // todo
                  return K;
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
                [&edge2OrderedVPAndAngles, &vpPositions, prevEdge, edge,
                 nextEdge, ntris](const int *varlabels, size_t nvar,
                                  FactorGraph::FactorCategoryId fcid,
                                  void *givenData) -> double {
                  assert(nvar == 3);
                  int vp1 =
                      varlabels[0] == edge2OrderedVPAndAngles[prevEdge].size()
                          ? -1
                          : (edge2OrderedVPAndAngles[prevEdge][varlabels[0]]
                                 .component);
                  if (vp1 == -1) {
                    return 0.0;
                  }
                  int vp2 = varlabels[1] == edge2OrderedVPAndAngles[edge].size()
                                ? -1
                                : (edge2OrderedVPAndAngles[edge][varlabels[1]]
                                       .component);
                  if (vp2 == -1) {
                    return 0.0;
                  }
                  int vp3 =
                      varlabels[2] == edge2OrderedVPAndAngles[nextEdge].size()
                          ? -1
                          : (edge2OrderedVPAndAngles[nextEdge][varlabels[2]]
                                 .component);
                  if (vp3 == -1) {
                    return 0.0;
                  }
                  if (vp1 == vp2 || vp2 == vp3 || vp1 == vp3) {
                    return 0.0;
                  }
                  double angle = AngleBetweenUndirected(
                      vpPositions[vp1] - vpPositions[vp2],
                      vpPositions[vp3] - vpPositions[vp2]);
                  assert(!IsInfOrNaN(angle));
                  const double K = 30.0 / ntris;
                  return (1.0 - Gaussian(angle, DegreesToRadians(10))) * K;
                },
                1.0);
            fg.addFactor({edge2vh[prevEdge], edge2vh[edge], edge2vh[nextEdge]},
                         fc);
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
      //      [&edge2OrderedVPAndAngles, vpid, &vp2edges](
      //          const int *varlabels, size_t nvar,
      //          FactorGraph::FactorCategoryId fcid, void *givenData) -> double
      //          {
      //        int bindedEdges = 0;
      //        auto &relatedEdges = vp2edges[vpid];
      //        assert(nvar == relatedEdges.size());
      //        for (int i = 0; i < relatedEdges.size(); i++) {
      //          int edge = relatedEdges[i];
      //          int edgeLabel = varlabels[i];
      //          if (edgeLabel == edge2OrderedVPAndAngles[edge].size()) {
      //            continue; // not bound to any vp
      //          }
      //          int edgeBindedVPId =
      //          edge2OrderedVPAndAngles[edge][edgeLabel].component;
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
    auto result =
        fg.solve(5, 1, [](int epoch, double energy, double denergy,
                          const FactorGraph::ResultTable &results) -> bool {
          Println("epoch: ", epoch, "  energy: ", energy);
          return true;
        });

    for (int edge = 0; edge < nedges; edge++) {
      int id = result[edge2vh[edge]];
      if (id == edge2OrderedVPAndAngles[edge].size()) {
        continue;
      }
      edge2vp[edge] = edge2OrderedVPAndAngles[edge][id].component;
      vp2edges[edge2vp[edge]].push_back(edge);
    }

    // invalidate the edge bindings for vps who have only 1 or 2 edges
    for (int vp = 0; vp < nvps; vp++) {
      if (vp2edges[vp].size() <= 2) {
        for (int edge : vp2edges[vp]) {
          edge2vp[edge] = -1;
        }
        vp2edges[vp].clear();
      }
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

    misc::SaveCache(objFile + "-" + camName, "edge2vp_vp2edges", edge2vp,
                    vp2edges);
  }

  assert(edge2vp.size() == nedges);
  assert(vp2edges.size() == nvps);

  // for each pp focal candidate
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

#if 0
    //// [Orient Edges As Much As Possible]
    HandledTable<FaceHandle, Vec3> fh2d2normal(mesh2d.internalFaces().size(),
                                               Origin());
    std::vector<Vec3> edge2dir(nedges, Origin());
    for (int edge = 0; edge < nedges; edge++) {
      int vp = edge2vp[edge];
      if (vp != -1) {
        edge2dir[edge] = normalize(vp2dir[vp]);
      }
    }
    while (true) {
      bool moreEdgesAreOriented = false;
      // update face orientations
      for (auto &f : mesh2d.faces()) {
        auto &hhs = f.topo.halfedges;
        std::vector<Vec3> dirs;
        for (auto hh : hhs) {
          int edge = hh2d2edge[hh];
          auto &dir = edge2dir[edge];
          if (dir == Origin()) {
            continue;
          }
          static const double theta = DegreesToRadians(3);
          if (std::all_of(dirs.begin(), dirs.end(), [&dir](const Vec3 &d) {
                return AngleBetweenUndirected(d, dir) > theta;
              })) {
            dirs.push_back(dir);
          }
        }
        if (dirs.size() >= 2) {
          // the orientation of this face can be determined
          Vec3 faceNormal = normalize(dirs[0].cross(dirs[1]));
          fh2d2normal[f.topo.hd] = faceNormal;

          // update other undetermined edges
          for (auto hh : hhs) {
            int edge = hh2d2edge[hh];
            auto &dir = edge2dir[edge];
            if (dir == Origin()) {
              auto &line = edge2line[edge];
              Vec3 normalOfLineProjPlane =
                  normalize(curCam.direction(line.first)
                                .cross(curCam.direction(line.second)));
              Vec3 edgeDir = normalize(faceNormal.cross(normalOfLineProjPlane));
              edge2dir[edge] = edgeDir;
              moreEdgesAreOriented = true;
            }
          }
        }
      }

      if (!moreEdgesAreOriented) {
        break;
      }
    }

    if (true) { // show determined edges
      std::vector<bool> edge2determined(nedges, false);
      for (int edge = 0; edge < nedges; edge++) {
        edge2determined[edge] = edge2dir[edge] != Origin();
      }
      Image3ub im(cam.screenSize(), Vec3ub(255, 255, 255));
      auto canvas = gui::MakeCanvas(im);
      canvas.thickness(2);
      for (int edge = 0; edge < nedges; edge++) {
        canvas.color(edge2determined[edge] ? gui::Black : gui::LightGray);
        canvas.add(edge2line[edge]);
      }
      canvas.show(0, "determined edges");
    }

    // optimize within each subMesh
    std::vector<std::map<VertHandle, double>> subMesh2vhProxyDepths(
        subMeshes.size());
    for (int i = 0; i < subMeshes.size(); i++) {
      auto &subMesh = subMeshes[i];
      if (subMesh.vhs.empty()) {
        continue;
      }
      std::map<VertHandle, double> vhProxy2depth;
      vhProxy2depth[*subMesh.vhs.begin()] = 1.0;

      while (true) {
        bool moreVertexIsAnchored = false; 

        // check each oriented edge for possible extensions
        for (auto &hhProxy : subMesh.hhs) {
          HalfHandle hh2d = meshProxy.data(hhProxy);
          if (hh2d.invalid()) {
            hh2d = meshProxy.data(meshProxy.topo(hhProxy).opposite);
            assert(hh2d.valid());
          }
          int edge = hh2d2edge[hh2d];
          if (edge2dir[edge] == Origin()) {
            continue;
          }
          const Vec3 &edgeDir = edge2dir[edge];
          VertHandle vhProxy1 = meshProxy.topo(hhProxy).from();
          Vec3 dir1 = normalize(
              curCam.direction(mesh2d.data(meshProxy.data(vhProxy1))));
          VertHandle vhProxy2 = meshProxy.topo(hhProxy).to();
          Vec3 dir2 = normalize(
              curCam.direction(mesh2d.data(meshProxy.data(vhProxy2))));
          if (Contains(vhProxy2depth, vhProxy1) &&
              !Contains(vhProxy2depth,
                        vhProxy2)) { // compute vh2's depth based on vh1's
            Point3 p1 = dir1 * vhProxy2depth.at(vhProxy1);
            Ray3 edgeRay(p1, edgeDir);
            Point3 p2 =
                DistanceBetweenTwoLines(edgeRay, Ray3(curCam.eye(), dir2))
                    .second.first;
            vhProxy2depth[vhProxy2] = norm(p2);
            moreVertexIsAnchored = true;
          }
          if (Contains(vhProxy2depth, vhProxy2) &&
              !Contains(vhProxy2depth,
                        vhProxy1)) { // compute vh1's depth based on vh2's
            Point3 p2 = dir2 * vhProxy2depth.at(vhProxy2);
            Ray3 edgeRay(p2, edgeDir);
            Point3 p1 =
                DistanceBetweenTwoLines(edgeRay, Ray3(curCam.eye(), dir1))
                    .second.first;
            vhProxy2depth[vhProxy1] = norm(p1);
            moreVertexIsAnchored = true;
          }
        }

        // check each oriented face for possible extensions
        for (auto &fhProxy : subMesh.fhs) {
          FaceHandle fh2d = meshProxy.data(fhProxy);
          if (fh2d.invalid()) { // TODO ! we should use the planarity of the
                                // cutting face
            continue;
          }
          if (fh2d2normal[fh2d] == Origin()) {
            continue;
          }
          const Vec3 &faceNormal = fh2d2normal[fh2d];
          // find any anchor
          Point3 anchor = Origin();
          auto &hhProxies = meshProxy.topo(fhProxy).halfedges;
          for (auto hhProxy : hhProxies) {
            VertHandle vhProxy = meshProxy.topo(hhProxy).to();
            if (Contains(vhProxy2depth, vhProxy)) {
              Vec3 dir = normalize(
                  curCam.direction(mesh2d.data(meshProxy.data(vhProxy))));
              anchor = dir * vhProxy2depth.at(vhProxy);
              break;
            }
          }
          if (anchor == Origin()) { // not found any determined anchor
            continue;
          }
          Plane3 plane(anchor, faceNormal);
          for (auto hhProxy : hhProxies) {
            VertHandle vhProxy = meshProxy.topo(hhProxy).to();
            if (!Contains(vhProxy2depth, vhProxy)) {
              Vec3 dir = normalize(
                  curCam.direction(mesh2d.data(meshProxy.data(vhProxy))));
              Point3 p = Intersection(Ray3(curCam.eye(), dir), plane);
              vhProxy2depth[vhProxy] = norm(p);
              moreVertexIsAnchored = true;
            }
          }
        }

        if (!moreVertexIsAnchored) {
          break;
        }
      }
    }

#endif

    // [Reconstruct]
    {
      // entity
      using SupportingPlane = PIConstraintGraph::Entity::SupportingPlane;
      struct EntityBase {
          SupportingPlane supportingPlane;
          EntityBase(const SupportingPlane &sp) : supportingPlane(sp) {}
          virtual FaceHandle fhProxy() const { return FaceHandle(); }
          virtual int edge() const { return -1; }
          bool isFace() const { return fhProxy().valid(); }
          bool isEdge() const { return edge() != -1; }
      };
      struct FaceEntity : EntityBase {
        FaceHandle fh;
        FaceEntity(FaceHandle fh, const Vec3 &center)
            : EntityBase(SupportingPlane(SegControl{-1, -1}, center, {})), fh(fh) {}
        virtual FaceHandle fhProxy() const override { return fh; }
      };
      struct EdgeEntity : EntityBase {
        int e;
        EdgeEntity(int e, const Classified<Line3> &line,
                   const std::vector<Vec3> &vp2dir)
            : EntityBase(SupportingPlane(line, vp2dir)), e(e) {}
        virtual int edge() const override { return e; }
      };

      std::vector<std::unique_ptr<EntityBase>> entities;

      std::vector<int> edge2ent(nedges, -1);
      HandledTable<FaceHandle, int> fhProxy2ent(
          meshProxy.internalFaces().size(), -1);

      // install supporting planes from edges and faces
      // from edges
      for (int edge = 0; edge < nedges; edge++) {
        auto &line2 = edge2line[edge];
        Line3 line3(normalize(curCam.direction(line2.first)),
                    normalize(curCam.direction(line2.second)));
        int vp = edge2vp[edge];
        entities.push_back(
            std::make_unique<EdgeEntity>(edge, ClassifyAs(line3, vp), vp2dir));
        int ent = entities.size() - 1;
        edge2ent[edge] = ent;
      }
      // from faces
      HandledTable<FaceHandle, Point2> fhProxy2center2d(
          meshProxy.internalFaces().size());
      for (auto &f : meshProxy.faces()) {
        Point2 center2 = Origin<2>();
        for (auto hh : f.topo.halfedges) {
          center2 += mesh2d.data(meshProxy.data(meshProxy.topo(hh).to()));
        }
        center2 /= double(f.topo.halfedges.size());
        fhProxy2center2d[f.topo.hd] = center2;
        entities.push_back(std::make_unique<FaceEntity>(
            f.topo.hd, normalize(curCam.direction(center2))));
        int ent = entities.size() - 1;
        fhProxy2ent[f.topo.hd] = ent;
      }

      // install connections
      std::map<std::pair<int, int>, std::vector<Vec3>> ents2anchors;
      for (auto &f : meshProxy.faces()) {
        for (auto hhProxy : f.topo.halfedges) {
          HalfHandle hh2d = meshProxy.data(hhProxy);
          if (hh2d.invalid()) {
            hh2d = meshProxy.data(meshProxy.topo(hhProxy).opposite);
          }
          assert(hh2d.valid());
          int edge = hh2d2edge[hh2d];

          auto &line2 = edge2line[edge];
          ents2anchors[std::make_pair(edge2ent[edge], fhProxy2ent[f.topo.hd])] =
              {normalize(curCam.direction(line2.first)),
               normalize(curCam.direction(line2.second))};
        }
      }

      // install face angles constraints
      std::vector<std::pair<int, int>> adjFaceEnts;
      std::vector<int> adjFace2SubMeshId;
      std::vector<bool> adjFaceOverlapInView;
      for (auto &half : meshProxy.halfedges()) {
        FaceHandle fh1 = half.topo.face;
        FaceHandle fh2 = meshProxy.topo(half.topo.opposite).face;

        int subMeshId1 = fhProxy2subMeshId[fh1];
        int subMeshId2 = fhProxy2subMeshId[fh2];
        assert(subMeshId1 != -1 && subMeshId2 != -1 &&
               subMeshId1 == subMeshId2);
        int faceEnt1 = fhProxy2ent[fh1];
        int faceEnt2 = fhProxy2ent[fh2];
        adjFaceEnts.emplace_back(faceEnt1, faceEnt2);
        adjFace2SubMeshId.push_back(subMeshId1);

        // judge the siding status of the two faces
        VertHandle vh1 = half.topo.from();
        VertHandle vh2 = half.topo.to();
        Line2 line(mesh2d.data(meshProxy.data(vh1)),
                   mesh2d.data(meshProxy.data(vh2)));
        // get a faceCorner which can form a triangle with line that is inside
        // the face polygon
        Point2 faceCorners[2] = {Origin<2>(), Origin<2>()};
        int id = 0;
        for (FaceHandle fh : {fh1, fh2}) {
          TriangulatePolygon(
              meshProxy.topo(fh).halfedges.begin(),
              meshProxy.topo(fh).halfedges.end(),
              [&meshProxy, &mesh2d](HalfHandle hhProxy) -> const Point2 & {
                return mesh2d.data(
                    meshProxy.data(meshProxy.topo(hhProxy).to()));
              },
              [&meshProxy, &mesh2d, &faceCorners, id, vh1,
               vh2](HalfHandle hhProxy1, HalfHandle hhProxy2,
                    HalfHandle hhProxy3) {
                VertHandle vhs[] = {meshProxy.topo(hhProxy1).to(),
                                    meshProxy.topo(hhProxy2).to(),
                                    meshProxy.topo(hhProxy3).to()};
                if (MakeOrderedPair(vh1, vh2) ==
                    MakeOrderedPair(vhs[0], vhs[1])) {
                  faceCorners[id] = mesh2d.data(meshProxy.data(vhs[2]));
                } else if (MakeOrderedPair(vh1, vh2) ==
                           MakeOrderedPair(vhs[1], vhs[2])) {
                  faceCorners[id] = mesh2d.data(meshProxy.data(vhs[0]));
                } else if (MakeOrderedPair(vh1, vh2) ==
                           MakeOrderedPair(vhs[2], vhs[0])) {
                  faceCorners[id] = mesh2d.data(meshProxy.data(vhs[1]));
                }
              });
          assert(faceCorners[id] != Origin<2>());
          id++;
        }

        bool onSameSide = IsOnLeftSide(faceCorners[0], line.first, line.second) ==
                          IsOnLeftSide(faceCorners[1], line.first, line.second);

        adjFaceOverlapInView.push_back(onSameSide);
      }
      int nadjFaces = adjFaceEnts.size();

      // Start Building Matrices
      // vert start position in variable vector
      std::vector<int> ent2varPosition(entities.size(), -1);
      std::vector<int> ent2nvar(entities.size(), -1);
      std::vector<DenseMatd> ent2matFromVarToPlaneCoeffs(entities.size());
      std::vector<int> var2ent;
      int nvars = 0;
      for (int ent = 0; ent < entities.size(); ent ++) {
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

      // A1 : [anchornum x (3 * nents)]
      // A2 : [anchornum x (3 * nents)]
      // A1 * P * X -> inversed anchor depths (left side)
      // A2 * P * X -> inversed anchor depths (right side)
      std::vector<SparseMatElementd> A1triplets, A2triplets;
      int anchorId = 0;
      for (auto &cons : ents2anchors) {
        int ent1 = cons.first.first;
        int ent2 = cons.first.second;

        auto &anchors = cons.second;
        for (auto &anchor : anchors) {
          Vec3 nanchor = normalize(anchor);
          for (int k = 0; k < 3; k++) {
            A1triplets.emplace_back(anchorId, k + ent1 * 3, nanchor[k]);
            A2triplets.emplace_back(anchorId, k + ent2 * 3, nanchor[k]);
          }
          anchorId++;
        }
      }
      int nanchors = anchorId;

      // J1: [(3 * nadjFaces) x (3 * nents)]
      // J1: [(3 * nadjFaces) x (3 * nents)]
      // J1 * P * X -> plane coefficients at edges (left side)
      // J2 * P * X -> plane coefficients at edges (right side)
      std::vector<SparseMatElementd> J1triplets, J2triplets;
      for (int i = 0; i < nadjFaces; i++) {
        int ent1 = adjFaceEnts[i].first;
        int ent2 = adjFaceEnts[i].second;
        for (int k = 0; k < 3; k++) {
          J1triplets.emplace_back(k + i * 3, k + ent1 * 3, 1);
          J2triplets.emplace_back(k + i * 3, k + ent2 * 3, 1);
        }
      }
      std::vector<double> J12Signs;
      // todo

      // solve!
      matlab << "clear()";
      matlab.setVar("P", MakeSparseMatFromElements(3 * nents, nvars,
                                                   Ptriplets.begin(),
                                                   Ptriplets.end()));
      matlab.setVar("A1", MakeSparseMatFromElements(nanchors, 3 * nents,
                                                    A1triplets.begin(),
                                                    A1triplets.end()));
      matlab.setVar("A2", MakeSparseMatFromElements(nanchors, 3 * nents,
                                                    A2triplets.begin(),
                                                    A2triplets.end()));
      matlab.setVar("J1", MakeSparseMatFromElements(3 * nadjFaces, 3 * nents,
                                                    J1triplets.begin(),
                                                    J1triplets.end()));
      matlab.setVar("J2", MakeSparseMatFromElements(3 * nadjFaces, 3 * nents,
                                                    J2triplets.begin(),
                                                    J2triplets.end()));

      matlab.setVar("nvars", nvars);
      matlab.setVar("nents", nents);
      matlab.setVar("nanchors", nanchors);
      matlab.setVar("nadjFaces", nadjFaces);

      double minE = std::numeric_limits<double>::infinity();

      std::vector<double> X;

      static const int maxIter = 10;
      matlab << "PHI = ones(nanchors, 1);";
      matlab << "DELTA = ones(3, nadjFaces);";
      matlab << "AvgCos = 0;";

      for (int t = 0; t < maxIter; t++) {
        const std::string objective =
            "1e5 * sum_square(((A1 - A2) * P * X) ./ PHI)"
            " + "
            "sum_square(dot(reshape(J1 * P * X, [3 nadjFaces]), "
            "DELTA) - AvgCos)";

        matlab << "cvx_begin quiet";
        matlab << "variable X(nvars);";
        matlab << "minimize " + objective + ";";
        matlab << "subject to";
        matlab << " ones(nanchors, 1) <= A1 * P * X;";
        matlab << " ones(nanchors, 1) <= A2 * P * X;";
        matlab << "cvx_end";

        matlab << "PHI = (A1 * P * X) .* (A2 * P * X);";
        matlab << "PHI = abs(PHI) / norm(PHI);";

        matlab << "Planes1 = reshape(J1 * P * X, [3 nadjFaces]);";
        matlab << "Planes2 = reshape(J2 * P * X, [3 nadjFaces]);";
        matlab << "DELTA = Planes2 ./ "
                  "sqrt(repmat(dot(Planes1, Planes1) .* "
                  "dot(Planes2, Planes2), "
                  "[3 1]));"; 
        matlab << "AvgCos = mean(dot(Planes1, Planes2));";

        matlab << "e = " + objective + ";";

        double curEnergy = matlab.var("e");
        Println("cur energy - ", curEnergy);
        if (IsInfOrNaN(curEnergy)) {
          break;
        }
      }
    }
  }
}
