#include <QtWidgets>

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
    std::vector<Point2> vps;
    vps.reserve(chain.size() * chain.size() / 4);
    for (int i = 0; i < chain.size() / 2; i++) {
      for (int j = i + 1; j < i + chain.size() / 2; j++) {
        vps.push_back(Intersection(Line2(chain.at(i), chain.at(j)).ray(),
                                   Line2(chain.at(i + chain.size() / 2),
                                         chain.at(j + chain.size() / 2))
                                       .ray()));
      }
    }
    return vps;
  } else {
    std::vector<Point2> vps;
    for (int i = 0; i < chain.size(); i++) {
      for (int j = i + 1; j < chain.size(); j++) {
        vps.push_back(Intersection(chain.edge(i).ray(), chain.edge(j).ray()));
      }
    }
    return vps;
  }
}

int main(int argc, char **argv) {
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");

  std::string name = "hex";
  std::string camName = "cam2";
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

  std::vector<std::pair<std::set<int>, PPFocalCandidate>> ppFocalIdGroups;
  { // dbscan
    std::vector<int> ppFocalId2group(ppFocalCandidates.size(), -1);
    int ngroups = 0;
    RTreeMap<Vec3, int> ppFocalIdTree;
    for (int i = 0; i < ppFocalCandidates.size(); i++) {
      Vec3 coordinate =
          cat(ppFocalCandidates[i].pp, ppFocalCandidates[i].focal);
      const double thres = scale / 50.0;
      ppFocalIdTree.search(BoundingBox(coordinate).expand(thres * 2),
                           [&coordinate, thres, i, &ppFocalId2group](
                               const std::pair<Vec3, int> &cand) {
                             if (Distance(cand.first, coordinate) <= thres) {
                               ppFocalId2group[i] =
                                   ppFocalId2group[cand.second];
                               return false;
                             }
                             return true;
                           });
      if (ppFocalId2group[i] == -1) {
        ppFocalId2group[i] = ngroups++;
      }
      ppFocalIdTree.emplace(coordinate, i);
    }

    ppFocalIdGroups.resize(ngroups);
    for (auto &g : ppFocalIdGroups) {
      g.second.focal = 0.0;
      g.second.pp = Point2();
    }
    for (int i = 0; i < ppFocalId2group.size(); i++) {
      auto &g = ppFocalIdGroups[ppFocalId2group[i]];
      g.first.insert(i);
      g.second.focal += ppFocalCandidates[i].focal;
      g.second.pp += ppFocalCandidates[i].pp;
    }
    for (auto &g : ppFocalIdGroups) {
      g.second.focal /= g.first.size();
      g.second.pp /= double(g.first.size());
    }

    std::sort(
        ppFocalIdGroups.begin(), ppFocalIdGroups.end(),
        [](auto &g1, auto &g2) { return g1.first.size() > g2.first.size(); });
  }



}