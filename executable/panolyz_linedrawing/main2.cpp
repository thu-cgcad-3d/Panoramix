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

int main(int argc, char **argv) {
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");

  std::string name = "towerx";
  std::string camName = "cam1";

  std::string objFile = "H:\\GitHub\\Panoramix\\data\\linedrawing\\" + name +
                        "\\" + name + ".obj";
  std::string camFile = "H:\\GitHub\\Panoramix\\data\\linedrawing\\" + name +
                        "\\" + name + ".obj." + camName + ".cereal";

  auto mesh = LoadFromObjFile(objFile);

  auto meshProxy = MakeMeshProxy(mesh);

  // decompose
  DecomposeAll(meshProxy,
               [](HalfHandle hh1, HalfHandle hh2) -> bool { return false; });
  auto subMeshes = ExtractSubMeshes(meshProxy,
                                    [](auto hhbegin, auto hhend) -> bool {
                                      return std::distance(hhbegin, hhend) <= 1;
                                    },
                                    10);
  Println("found ", subMeshes.size(), " subMeshes");

  PerspectiveCamera cam;
  if (!LoadFromDisk(camFile, cam)) {
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
    if (false) {
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

  {
    Image3ub im(cam.screenSize(), Vec3ub(255, 255, 255));
    auto canvas = gui::MakeCanvas(im);
    canvas.color(gui::Black);
    canvas.thickness(2);
    for (auto &h : mesh2d.halfedges()) {
      auto p1 = mesh2d.data(h.topo.from());
      auto p2 = mesh2d.data(h.topo.to());
      canvas.add(Line2(p1, p2));
    }
    canvas.show(0, "mesh2d");
  }

  // build the graph
}