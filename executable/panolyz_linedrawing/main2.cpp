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

template <class HalfColorerFunT, class FaceColorerFunT>
void AddToScene(gui::SceneBuilder &sb, const Mesh<Point3> &m,
                HalfColorerFunT colorHalf, FaceColorerFunT colorFace) {
  sb.installingOptions().lineWidth = 1.0;
  sb.installingOptions().defaultShaderSource =
      gui::OpenGLShaderSourceDescriptor::XLines;
  HandledTable<HalfHandle, int> added(m.internalHalfEdges().size(), false);
  for (auto &h : m.halfedges()) {
    Line3 line(m.data(h.topo.from()), m.data(h.topo.to()));
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
    Polygon3 poly;
    for (auto h : f.topo.halfedges) {
      auto v = m.topo(h).to();
      poly.corners.push_back(m.data(v));
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

  std::string objFile = "H:\\GitHub\\Panoramix\\data\\linedrawing\\tower.obj";
  std::string camFile =
      "H:\\GitHub\\Panoramix\\data\\linedrawing\\tower.obj.cam1.cereal";

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

  PerspectiveCamera cam;
  if (!LoadFromDisk(camFile, cam)) {
    auto sphere = BoundingBoxOfContainer(mesh.vertices()).outerSphere();
    PerspectiveCamera projCam(500, 500, Point2(250, 250), 200,
                              sphere.center + Vec3(1, 2, 3) * sphere.radius * 2,
                              sphere.center);
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XTriangles;
    sb.installingOptions().discretizeOptions.color(gui::Black);
    sb.installingOptions().lineWidth = 0.01;
    const gui::ColorTable ctable =
        gui::CreateRandomColorTableWithSize(subMeshes.size());
    HandledTable<HalfHandle, int> hh2subMeshId(mesh.internalHalfEdges().size(),
                                               -1);
    for (auto &h : mesh.halfedges()) {
      int &id = hh2subMeshId[h.topo.hd];
      for (int i = 0; id == -1 && i < subMeshes.size(); i++) {
        for (auto hh : subMeshes[i].hhs) {
          if (meshProxy.data(hh) == h.topo.hd) {
            id = i;
            break;
          }
        }
      }
    }
    HandledTable<FaceHandle, int> fh2subMeshId(mesh.internalFaces().size(), -1);
    for (auto &f : mesh.faces()) {
      int &id = fh2subMeshId[f.topo.hd];
      for (int i = 0; id == -1 && i < subMeshes.size(); i++) {
        for (auto fh : subMeshes[i].fhs) {
          if (meshProxy.data(fh) == f.topo.hd) {
            id = i;
            break;
          }
        }
      }
    }
    AddToScene(sb, mesh,
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
    SaveToDisk(camFile, cam);
  }

  // convert to 2d
  auto mesh2d = Transform(
      mesh, [&cam](const Point3 &p) -> Point2 { return cam.toScreen(p); });


}