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
  sb.installingOptions().lineWidth = 10.0;
  sb.installingOptions().defaultShaderSource =
      gui::OpenGLShaderSourceDescriptor::XLines;
  for (auto &h : m.halfedges()) {
    Line3 line(m.data(h.topo.from()), m.data(h.topo.to()));
    auto hh = h.topo.hd;
    auto fh = h.topo.face;
    auto oppohh = h.topo.opposite;
    auto oppofh = oppohh.valid() ? m.topo(oppohh).face : FaceHandle();
    bool hasFace = fh.valid();
    bool hasOppo = oppohh.valid();
    gui::Color color = colorHalf(hh);
    if (!hasFace && hasOppo) {
      color = gui::Red;
    } else if (hasFace && !hasOppo) {
      color = gui::Blue;
    } else if (!hasFace && !hasOppo) {
      color = gui::Yellow;
    }
    sb.add(gui::ColorAs(line, color), [hh, fh, oppohh, oppofh](auto &...) {
      std::cout << "halfedge id: " << hh.id
                << ", opposite halfedge id: " << oppohh.id
                << ", face id: " << fh.id << ", opposite face id: " << oppofh.id
                << '\n';
    });
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
    sb.add(gui::ColorAs(poly, colorFace(f.topo.hd)),
           [fh](auto &...) { std::cout << "face id: " << fh.id << '\n'; });
  }
}

int main(int argc, char **argv) {
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");
  misc::Matlab matlab;

  std::string name = "tritower";
  std::string folder = "F:\\DataSets\\linedrawing_"
                       "dataset\\DatabaseFile\\" +
                       name + "\\";
  auto drawing =
      LoadLineDrawing(folder + name + "_cy.3dr", folder + name + ".gt");
  auto sphere = BoundingBox(drawing).outerSphere();

  PerspectiveCamera cam(500, 500, Point2(250, 250), 200,
                        sphere.center + Vec3(1, 1, 1) * sphere.radius * 2,
                        sphere.center);
  {
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XTriangles;
    sb.installingOptions().discretizeOptions.color(gui::Black);
    sb.installingOptions().lineWidth = 3.0;
    // todo, debug the mesh faces
    auto mesh = ToMesh(drawing);
    const gui::ColorTable ctable =
        gui::CreateRandomColorTableWithSize(mesh.internalFaces().size());
    AddToScene(sb, mesh, ConstantFunctor<gui::Color>(gui::Black),
               [&ctable](FaceHandle fh) { return ctable[fh.id]; });
    cam = sb.show(true, true, gui::RenderOptions()
                                  .backgroundColor(gui::White)
                                  .renderMode(gui::Lines)
                                  .bwTexColor(0.0)
                                  .bwColor(1.0)
                                  .fixUpDirectionInCameraMove(false)
                                  .cullBackFace(false)
                                  .cullFrontFace(false))
              .camera();
  }

  // generate perspective 2d line drawing
  auto drawing2d = Transform(drawing, [&cam](const Point3 &p) -> Point3 {
    return cat(cam.toScreen(p), 0.0);
  });

  auto mesh2d = ToMesh(drawing2d);
  int ncc = ConnectedComponents(mesh2d, [](auto &&...) {});
  DecomposeAll(mesh2d, [&mesh2d](HalfHandle h1, HalfHandle h2) -> bool {
    if (h1 == h2) {
      return true;
    }
    return false;
   /* Line3 line1(mesh2d.data(mesh2d.topo(h1).from()),
                mesh2d.data(mesh2d.topo(h1).to()));
    Line3 line2(mesh2d.data(mesh2d.topo(h2).from()),
                mesh2d.data(mesh2d.topo(h2).to()));
    return DistanceBetweenTwoLines(line1, line2).first == 0;*/
  });
  std::map<VertHandle, int> vh2ccid;
  int ncc2 =
      ConnectedComponents(mesh2d, [&vh2ccid](const auto &mesh, VertHandle vh,
                                             int ccid) { vh2ccid[vh] = ccid; });
  std::cout << ncc2 << " connected components after decomposition" << '\n';

  {
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XLines;
    sb.installingOptions().discretizeOptions.color(gui::Black);
    sb.installingOptions().lineWidth = 1.0;
    gui::ColorTable ctable = gui::CreateRandomColorTableWithSize(ncc2);
    ctable.exceptionalColor() = gui::Black;
    AddToScene(sb, mesh2d,
               [&mesh2d, &vh2ccid, &ctable](HalfHandle hh) {
                 return ctable[vh2ccid.at(mesh2d.topo(hh).from())];
               },
               [](FaceHandle fh) { return gui::Transparent; });
    cam = sb.show(true, true, gui::RenderOptions()
                                  .backgroundColor(gui::White)
                                  .renderMode(gui::Lines)
                                  .bwTexColor(0.0)
                                  .bwColor(1.0)
                                  .fixUpDirectionInCameraMove(false)
                                  .cullBackFace(false)
                                  .cullFrontFace(false))
              .camera();
  }

  //
  std::vector<Classified<Line2>> lines(drawing2d.line2corners.size());
  for (int i = 0; i < drawing2d.line2corners.size(); i++) {
    auto &p1 = drawing2d.corners[drawing2d.line2corners[i].first];
    auto &p2 = drawing2d.corners[drawing2d.line2corners[i].second];
    lines[i].component.first = Point2(p1[0], p1[1]);
    lines[i].component.second = Point2(p2[0], p2[1]);
    lines[i].claz = -1;
  }

  VanishingPointsDetector vpd;
  vpd.params().algorithm = VanishingPointsDetector::TardifSimplified;
  std::vector<HPoint2> vps;
  double focal = 0;
  auto result = vpd(lines, cam.screenSize());
  if (result.failed()) {
    return 0;
  }
  std::tie(vps, focal) = result.unwrap();

  {
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XLines;
    sb.installingOptions()
        .discretizeOptions.color(gui::Black)
        .colorTable(gui::RGB);
    sb.installingOptions().discretizeOptions.colorTable().exceptionalColor() =
        gui::Black;
    sb.installingOptions().lineWidth = 3.0;
    std::vector<Classified<Line3>> line3s(lines.size());
    for (int i = 0; i < lines.size(); i++) {
      line3s[i].component.first = cat(lines[i].component.first, 0.0);
      line3s[i].component.second = cat(lines[i].component.second, 0.0);
      line3s[i].claz = lines[i].claz;
    }
    sb.add(line3s);
    sb.show(true, true, gui::RenderOptions()
                            .backgroundColor(gui::White)
                            .renderMode(gui::Lines)
                            .fixUpDirectionInCameraMove(false));
  }

  return 0;
}
