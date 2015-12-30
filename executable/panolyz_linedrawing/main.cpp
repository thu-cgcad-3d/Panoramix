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

#define CONCAT_IMPL(x, y) x##y
#define MACRO_CONCAT(x, y) CONCAT_IMPL(x, y)
#define DISABLED_main MACRO_CONCAT(main_, __COUNTER__)

int main(int argc, char **argv) {
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Cache\\Panoramix\\LineDrawing\\");
  misc::Matlab matlab;

  std::string name = "tritower";
  std::string folder = "F:\\DataSets\\linedrawing_"
                       "dataset\\DatabaseFile\\" +
                       name + "\\";
  auto drawing =
      LoadLineDrawing(folder + name + "_cy.3dr", folder + name + ".gt");
  auto sphere = BoundingBoxOfContainer(drawing.vertPositions).outerSphere();
  double scale = sphere.radius;

  {
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XLines;
    sb.installingOptions().discretizeOptions.color = gui::Black;
    sb.installingOptions().lineWidth = 3.0;
    for (auto &l : drawing.line2verts) {
      Line3 line3(drawing.vertPositions[l.first],
                  drawing.vertPositions[l.second]);
      sb.add(line3 / scale);
    }
    sb.show(true, true, gui::RenderOptions()
                            .backgroundColor(gui::White)
                            .renderMode(gui::Lines));
  }

  // generate perspective 2d line drawing
  PerspectiveCamera cam(500, 500, Point2(250, 250), 200,
                        sphere.center + Vec3(1, 1, 1) * sphere.radius * 2,
                        sphere.center);
  auto drawing2d = Transform(drawing, [&cam](const Point3 &p) -> Point3 {
    return cat(cam.toScreen(p), 0.0);
  });

  //
  std::vector<Classified<Line2>> lines(drawing2d.line2verts.size());
  for (int i = 0; i < drawing2d.line2verts.size(); i++) {
    auto &p1 = drawing2d.vertPositions[drawing2d.line2verts[i].first];
    auto &p2 = drawing2d.vertPositions[drawing2d.line2verts[i].second];
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
    double scale =
        BoundingBoxOfContainer(drawing2d.vertPositions).outerSphere().radius;
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XLines;
    sb.installingOptions().discretizeOptions.color = gui::Black;
    sb.installingOptions().lineWidth = 3.0;
    sb.installingOptions().discretizeOptions.colorTable = gui::RGB;
    sb.installingOptions().discretizeOptions.colorTable.exceptionalColor() =
        gui::Black;
    std::vector<Classified<Line3>> line3s(lines.size());
    for (int i = 0; i < lines.size(); i++) {
      line3s[i].component.first = cat(lines[i].component.first / scale, 0.0);
      line3s[i].component.second = cat(lines[i].component.second / scale, 0.0);
      line3s[i].claz = lines[i].claz;
    }
    sb.add(line3s);
    sb.show(true, true, gui::RenderOptions()
                            .backgroundColor(gui::White)
                            .renderMode(gui::Lines));
  }

  return 0;
}
