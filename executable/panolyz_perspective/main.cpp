#include "../../src/core/parallel.hpp"
#include "../../src/misc/cache.hpp"
#include "../../src/misc/clock.hpp"

#include "../../src/gui/canvas.hpp"
#include "../../src/gui/singleton.hpp"
#include "../../src/gui/utility.hpp"

#include "../../src/experimental/pi_graph_annotation.hpp"
#include "../../src/experimental/pi_graph_cg.hpp"
#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
#include "../../src/experimental/pi_graph_optimize.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

struct PerspectiveOptions {
  bool refresh_preparation;
  bool refresh_gc;
  bool refresh_mg_init;
};

int main(int argc, char **argv) {

  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\Perspective\\");
  misc::Matlab matlab;

  PerspectiveOptions options;
  options.refresh_preparation = false;
  options.refresh_gc = false;
  options.refresh_mg_init = false;

  std::vector<std::string> dirs;
  auto images = gui::PickImages("H:\\GitHub\\Panoramix\\data\\normal", &dirs);

  for (int i = 0; i < dirs.size(); i++) {
    Image3ub original = images[i];
    std::string dir = dirs[i];

    Println("processing ", dir);

    Image3ub im = original.clone();
    ResizeToArea(im, 500 * 500);

    PerspectiveView view;
    std::vector<Classified<Line3>> line3s;
    std::vector<Classified<Line2>> line2s;
    std::vector<Vec3> vps;
    double focal = 0.0;
    Imagei segs;
    int nsegs;

    if (options.refresh_preparation ||
        !misc::LoadCache(dir, "preparation", view, line3s, line2s, vps, focal,
                         segs, nsegs)) {
      VanishingPointsDetector vpd;
      vpd.params().algorithm = VanishingPointsDetector::TardifSimplified;
      auto result =
          CreatePerspectiveView(im, Origin(), X(), -Z(), LineSegmentExtractor(),
                                vpd, &line3s, &line2s, &vps, &focal);
      if (result.failed()) {
        Println("failed calibrating the camera");
        continue;
      }
      view = result.unwrap();

      if (false) {
        gui::AsCanvas(im)
            .colorTable(gui::ColorTableDescriptor::RGBGreys)
            .thickness(2.0)
            .add(line2s)
            .show();
      }

      SegmentationExtractor segmenter;
      segmenter.params().algorithm = SegmentationExtractor::GraphCut;
      std::tie(segs, nsegs) = segmenter(im, false);

      misc::SaveCache(dir, "preparation", view, line3s, line2s, vps, focal,
                      segs, nsegs);
    }

    Image5d gc;
    if (options.refresh_gc || !misc::LoadCache(dir, "gc", gc)) {
      gc = ComputeGeometricContext(matlab, im, false, true);
      misc::SaveCache(dir, "gc", gc);
    }

    // build pi graph
    PIGraph<PerspectiveCamera> mg;
    if (options.refresh_mg_init || !misc::LoadCache(dir, "mg_init", mg)) {
        //mg = BuildPIGraph()
    }

  }
}