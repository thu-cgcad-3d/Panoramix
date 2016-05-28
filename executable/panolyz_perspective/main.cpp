#include "parallel.hpp"
#include "cache.hpp"
#include "clock.hpp"

#include "canvas.hpp"
#include "singleton.hpp"
#include "gui_util.hpp"

#include "pi_graph_annotation.hpp"
#include "pi_graph_cg.hpp"
#include "pi_graph_control.hpp"
#include "pi_graph_occlusion.hpp"
#include "pi_graph_optimize.hpp"
#include "pi_graph_vis.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

struct PerspectiveOptions {
  bool refresh_preparation;
  bool refresh_gc;
  bool refresh_mg_init;
  bool refresh_mg_oriented;
  bool refresh_mg_reconstructed;

  bool usePrincipleDirectionPrior;
  bool useWallPrior;
  bool useGeometricContextPrior;
};

int main(int argc, char **argv) {

  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\Perspective\\");
  misc::Matlab matlab;

  PerspectiveOptions options;
  options.refresh_preparation = false;
  options.refresh_gc = false;
  options.refresh_mg_init = false;
  options.refresh_mg_oriented = true;
  options.refresh_mg_reconstructed = true;

  options.usePrincipleDirectionPrior = true;
  options.useGeometricContextPrior = true;
  options.useWallPrior = true;

  std::vector<std::string> dirs;
  auto images = gui::PickImages("F:\\DataSets\\CVPR2016EXT", &dirs);

  for (int i = 0; i < dirs.size(); i++) {
    Image3ub original = images[i];
    std::string dir = dirs[i];

    Println("processing ", dir);

    Image3ub im = original.clone();
    ResizeToArea(im, 500 * 500);

    View<PerspectiveCamera, Image3ub> perspView;
    View<PanoramicCamera, Image3ub> view;
    Imageb maskInPano;

    std::vector<Classified<Line3>> line3s;
    std::vector<Classified<Line2>> line2s;
    std::vector<Vec3> vps;
    double focal = 0.0;
    Imagei segs;
    int nsegs;

    // view, line, vps, focal, segs
    if (options.refresh_preparation ||
        !misc::LoadCache(dir, "preparation", perspView, view, maskInPano,
                         line3s, line2s, vps, focal, segs, nsegs)) {
      VanishingPointsDetector vpd;
      vpd.params().algorithm = VanishingPointsDetector::TardifSimplified;
      auto result =
          CreatePerspectiveView(im, Origin(), X(), -Z(), LineSegmentExtractor(),
                                vpd, &line3s, &line2s, &vps, &focal);
      if (result.failed()) {
        Println("failed calibrating the camera");
        continue;
      }
      perspView = result.unwrap();
      OrderVanishingPoints(vps, perspView.camera.upward(), line2s, line3s);

      // get panoramic view
      view = CreatePanoramicView(Image3ub(700, 1400, Vec3ub()));
      auto sampler = MakeCameraSampler(view.camera, perspView.camera);
      view.image = sampler(perspView.image, cv::BORDER_CONSTANT, Vec3ub());
      maskInPano = sampler(Imageb(perspView.image.size(), true),
                           cv::BORDER_CONSTANT, false);

      // segs
      // estimate segs
      nsegs = SegmentationForPIGraph(view, line3s, segs, DegreesToRadians(1));
      RemoveThinRegionInSegmentation(segs, 1, false);
      RemoveEmbededRegionsInSegmentation(segs, false);
      nsegs = DensifySegmentation(segs, false);
      assert(IsDenseSegmentation(segs));

      if (true) {
        auto ctable = gui::CreateGreyColorTableWithSize(nsegs);
        ctable.randomize();
        gui::ColorTable rgb = gui::RGBGreys;
        auto canvas = gui::MakeCanvas(view.image).alpha(0.9).add(ctable(segs));
        for (auto &l : line3s) {
          static const double sampleAngle = M_PI / 100.0;
          auto &line = l.component;
          double spanAngle = AngleBetweenDirected(line.first, line.second);
          std::vector<Point2> ps;
          ps.reserve(spanAngle / sampleAngle);
          for (double angle = 0.0; angle <= spanAngle; angle += sampleAngle) {
            Vec3 dir = RotateDirection(line.first, line.second, angle);
            ps.push_back(view.camera.toScreen(dir));
          }
          for (int i = 1; i < ps.size(); i++) {
            auto &p1 = ps[i - 1];
            auto &p2 = ps[i];
            if (Distance(p1, p2) >= view.image.cols / 2) {
              continue;
            }
            canvas.thickness(2);
            canvas.colorTable(rgb).add(gui::ClassifyAs(Line2(p1, p2), l.claz));
          }
        }
        canvas.show();
      }

      misc::SaveCache(dir, "preparation", perspView, view, maskInPano, line3s,
                      line2s, vps, focal, segs, nsegs);
    }

    // gc
    Image5d gc;
    if (options.refresh_gc || !misc::LoadCache(dir, "gc", gc)) {
      gc = ComputeGeometricContext(matlab, im, false, true);
      misc::SaveCache(dir, "gc", gc);
    }

    // panoramic view

    // build pi graph
    static const double thetaTiny = DegreesToRadians(2);
    static const double thetaMid = DegreesToRadians(5);
    static const double thetaLarge = DegreesToRadians(15);

    int vertVPId = 0;

    PIGraph<PanoramicCamera> mg;
    if (options.refresh_mg_init || !misc::LoadCache(dir, "mg_init", mg)) {
      mg = BuildPIGraph(view, vps, vertVPId, segs, line3s, DegreesToRadians(1),
                        DegreesToRadians(1), DegreesToRadians(1),
                        ///!!!0.04,
                        thetaTiny, thetaLarge, thetaTiny);
      mg.seg2control[segs(Pixel(0, 0))].used = false;
      misc::SaveCache(dir, "mg_init", mg);
    }

    // attach orientation constraints
    if (options.refresh_mg_oriented ||
        !misc::LoadCache(dir, "mg_oriented", mg)) {
      std::cout << "########## refreshing mg oriented ###########" << std::endl;
      if (options.usePrincipleDirectionPrior) {
        AttachPrincipleDirectionConstraints(mg);
      }
      if (options.useWallPrior) {
        AttachWallConstraints(mg, thetaTiny);
      }
      if (options.useGeometricContextPrior) {
        AttachGCConstraints(mg, MakeView(gc, perspView.camera), 0.7, 0.7, true);
      }
      misc::SaveCache(dir, "mg_oriented", mg);
    }

    const auto printPIGraph = [&mg, &dir](int delay,
                                          const std::string &saveAs) {
      static const gui::ColorTable randColors =
          gui::CreateRandomColorTableWithSize(mg.nsegs);
      static gui::ColorTable rgbColors = gui::RGBGreys;
      rgbColors.exceptionalColor() = gui::Gray;
      auto pim = Print2(
          mg,
          [&mg](int seg, Pixel pos) -> gui::Color {
            static const gui::ColorTable ctable =
                gui::ColorTableDescriptor::RGBGreys;
            auto &c = mg.seg2control[seg];
            if (!c.used) {
              return gui::Black;
            }
            /*if (c.orientationClaz == 0) {
                return gui::Red;
            }*/
            if (c.orientationClaz != -1) {
              return ctable[c.orientationClaz].blendWith(gui::White, 0.3);
            }
            if (c.orientationNotClaz != -1) {
              static const int w = 10;
              if (IsBetween((pos.x + pos.y) % w, 0, w / 2)) {
                return ctable[c.orientationNotClaz].blendWith(gui::White, 0.3);
              } else {
                return gui::White;
              }
            }
            return gui::White;
          },
          [&mg](int lp) { return gui::Transparent; },
          [&mg](int bp) -> gui::Color { return gui::Black; }, 1, 0);
      auto canvas = gui::AsCanvas(pim);
      for (auto &l : mg.lines) {
        static const double sampleAngle = M_PI / 100.0;
        auto &line = l.component;
        int claz = l.claz;
        if (claz >= mg.vps.size()) {
          claz = -1;
        }
        double spanAngle = AngleBetweenDirected(line.first, line.second);
        std::vector<Point2> ps;
        ps.reserve(spanAngle / sampleAngle);
        for (double angle = 0.0; angle <= spanAngle; angle += sampleAngle) {
          Vec3 dir = RotateDirection(line.first, line.second, angle);
          ps.push_back(mg.view.camera.toScreen(dir));
        }
        for (int i = 1; i < ps.size(); i++) {
          auto p1 = ToPixel(ps[i - 1]);
          auto p2 = ToPixel(ps[i]);
          if (Distance(p1, p2) >= mg.view.image.cols / 2) {
            continue;
          }
          gui::Color color = rgbColors[claz];
          cv::clipLine(cv::Rect(0, 0, canvas.image().cols, canvas.image().rows),
                       p1, p2);
          cv::line(canvas.image(), p1, p2, (cv::Scalar)color / 255.0, 2);
        }
      }
      canvas.show(delay, "pi graph");
      if (saveAs != "") {
        cv::imwrite(saveAs, Image3ub(canvas.image() * 255));
      }
      return canvas.image();
    };

    printPIGraph(0, "");

    // connect all
    {
      for (int bp = 0; bp < mg.nbndPieces(); bp++) {
        mg.bndPiece2segRelation[bp] = SegRelation::Connected;
      }
      for (int lp = 0; lp < mg.nlinePieces(); lp++) {
        mg.linePiece2segLineRelation[lp] = SegLineRelation::Attached;
      }
      for (int lr = 0; lr < mg.nlineRelations(); lr++) {
        mg.lineRelations[lr] = LineRelation::Attached;
      }
    }

    PIConstraintGraph cg;
    PICGDeterminablePart dp;
    if (options.refresh_mg_reconstructed ||
        !misc::LoadCache(dir, "mg_reconstructed", mg, cg, dp)) {
      std::cout << "########## refreshing mg reconstructed ###########"
                << std::endl;
      cg = BuildPIConstraintGraph(mg, DegreesToRadians(1), 0.01);
      dp = LocateDeterminablePart(cg, DegreesToRadians(3), false);
      auto start = std::chrono::system_clock::now();
      double energy = Solve(dp, cg, matlab, 5, 1e6);
      if (IsInfOrNaN(energy)) {
        std::cout << "solve failed" << std::endl;
        continue;
      }
      misc::SaveCache(dir, "mg_reconstructed", mg, cg, dp);
    }

    VisualizeReconstructionCompact(view.image, dp, cg, mg, true);
  }
}