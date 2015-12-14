#include "../src/core/cameras.hpp"
#include "../src/core/feature.hpp"
#include "../src/core/utility.hpp"
#include "../src/gui/canvas.hpp"
#include "../src/gui/utility.hpp"

#include "../src/experimental/pi_graph.hpp"

#include "config.hpp"

using namespace pano;
using namespace test;

TEST(Feature, SegmentationExtractor) {
  core::Image3ub im = gui::PickAnImage();
  if (im.empty())
    return;

  core::ResizeToMakeHeightUnder(im, 600);
  {
    core::SegmentationExtractor::Params p;
    p.c = 5;
    p.minSize = 400;
    p.sigma = 1;
    core::SegmentationExtractor seg(p);
    gui::AsCanvas(im).show();
    auto segs = seg(im);
    gui::AsCanvas(gui::CreateRandomColorTableWithSize(segs.second)(segs.first))
        .show();
  }
  {
    core::SegmentationExtractor::Params p;
    p.c = 5;
    p.minSize = 400;
    p.sigma = 1;
    core::SegmentationExtractor seg(p);
    auto segs = seg(
        im, {core::Line2({0.0, 0.0}, core::Point2(im.cols, im.rows)),
             core::Line2(core::Point2(im.cols, 0), core::Point2(0, im.rows))});
    gui::AsCanvas(gui::CreateRandomColorTableWithSize(segs.second)(segs.first))
        .show();
  }
  {
    core::SegmentationExtractor::Params p;
    p.algorithm = core::SegmentationExtractor::SLIC;
    p.superpixelSizeSuggestion = 3000;
    core::SegmentationExtractor seg(p);
    gui::AsCanvas(im).show();
    auto segs = seg(im);
    gui::AsCanvas(gui::CreateRandomColorTableWithSize(segs.second)(segs.first))
        .show();
  }
  {
    core::SegmentationExtractor::Params p;
    p.algorithm = core::SegmentationExtractor::QuickShiftCPU;
    core::SegmentationExtractor seg(p);
    gui::AsCanvas(im).show();
    auto segs = seg(im);
    gui::AsCanvas(gui::CreateRandomColorTableWithSize(segs.second)(segs.first))
        .show();
  }
}

TEST(Feature, SegmentationBoundaryJunction) {
  core::Image im =
      core::ImageRead(ProjectDataDirStrings::PanoramaOutdoor + "/univ0.jpg");
  core::ResizeToMakeHeightUnder(im, 800);
  core::SegmentationExtractor::Params p;
  auto segs = core::SegmentationExtractor(p)(im, true);
  auto junctions = core::ExtractBoundaryJunctions(segs.first);
  std::vector<core::Pixel> ps;
  for (auto &j : junctions) {
    ps.insert(ps.end(), j.second);
  }
  auto ctable = gui::CreateRandomColorTableWithSize(segs.second);
  auto image = ctable(segs.first);
  gui::AsCanvas(image).color(gui::Black).add(ps).show();
}

TEST(Feature, SegmentationExtractorInPanorama) {
  // core::Image im = core::ImageRead(ProjectDataDirStrings::PanoramaOutdoor +
  // "/univ0.jpg");
  core::Image3ub im = gui::PickAnImage(ProjectDataDirStrings::PanoramaIndoor);
  auto cam = core::CreatePanoramicCamera(im);
  if (im.empty()) {
    return;
  }
  core::ResizeToMakeHeightUnder(im, 700);
  core::Imagei segs;
  experimental::SegmentationForPIGraph(
      core::CreatePanoramicView(im),
      std::vector<core::Classified<core::Line3>>(), segs);
  int nsegs = core::MinMaxValOfImage(segs).second + 1;

  core::PerspectiveCamera cam1(500, 500, core::Point2(250, 250), 200,
                               core::Origin(), core::Z(), core::Y()),
      cam2(500, 500, core::Point2(250, 250), 200, core::Origin(), -core::Z(),
           core::Y());

  static const std::string folder =
      "F:\\GitHub\\write-papers\\papers\\a\\figure\\supp\\";

  double alpha = 0.7;
  double beta = 1 - alpha;

  core::Image3ub imsegs =
      gui::CreateRandomColorTableWithSize(nsegs, gui::Transparent)(segs);
  imsegs = imsegs * alpha + im * beta;

  gui::AsCanvas(im).saveAs(folder + "pano.png");
  gui::AsCanvas(core::MakeView(im, cam).sampled(cam1).image)
      .show()
      .saveAs(folder + "pano_polar1.png");
  gui::AsCanvas(core::MakeView(im, cam).sampled(cam2).image)
      .show()
      .saveAs(folder + "pano_polar2.png");

  // gui::AsCanvas(imsegs).show().saveAs(folder + "panosegs.png");
  // gui::AsCanvas(core::MakeView(imsegs,
  // cam).sampled(cam1).image).show().saveAs(folder + "panosegs_polar1.png");
  // gui::AsCanvas(core::MakeView(imsegs,
  // cam).sampled(cam2).image).show().saveAs(folder + "panosegs_polar2.png");

  core::SegmentationExtractor segmenter;
  segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
  segmenter.params().sigma = 10.0;
  segmenter.params().c = 1.0;
  segmenter.params().minSize = 200;
  auto segs2 = segmenter(im).first;
  int nsegs2 = core::MinMaxValOfImage(segs2).second + 1;
  core::Image3ub imsegs2 =
      gui::CreateRandomColorTableWithSize(nsegs2, gui::Transparent)(segs2);
  imsegs2 = imsegs2 * alpha + im * beta;

  // gui::AsCanvas(imsegs2).show().saveAs(folder + "originalsegs.png");
  // gui::AsCanvas(core::MakeView(imsegs2,
  // cam).sampled(cam1).image).show().saveAs(folder +
  // "originalsegs_polar1.png");
  // gui::AsCanvas(core::MakeView(imsegs2,
  // cam).sampled(cam2).image).show().saveAs(folder +
  // "originalsegs_polar2.png");
}

TEST(Feature, RemoveSmallRegionInSegmentation) {
  core::Image3ub im = gui::PickAnImage(ProjectDataDirStrings::PanoramaIndoor);
  core::ResizeToMakeHeightUnder(im, 800);
  core::SegmentationExtractor segmenter;
  segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
  segmenter.params().sigma = 10.0;
  segmenter.params().c = 1.0;
  segmenter.params().superpixelSizeSuggestion = 2000;
  core::Imagei segs;
  int nsegs = 0;
  std::tie(segs, nsegs) = segmenter(im, true);
  gui::AsCanvas(gui::CreateRandomColorTableWithSize(nsegs)(segs)).show();

  EXPECT_TRUE(core::IsDenseSegmentation(segs));
  nsegs = core::RemoveSmallRegionInSegmentation(segs, 1000, true);
  EXPECT_TRUE(core::IsDenseSegmentation(segs));

  gui::AsCanvas(gui::CreateRandomColorTableWithSize(nsegs)(segs)).show();
}

TEST(Feature, RemoveThinRegionInSegmentation) {
  core::Image3ub im = gui::PickAnImage(ProjectDataDirStrings::PanoramaIndoor);
  core::ResizeToMakeHeightUnder(im, 800);
  core::SegmentationExtractor segmenter;
  segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
  segmenter.params().sigma = 10.0;
  segmenter.params().c = 1.0;
  segmenter.params().superpixelSizeSuggestion = 2000;
  core::Imagei segs;
  int nsegs = 0;
  std::tie(segs, nsegs) = segmenter(im, true);
  core::Imagei segs2 = segs.clone();
  core::RemoveThinRegionInSegmentation(segs2, 1, true);
  int nsegs2 = core::DensifySegmentation(segs2);
  EXPECT_TRUE(core::IsDenseSegmentation(segs2));
  gui::AsCanvas(gui::CreateRandomColorTableWithSize(nsegs)(segs)).show();
  gui::AsCanvas(gui::CreateRandomColorTableWithSize(nsegs2)(segs2)).show();
}

TEST(Feature, LineSegmentExtractor) {
  core::LineSegmentExtractor lineseg;
  core::Image3ub im = core::ImageRead(ProjectDataDirStrings::LocalManhattan +
                                      "/buildings2.jpg");
  core::LineSegmentExtractor::Params params;
  params.algorithm = core::LineSegmentExtractor::LSD;
  core::LineSegmentExtractor lineseg2(params);
  gui::AsCanvas(im)
      .color(gui::ColorTag::Yellow)
      .thickness(2)
      .add(lineseg2(im))
      .show();
}

TEST(Feature, VanishingPointsDetector) {

  std::vector<std::string> filenames = {
      //"buildings.jpg",
      //"room.png",
      //"room2e.jpg",
      //"room3.jpg",
      /* "room4.jpg",
       "room5.jpg",
       "room6.jpg",*/
      "room7.jpg",  "room8.jpg",  "room10.jpg", "room11.jpg", "room12.png",
      "room13.jpg", "room14.jpg", "room15.jpg", "room16.jpg", "room17.jpg",
      "room18.jpg", "room19.jpg", "room20.jpg", "room21.jpg", "room22.jpg"};

  auto images = gui::PickImages("F:\\DataSets\\YorkUrbanDB\\data\\");

  core::LineSegmentExtractor::Params lsParams;
  lsParams.minLength = 20;
  lsParams.xBorderWidth = lsParams.yBorderWidth = 20;
  core::LineSegmentExtractor lineseg(lsParams);
  core::VanishingPointsDetector vpdetector;

  std::vector<std::string> failedFileNames;

  for (auto &image : images) {
    // std::cout << "testing image file: " << filename << std::endl;
    // core::Image3ub im = core::ImageRead(ProjectDataDirStrings::Normal + "/" +
    // filename);
    core::Image3ub im = image;
    core::ResizeToMakeWidthUnder(im, 400);

    std::vector<core::HPoint2> vps;
    double focalLength;

    std::vector<core::Classified<core::Line2>> classifiedLines =
        core::ClassifyEachAs(lineseg(im, 3), -1);
    vpdetector.params().algorithm =
        core::VanishingPointsDetector::TardifSimplified;
    auto result = vpdetector(classifiedLines, im.size());
    if (result.null()) {
      std::cout << "failed to find vanishing points!" << std::endl;
      // failedFileNames.push_back(filename);
      continue;
    }
    std::tie(vps, focalLength) = result.unwrap();

    std::vector<core::Classified<core::Ray2>> vpRays;
    for (int i = 0; i < 3; i++) {
      std::cout << "vp[" << i << "] = " << vps[i].value() << std::endl;
      for (double a = 0; a <= M_PI * 2.0; a += 0.1) {
        core::Point2 p = core::Point2(im.cols / 2, im.rows / 2) +
                         core::Vec2(cos(a), sin(a)) * 1000.0;
        vpRays.push_back(core::ClassifyAs(
            core::Ray2(p, (vps[i] - core::HPoint2(p, 1.0)).numerator), i));
      }
    }
    gui::AsCanvas(im)
        .colorTable(gui::ColorTable(gui::ColorTableDescriptor::RGB)
                        .appendRandomizedGreyColors(vps.size() - 3))
        .thickness(1)
        .add(vpRays)
        .thickness(2)
        .add(classifiedLines)
        .show(0);
  }

  for (auto &filename : failedFileNames) {
    std::cout << "failed file: " << filename << std::endl;
  }
}

TEST(Feature, LocalManhattanVanishingPointDetector) {

  using namespace core;

  // forged experiment for panorama
  core::Image3ub im =
      core::ImageRead(ProjectDataDirStrings::PanoramaIndoor + "/14.jpg");
  core::ResizeToMakeWidthUnder(im, 2000);
  // gui::Visualizer2D(im) << gui::manip2d::Show();

  core::PanoramicCamera ocam(im.cols / M_PI / 2.0);
  core::PerspectiveCamera cam(800, 800, core::Point2(400, 400), ocam.focal(),
                              {0, 0, 0}, {-2, 0, -0.5}, {0, 0, -1});

  auto pim = core::MakeCameraSampler(cam, ocam)(im);
  // gui::Visualizer2D(pim) << gui::manip2d::Show();

  core::LineSegmentExtractor lineseg;
  lineseg.params().minLength = 5;
  lineseg.params().algorithm = core::LineSegmentExtractor::LSD;
  std::vector<Line2> line2s = lineseg(pim);

  Vec3 vp1 = {0, 0, 1};
  std::vector<Classified<Line3>> line3s(line2s.size());
  std::vector<Vec3> line3norms(line2s.size());
  for (int i = 0; i < line2s.size(); i++) {
    line3s[i].component.first = normalize(cam.toSpace(line2s[i].first));
    line3s[i].component.second = normalize(cam.toSpace(line2s[i].second));
    line3norms[i] = line3s[i].component.first.cross(line3s[i].component.second);
    line3s[i].claz = abs(line3norms[i].dot(vp1)) < 0.006 ? 0 : -1;
  }

  std::vector<std::pair<int, int>> pairs;
  for (int i = 0; i < line2s.size(); i++) {
    if (line3s[i].claz == 0)
      continue;
    if (abs(line3norms[i].dot(vp1)) < 0.01)
      continue;
    for (int j = i + 1; j < line2s.size(); j++) {
      if (line3s[i].claz == 0)
        continue;
      if (abs(line3norms[j].dot(vp1)) < 0.01)
        continue;
      double dist = DistanceBetweenTwoLines(line2s[i], line2s[j]).first;
      auto &n1 = line3norms[i];
      auto &n2 = line3norms[j];
      auto inter = n1.cross(n2);
      auto interp = cam.toScreen(inter);
      double dd = 40;
      if (dist < dd && DistanceFromPointToLine(interp, line2s[i]).first < dd &&
          DistanceFromPointToLine(interp, line2s[j]).first < dd) {
        pairs.emplace_back(i, j);
      }
    }
  }

  std::vector<std::pair<int, int>> orthoPairs;
  for (auto &p : pairs) {
    auto &n1 = line3norms[p.first];
    auto &n2 = line3norms[p.second];
    auto p1 = normalize(n1.cross(vp1));
    auto p2 = normalize(n2.cross(vp1));
    if (abs(p1.dot(p2)) < 0.02)
      orthoPairs.push_back(p);
  }

  auto viz = gui::AsCanvas(pim).thickness(2);
  for (int i = 0; i < line2s.size(); i++) {
    if (line3s[i].claz == 0) {
      viz.color(gui::ColorTag::Red).add(line2s[i]);
    } else {
      // viz << gui::manip2d::SetColor(gui::ColorTag::Black) << line2s[i];
    }
  }
  for (auto &op : orthoPairs) {
    auto &n1 = line3norms[op.first];
    auto &n2 = line3norms[op.second];
    auto inter = n1.cross(n2);
    auto interp = cam.toScreen(inter);
    viz.color(gui::ColorTag::LightGray)
        .thickness(1)
        .add(Line2(line2s[op.first].center(), interp))
        .add(Line2(line2s[op.second].center(), interp));
    viz.color(gui::ColorTag::White)
        .thickness(2)
        .add(line2s[op.first])
        .add(line2s[op.second]);
  }

  viz.show();
}

TEST(Feature, FeatureExtractor) {
  core::SegmentationExtractor segmenter;
  core::LineSegmentExtractor::Params params;
  params.algorithm = core::LineSegmentExtractor::LSD;
  core::LineSegmentExtractor lineSegmentExtractor(params);

  for (int i = 0; i < 4; i++) {
    std::string name = ProjectDataDirStrings::Normal + "/" + "sampled_" +
                       std::to_string(i) + ".png";
    core::Image3ub im = core::ImageRead(name);
    auto segs = segmenter(im);
    gui::AsCanvas(im)
        .colorTable(gui::CreateRandomColorTableWithSize(segs.second))
        .add(segs.first)
        .add(lineSegmentExtractor(im))
        .show();
  }
}

TEST(Feature, OcclusionDetection) {

  auto images = gui::PickImages("H:\\GitHub\\Panoramix\\data\\normal");
  misc::Matlab matlab;
  for (core::Image3ub &&im : images) {
    core::ResizeToMakeHeightUnder(im, 400);
    auto bd = core::DetectOcclusionBoundary(matlab, im);
    gui::MakeCanvas(im).color(gui::Red).thickness(2.0).add(bd).show();
  }
}

TEST(Feature, ExtractSegmentationTopology) {
  using namespace core;

  Image3ub im = gui::PickAnImage("H:\\GitHub\\Panoramix\\data\\normal");
  if (im.empty())
    return;

  core::ResizeToMakeHeightUnder(im, 600);

  core::SegmentationExtractor::Params p;
  p.c = 5;
  p.minSize = 400;
  p.sigma = 1;
  core::SegmentationExtractor seg(p);
  gui::AsCanvas(im).show();
  auto segs = seg(im, true);
  gui::AsCanvas(gui::CreateRandomColorTableWithSize(segs.second)(segs.first))
      .show();

  std::vector<std::vector<Pixel>> bndpixels;
  std::vector<Pixel> juncpositions;
  std::vector<std::vector<int>> seg2bnds;
  std::vector<std::pair<int, int>> bnd2segs;
  std::vector<std::vector<int>> seg2juncs;
  std::vector<std::vector<int>> junc2segs;
  std::vector<std::pair<int, int>> bnd2juncs;
  std::vector<std::vector<int>> junc2bnds;
  core::ExtractSegmentationTopology(segs.first, bndpixels, juncpositions,
                                    seg2bnds, bnd2segs, seg2juncs, junc2segs,
                                    bnd2juncs, junc2bnds, true);

  gui::MakeCanvas(gui::CreateRandomColorTableWithSize(segs.second)(segs.first))
      .color(gui::Red)
      .thickness(2)
      .add(bndpixels)
      .show();
}