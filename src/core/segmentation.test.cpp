#include "cameras.hpp"
#include "feature.hpp"
#include "utility.hpp"
#include "../gui/canvas.hpp"
#include "../gui/gui_util.hpp"

#include "../experimental/pi_graph.hpp"

#include "../panoramix.unittest.hpp"

using namespace pano;
using namespace test;

TEST(SegmentationTest, SegmentationExtractor) {
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

TEST(SegmentationTest, SegmentationBoundaryJunction) {
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

TEST(SegmentationTest, SegmentationExtractorInPanorama) {
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
      "H:\\GitHub\\write-papers\\papers\\a\\figure\\supp\\";

  double alpha = 0.7;
  double beta = 1 - alpha;

  core::Image3ub imsegs =
      gui::CreateRandomColorTableWithSize(nsegs, gui::Transparent)(segs);
  imsegs = imsegs * alpha + im * beta;

  //gui::AsCanvas(im).saveAs(folder + "pano.png");
  //gui::AsCanvas(core::MakeView(im, cam).sampled(cam1).image)
  //    .show()
  //    .saveAs(folder + "pano_polar1.png");
  //gui::AsCanvas(core::MakeView(im, cam).sampled(cam2).image)
  //    .show()
  //    .saveAs(folder + "pano_polar2.png");

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

TEST(SegmentationTest, RemoveSmallRegionInSegmentation) {
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

TEST(SegmentationTest, RemoveThinRegionInSegmentation) {
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

TEST(SegmentationTest, ExtractSegmentationTopology) {
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