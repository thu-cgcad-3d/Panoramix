#include "cameras.hpp"
#include "feature.hpp"
#include "utility.hpp"
#include "canvas.hpp"
#include "gui_util.hpp"

#include "../panoramix.unittest.hpp"

using namespace pano;
using namespace test;

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