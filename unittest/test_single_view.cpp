#include "../src/core/cameras.hpp"
#include "../src/core/feature.hpp"
#include "../src/core/single_view.hpp"
#include "../src/core/utility.hpp"
#include "../src/gui/canvas.hpp"
#include "../src/gui/utility.hpp"

#include "config.hpp"

using namespace pano;
using namespace core;
using namespace test;

TEST(SingleView, ComputeSpatialRegionProperties) {

  Image3ub im = gui::PickAnImage();
  SegmentationExtractor segmenter;
  segmenter.params().algorithm = SegmentationExtractor::GraphCut;

  Imagei segs;
  int nsegs;
  std::tie(segs, nsegs) = segmenter(im);

  auto result = CreatePerspectiveView(im);
  if (result.failed()) {
    std::cout << "failed creating perspective view from this image"
              << std::endl;
    return;
  }

  auto view = result.unwrap();

  std::vector<std::vector<std::vector<Vec3>>> contours;
  ComputeSpatialRegionProperties(segs, view.camera, &contours);
}
