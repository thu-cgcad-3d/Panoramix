#include "../core/cameras.hpp"
#include "../core/feature.hpp"
#include "../core/single_view.hpp"
#include "../core/utility.hpp"
#include "../gui/canvas.hpp"
#include "../gui/utility.hpp"

#include "../panoramix.unittest.hpp"

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
