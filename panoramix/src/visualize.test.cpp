#include "visualize.hpp"

#include "../panoramix.unittest.hpp"

using namespace pano;
using namespace pano::core;

TEST(Visualize, Box) {

  gui::VisMesh mesh;
  gui::VisMeshBuilder(mesh) << Sphere3(Point3(1, 1, 1), 1.5)
                            << Box3(Point3(2, 2, -1), Point3(1, 1, -2));

  //auto geo = std::make_shared<VisGeometry>()

}