#include "../gui/scene.hpp"
#include "../panoramix.unittest.hpp"

#include "mesh.hpp"
#include "mesh_util.hpp"

using namespace pano;
using namespace pano::core;

TEST(Mesh, Basic) {

  Mesh<Point3> mesh;
  MakeQuadFacedCube(mesh);
  SearchAndAddFaces(mesh);
  gui::SceneBuilder sb;
  Image3f tex(200, 100);
  
  sb.begin(mesh)
      .shaderSource(gui::OpenGLShaderSourceDescriptor::XLines)
      .end()
      .show(true, true);
}
