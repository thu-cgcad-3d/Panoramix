#include "mesh.hpp"
#include "../panoramix.unittest.hpp"

using namespace pano::core;

TEST(Mesh, A) {

    Mesh<Point3> mesh;
    MakeIcosahedron(mesh);


}
