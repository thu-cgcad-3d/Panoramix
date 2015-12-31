#pragma once

#include "../core/cameras.hpp"
#include "../core/forest.hpp"
#include "../core/meta.hpp"
#include "../core/utility.hpp"

#include "color.hpp"
#include "shader.hpp"
#include "discretization.hpp"
#include "resource.hpp"

namespace pano {
namespace gui {
using namespace pano::core;

class Renderable {
public:
  ~Renderable();

protected:
  Mat4f _modelMat;
  TriMesh _mesh;
};

class Light {
public:
  enum Shape { PointLight, ParallelLight, ConicLight };

public:
  Vec4f position;
  Vec3f direction;
  Vec4f color;
  Shape shape;
};

class Scene {
public:
public:
  std::vector<Light> lights;
  std::vector<PerspectiveCamera> cams;
};
}
}
