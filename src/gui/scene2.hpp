#pragma once

#include "../core/cameras.hpp"
#include "../core/forest.hpp"
#include "../core/meta.hpp"
#include "../core/utility.hpp"

#include "color.hpp"
#include "discretization.hpp"
#include "resource.hpp"
#include "shader.hpp"

namespace pano {
namespace gui2 {

using namespace pano::core;
using namespace pano::gui;

class Geometry {
public:
  virtual ~Geometry() {}
};

class RenderableObject {
public:
  virtual ~RenderableObject() {}

public:
  Mat4f modelMat;
};

class Light {
public:
  enum Shape { PointLight, ParallelLight, ConicLight };

public:
  Vec4f position;
  Vec3f direction;
  Vec4f color;
  Shape shape;
  bool enabled;
};

class Scene {
public:
  Scene() {}

public:
  std::vector<std::unique_ptr<RenderableObject>> objects;
  std::vector<std::unique_ptr<Light>> lights;
  std::vector<std::unique_ptr<PerspectiveCamera>> cameras;
};
}
}
