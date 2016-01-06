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

class RenderableObject {
public:
  ~RenderableObject();

public:
  Mat4f modelMat;
  TriMesh mesh;
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
  Scene() {}
  

public:
  std::vector<std::unique_ptr<RenderableObject>> objects;
  std::vector<std::unique_ptr<Light>> lights;
  std::unique_ptr<PerspectiveCamera> cam;
};
}
}
