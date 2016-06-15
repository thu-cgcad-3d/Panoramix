#pragma once

#include "basic_types.hpp"
#include "cameras.hpp"
#include "utility.hpp"

#include "color.hpp"
#include "shader.hpp"

namespace pano {
namespace experimental {

using namespace pano::core;

class PiEngine {
public:
  class NodeBase {
  public:
    virtual DenseMatd matP() const = 0; // returns [3 x nparams]
    inline int nparams() const { return matP().cols; }
  };
  class EdgeBase {
  public:
    virtual DenseMatd invDirs() const = 0; // returns [k x 3]
    inline int ndirections() const { return invDirs().rows; }
    int vert1, vert2;
    double weight;
  };

public:
  PiEngine();
  PiEngine(const PanoramicCamera &cam, const Imagei &segs,
           const std::vector<Classified<Line3>> &lines,
           const std::vector<Vec3> &vps);
  PiEngine(const PerspectiveCamera &cam, const Imagei &segs,
           const std::vector<Classified<Line3>> &lines,
           const std::vector<Vec3> &vps);

  std::vector<Plane3> reconstruct() const;

public:
  std::vector<std::unique_ptr<NodeBase>> verts;
  std::vector<std::unique_ptr<EdgeBase>> edges;
};
}
}