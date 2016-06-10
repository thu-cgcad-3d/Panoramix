#include "pch.hpp"

#include "pi_engine.hpp"

namespace pano {
namespace experimental {
PiEngine::PiEngine() {}
PiEngine::PiEngine(const PanoramicCamera &cam, const Imagei &segs,
                   const std::vector<Classified<Line3>> &lines,
                   const std::vector<Vec3> &vps) {}
PiEngine::PiEngine(const PerspectiveCamera &cam, const Imagei &segs,
                   const std::vector<Classified<Line3>> &lines,
                   const std::vector<Vec3> &vps) {}
std::vector<Plane3> PiEngine::reconstruct() const {
  return std::vector<Plane3>();
}
}
}
