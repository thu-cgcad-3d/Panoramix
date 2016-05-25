#pragma once

#include "basic_types.hpp"
#include "cameras.hpp"

namespace pano {
namespace core {
std::vector<int> ComputeSpatialRegionProperties(
    const Imagei &segmentedRegions, const PerspectiveCamera &cam,
    std::vector<std::vector<std::vector<Vec3>>> *ncontoursPtr = nullptr,
    std::vector<Vec3> *ncentersPtr = nullptr,
    std::vector<double> *areasPtr = nullptr);

std::vector<int> ComputeSpatialRegionProperties(
    const Imagei &segmentedRegions, const PanoramicCamera &cam,
    std::vector<std::vector<std::vector<Vec3>>> *ncontoursPtr = nullptr,
    std::vector<Vec3> *ncentersPtr = nullptr,
    std::vector<double> *areasPtr = nullptr);

std::vector<int> ComputeSpatialRegionProperties(
    const Imagei &segmentedRegions, const PartialPanoramicCamera &cam,
    std::vector<std::vector<std::vector<Vec3>>> *ncontoursPtr = nullptr,
    std::vector<Vec3> *ncentersPtr = nullptr,
    std::vector<double> *areasPtr = nullptr);

View<PartialPanoramicCamera, Imageub>
PerfectRegionMaskView(const std::vector<std::vector<Vec3>> &contours,
                      const Vec3 &center, double focal);
inline View<PartialPanoramicCamera, Imageub>
PerfectRegionMaskView(const std::vector<Vec3> &contours, const Vec3 &center,
                      double focal) {
  return PerfectRegionMaskView(std::vector<std::vector<Vec3>>{contours}, center,
                               focal);
}
}
}
