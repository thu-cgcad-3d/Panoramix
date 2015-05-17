#ifndef PANORAMIX_CORE_SINGLE_VIEW_HPP
#define PANORAMIX_CORE_SINGLE_VIEW_HPP

#include "basic_types.hpp"
#include "cameras.hpp"

namespace panoramix {
    namespace core {

        std::vector<int> ComputeSpatialRegionProperties(
            const Imagei & segmentedRegions, const PanoramicCamera & cam,
            std::vector<std::vector<std::vector<Vec3>>> * ncontoursPtr = nullptr,
            std::vector<Vec3> * ncentersPtr = nullptr,
            std::vector<double> * areasPtr = nullptr);

        View<PartialPanoramicCamera, Imageub> PerfectRegionMaskView(
            const std::vector<std::vector<Vec3>> & contours,
            const Vec3 & center, double focal);





    }
}


#endif