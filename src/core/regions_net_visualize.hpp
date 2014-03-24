#ifndef PANORAMIX_CORE_REGIONS_NET_VISUALIZE_HPP
#define PANORAMIX_CORE_REGIONS_NET_VISUALIZE_HPP

#include "feature_visualize.hpp"
#include "regions_net.hpp"
 
namespace panoramix {
    namespace core {

        ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const RegionsNet & net);

    }
}
 
#endif