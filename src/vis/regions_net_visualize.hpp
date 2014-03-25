#ifndef PANORAMIX_VIS_REGIONS_NET_VISUALIZE_HPP
#define PANORAMIX_VIS_REGIONS_NET_VISUALIZE_HPP

#include "feature_visualize.hpp"
#include "../core/regions_net.hpp"
 
namespace panoramix {
    namespace vis {

        ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const core::RegionsNet & net);

    }
}
 
#endif