#ifndef PANORAMIX_VIS_REGIONS_NET_VISUALIZE_HPP
#define PANORAMIX_VIS_REGIONS_NET_VISUALIZE_HPP

#include "visualize2d.hpp"
#include "../core/regions_net.hpp"
 
namespace panoramix {
    namespace vis {

        Visualizer2D operator << (Visualizer2D viz, const core::RegionsNet & net);

    }
}
 
#endif