#ifndef PANORAMIX_REC_REGIONS_NET_VISUALIZE_HPP
#define PANORAMIX_REC_REGIONS_NET_VISUALIZE_HPP

#include "../vis/visualize2d.hpp"
#include "regions_net.hpp"
 
namespace panoramix {
    namespace rec {

        using vis::Visualizer2D;

        Visualizer2D operator << (Visualizer2D viz, const RegionsNet & net);

    }
}
 
#endif