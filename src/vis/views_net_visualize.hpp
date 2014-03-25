#ifndef PANORAMIX_VIS_VIEWS_NET_VISUALIZE_HPP
#define PANORAMIX_VIS_VIEWS_NET_VISUALIZE_HPP

#include "visualize2d.hpp"
#include "../core/views_net.hpp"
 
namespace panoramix {
    namespace vis {

        Visualizer2D operator << (Visualizer2D viz, const core::ViewsNet::VertData & vd);

        Visualizer2D operator << (Visualizer2D viz, const core::ViewsNet::GlobalData & netgb);

        Visualizer2D operator << (Visualizer2D viz, const core::ViewsNet & net);

    }
}
 
#endif