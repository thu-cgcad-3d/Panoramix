#ifndef PANORAMIX_VIS_VIEWS_NET_VISUALIZE_HPP
#define PANORAMIX_VIS_VIEWS_NET_VISUALIZE_HPP

#include "visualize2d.hpp"
#include "visualize3d.hpp"
#include "../core/views_net.hpp"
 
namespace panoramix {
    namespace vis {

        // 2d vis for vertdata
        Visualizer2D operator << (Visualizer2D viz, const core::ViewsNet::VertData & vd);

        // 3d vis for globaldata
        Visualizer3D operator << (Visualizer3D viz, const core::ViewsNet::GlobalData & netgb);

        // 2d for whole net
        Visualizer2D operator << (Visualizer2D viz, const core::ViewsNet & net);

    }
}
 
#endif