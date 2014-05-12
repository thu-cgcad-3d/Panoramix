#ifndef PANORAMIX_VIS_VIEWS_NET_VISUALIZE_HPP
#define PANORAMIX_VIS_VIEWS_NET_VISUALIZE_HPP

#include "../vis/visualize2d.hpp"
#include "../vis/visualize3d.hpp"
#include "views_net.hpp"
 
namespace panoramix {
    namespace rec {

        using vis::Visualizer2D;
        using vis::Visualizer3D;

        // 2d vis for vertdata
        Visualizer2D operator << (Visualizer2D viz, const ViewsNet::VertData & vd);


        // 3d vis for globaldata
        Visualizer3D operator << (Visualizer3D viz, const ViewsNet::GlobalData & netgb);

        // 2d for whole net
        Visualizer2D operator << (Visualizer2D viz, const ViewsNet & net);

    }
}
 
#endif