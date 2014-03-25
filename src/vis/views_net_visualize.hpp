#ifndef PANORAMIX_VIS_VIEWS_NET_VISUALIZE_HPP
#define PANORAMIX_VIS_VIEWS_NET_VISUALIZE_HPP

#include "feature_visualize.hpp"
#include "../core/views_net.hpp"
 
namespace panoramix {
    namespace vis {

        ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const core::ViewsNet::VertData & vd);

        ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const core::ViewsNet::GlobalData & netgb);

        ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const core::ViewsNet & net);

    }
}
 
#endif