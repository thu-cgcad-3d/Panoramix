#ifndef PANORAMIX_CORE_VIEWS_NET_VISUALIZE_HPP
#define PANORAMIX_CORE_VIEWS_NET_VISUALIZE_HPP

#include "feature_visualize.hpp"
#include "views_net.hpp"
 
namespace panoramix {
    namespace core {

        ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const ViewsNet::VertData & vd);

        ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const ViewsNet::GlobalData & netgb);

        ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const ViewsNet & net);

    }
}
 
#endif