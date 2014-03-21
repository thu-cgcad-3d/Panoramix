#include "views_net_visualize.hpp"

namespace panoramix {
    namespace core {

        ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const ViewsNet::VertData & vd) {
            static const Color colors[] = { Color(255, 0, 0), Color(0, 255, 0), Color(0, 0, 255) };
            viz.image() = vd.image;
            viz.params.thickness = 2;
            for (size_t i = 0; i < vd.lineSegments.size(); ++i){
                viz.params.color = vd.lineSegmentClasses[i] < 0 ? Color(255, 255, 255) : colors[vd.lineSegmentClasses[i]];
                viz << vd.lineSegments[i];
            }
            viz.params.thickness = 1;
            viz << vd.lineSegmentIntersections;
            return viz;
        }

        ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const ViewsNet::GlobalData & netgb) {
            
            return viz;
        }

        ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const ViewsNet & net) {

            return viz;
        }

    }
}