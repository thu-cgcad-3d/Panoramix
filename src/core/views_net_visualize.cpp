#include "views_net_visualize.hpp"

namespace panoramix {
    namespace core {

        ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const ViewsNet::VertData & vd) {
            static const Color colors[] = { Color(255, 0, 0), Color(0, 255, 0), Color(0, 0, 255) };
            viz.setImage(vd.image);
            
            viz.params.thickness = 2;
            viz.params.colorTableDescriptor = ColorTableDescriptor::RGB;
            viz << vd.lineSegments;

            viz.params.thickness = 1;
            viz.params.color = ColorFromTag(ColorTag::White);
            viz << vd.lineSegmentIntersections;
            
            return viz;
        }

        ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const ViewsNet::GlobalData & netgb) {
            
            return viz;
        }

        ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const ViewsNet & net) {

            return viz;
        }

    }
}