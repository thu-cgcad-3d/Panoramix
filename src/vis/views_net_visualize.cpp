#include "views_net_visualize.hpp"

namespace panoramix {
    namespace vis {

        using namespace core;

        Visualizer2D operator << (Visualizer2D viz, const ViewsNet::VertData & vd) {
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

        Visualizer2D operator << (Visualizer2D viz, const ViewsNet::GlobalData & netgb) {
            
            return viz;
        }

        Visualizer3D operator << (Visualizer3D viz, const core::ViewsNet::GlobalData & netgb) {
            viz << netgb.spatialLineSegments;
            return viz;
        }

        Visualizer2D operator << (Visualizer2D viz, const ViewsNet & net) {

            return viz;
        }

    }
}