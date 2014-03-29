#include "views_net_visualize.hpp"

namespace panoramix {
    namespace vis {

        using namespace core;

        Visualizer2D operator << (Visualizer2D viz, const ViewsNet::VertData & vd) {
            viz.setImage(vd.image);
            
            viz.params.thickness = 2;
            viz.params.colorTableDescriptor = ColorTableDescriptor::RGB;
            viz << vd.lineSegments;

            viz.params.thickness = 1;
            viz.params.color = ColorFromTag(ColorTag::White);
            viz << vd.lineSegmentIntersections;
            
            return viz;
        }

        Visualizer3D operator << (Visualizer3D viz, const core::ViewsNet::GlobalData & netgb) {
            return viz
                << vis::manip3d::SetDefaultColor(core::ColorFromTag(core::ColorTag::Black))
                << vis::manip3d::SetColorTableDescriptor(core::ColorTableDescriptor::RGB)
                << vis::manip3d::SetLineWidth(2.0f)
                << netgb.spatialLineSegments;
        }

        Visualizer2D operator << (Visualizer2D viz, const ViewsNet & net) {

            return viz;
        }

    }
}