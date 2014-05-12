#include "views_net_visualize.hpp"

#include "../core/utilities.hpp"

namespace panoramix {
    namespace rec {

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

        Visualizer3D operator << (Visualizer3D viz, const ViewsNet::GlobalData & netgb) {
            std::vector<Line3> consLines;
            std::vector<Point3> consPoints;
            consLines.reserve(netgb.constraints.size());
            consPoints.reserve(netgb.constraints.size());
            for (auto & cons : netgb.refinedConstraints){
                auto & line1 = netgb.mergedSpatialLineSegments[cons.mergedSpatialLineSegmentIds[0]].component;
                auto & line2 = netgb.mergedSpatialLineSegments[cons.mergedSpatialLineSegmentIds[1]].component;
                auto pp = DistanceBetweenTwoLines(line1, line2);
                consLines.push_back(Line3(pp.second.first.position, pp.second.second.position));
                consPoints.push_back(cons.position);
            }
            viz << vis::manip3d::SetDefaultColor(ColorTag::Yellow)
                << vis::manip3d::SetPointSize(10.0)
                << vis::manip3d::SetLineWidth(1.0f)
                << consLines;
            return viz
                << vis::manip3d::SetDefaultColor(ColorTag::Black)
                << vis::manip3d::SetColorTableDescriptor(core::ColorTableDescriptor::RGB)
                << vis::manip3d::SetLineWidth(2.0f)
                << netgb.mergedSpatialLineSegments;
        }

        Visualizer2D operator << (Visualizer2D viz, const ViewsNet & net) {

            return viz;
        }

    }
}