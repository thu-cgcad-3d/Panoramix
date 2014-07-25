#include "reconstruction_engine_visualize.hpp"

#include "../core/utilities.hpp"
#include "regions_net_visualize.hpp"

namespace panoramix {
    namespace rec {

        using namespace core;

        Visualizer2D operator << (Visualizer2D viz, const ReconstructionEngine::ViewData & vd) {
            viz.setImage(vd.image);
            
            if (vd.regionNet){
                viz = viz << *vd.regionNet;
            }
            if (vd.lineNet){
                viz.params.thickness = 2;
                viz.params.colorTableDescriptor = vis::ColorTableDescriptor::RGB;
                viz = viz << vd.lineNet->lineSegments();

                viz.params.thickness = 1;
                viz.params.color = vis::ColorFromTag(vis::ColorTag::White);
                viz = viz << vd.lineNet->lineSegmentIntersections();
            }
            return viz;
        }




        Visualizer3D operator << (Visualizer3D viz, const ReconstructionEngine::GlobalData & netgb) {
            std::vector<Line3> consLines;
            std::vector<Point3> consPoints;
            /*consLines.reserve(netgb.constraints.size());
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
                << consLines;*/
           /* return viz
                << vis::manip3d::SetDefaultColor(vis::ColorTag::Black)
                << vis::manip3d::SetColorTableDescriptor(vis::ColorTableDescriptor::AllColors)
                << vis::manip3d::SetLineWidth(2.0f)
                << netgb.mergedSpatialLineSegmentsClassifiedWithStructureIds;*/
            return viz;
        }

        Visualizer2D operator << (Visualizer2D viz, const ReconstructionEngine & net) {

            return viz;
        }

    }
}