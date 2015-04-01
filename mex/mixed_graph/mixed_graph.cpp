#include "mixed_graph.hpp"

namespace panoramix {
    namespace experimental {


        MixedGraph::MixedGraph(const Image & image,
            const HandledTable<RegionHandle, Vec7> & gcResp,
            const HandledTable<RegionBoundaryHandle, bool> & occlusionResp,
            const Point2 & cameraCenterPosition,
            double cameraFocal) : gcResponse(gcResp), occlusionResponse(occlusionResp) {

            core::LineSegmentExtractor lineExtractor;
            lineExtractor.params().algorithm = core::LineSegmentExtractor::LSD;
            lines = core::ClassifyEachAs(lineExtractor(image, 2), -1);

            view.camera = core::PerspectiveCamera(image.cols, image.rows, cameraCenterPosition, cameraFocal,
                core::Point3(0, 0, 0), core::Point3(1, 0, 0), core::Point3(0, 0, -1));
            view.image = image;

            // classify lines
            vps = core::EstimateVanishingPointsAndClassifyLines(view.camera, lines);
            // remove non-manhattan vps
            vps.erase(vps.begin() + 3, vps.end());
            for (auto & l : lines){
                if (l.claz >= 3){
                    l.claz = -1;
                }
            }

            // append lines
            core::AppendLines(regionLineGraph, lines, view.camera, vps);

            // append regions
            core::SegmentationExtractor segmenter;
            segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
            segmenter.params().c = 100.0;
            std::vector<core::Line2> pureLines(lines.size());
            for (int i = 0; i < lines.size(); i++)
                pureLines[i] = lines[i].component;
            int segmentsNum = 0;
            std::tie(regions, segmentsNum) = segmenter(image, pureLines, image.cols / 100.0);

            regionIds2Handles = core::AppendRegions(regionLineGraph, regions, view.camera, 0.001, 0.001, 3, 1);


            // bnd graph


        }



        void MixedGraph::showSegmentations() const{

        }

        void MixedGraph::showVanishingPointsAndLines() const{
            
        }

        void MixedGraph::showRegionOrientationConstriants() const{

        }

        HandledTable<RegionHandle, Plane3> MixedGraph::solve() const {

            HandledTable<RegionHandle, Plane3> results;





            return results;

        }

    }
}