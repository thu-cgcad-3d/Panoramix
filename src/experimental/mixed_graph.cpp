#include "../core/utilities.hpp"
#include "../core/containers.hpp"
#include "../ml/factor_graph.hpp"
#include "../gui/visualize2d.hpp"
#include "../gui/visualizers.hpp"
#include "mixed_graph.hpp"

namespace panoramix {
    namespace experimental {


        MixedGraph::MixedGraph(const Image & image,
            const Point2 & cameraCenterPosition,
            double cameraFocal) {

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
            segmenter.params().sigma = 0.5;
            segmenter.params().c = 30.0;
            segmenter.params().superpixelSizeSuggestion = -1;
            segmenter.params().superpixelNumberSuggestion = 500;
            std::vector<core::Line2> pureLines(lines.size());
            for (int i = 0; i < lines.size(); i++)
                pureLines[i] = lines[i].component;
            int segmentsNum = 0;
            std::tie(regions, segmentsNum) = segmenter(image, pureLines, image.cols / 100.0);

            regionIds2Handles = core::AppendRegions(regionLineGraph, regions, view.camera, 0.001, 0.001, 3, 1);
            

            // boundary junctions
            // add region boundary handle
            auto bjunctions = ExtractBoundaryJunctions(regions);
            for (auto & junc : bjunctions){
                auto & regionIds = junc.first;
                auto & ps = junc.second;
                std::set<RegionHandle> rhs;
                for (int regionId : regionIds){
                    RegionHandle rh = regionIds2Handles[regionId];
                    if (rh.invalid())
                        continue;
                    rhs.insert(rh);
                }
                if (rhs.size() < 3)
                    continue;
                // locate boundary handles among rhs
                std::set<RegionBoundaryHandle> bhs;
                for (RegionHandle rh : rhs){
                    auto & relatedBhs = regionLineGraph.topo(rh).constraints<RegionBoundaryData>();
                    for (RegionBoundaryHandle bh : relatedBhs){
                        auto anotherRh = regionLineGraph.topo(bh).component<0>();
                        if (anotherRh == rh){
                            anotherRh = regionLineGraph.topo(bh).component<1>();
                        }
                        if (Contains(rhs, anotherRh)){
                            bhs.insert(bh);
                        }
                    }
                }
                std::vector<Vec3> positions(ps.size());
                for (int i = 0; i < ps.size(); i++){
                    positions[i] = normalize(view.camera.spatialDirection(ps[i]));
                }
                boundaryJunctions.push_back(BoundaryJunction{ std::move(bhs), std::move(positions) });
            }

        }

        void MixedGraph::installGCResponse(const Imaged7 & gc){
            assert(gc.size() == regions.size());
            gcResponse = regionLineGraph.createComponentTable<RegionData>(Vec7());
            auto regionAreas = regionLineGraph.createComponentTable<RegionData>(0.0);
            for (auto it = gc.begin(); it != gc.end(); ++it){
                auto p = it.pos();
                const Vec7 & v = gc(p);
                RegionHandle rh(regions(p));
                if (rh.invalid())
                    continue;
                gcResponse[rh] += v;
                regionAreas[rh] += 1.0;
            }
            for (int i = 0; i < gcResponse.data.size(); i++){
                gcResponse.data[i] /= regionAreas.data[i];
                double sum = std::accumulate(
                    std::begin(gcResponse.data[i].val), 
                    std::end(gcResponse.data[i].val), 
                    0.0);
                assert(IsFuzzyZero(sum - 1.0, 0.1));
            }
        }


        double Distance(const Point2 & p, const std::vector<PixelLoc> & edge){
            double mind = std::numeric_limits<double>::max();
            for (const auto & c : edge){
                auto cc = core::vec_cast<double>(c);
                double d = core::Distance(p, cc);
                if (d < mind){
                    mind = d;
                }
            }
            return mind;
        }


        void MixedGraph::installOcclusionResponce(const std::vector<std::vector<PixelLoc>> & edges,
            const std::vector<double> & scores){
            assert(edges.size() == scores.size());
            occlusionResponse = regionLineGraph.createConstraintTable<RegionBoundaryData>(0.0);
            auto getBB = [&edges](int eid)->Box2 {
                return BoundingBoxOfContainer(edges.at(eid));
            };
            std::vector<int> eids(edges.size());
            std::iota(eids.begin(), eids.end(), 0);
            RTreeWrapper<int, decltype(getBB)> rtree(eids.begin(), eids.end(), getBB);
            for (auto & bd : regionLineGraph.constraints<RegionBoundaryData>()){
                double & resp = occlusionResponse[bd.topo.hd];
                double allSampleScoreSum = 0;
                int allSampleNum = 0;
                for (auto & ps : bd.data.normalizedEdges){
                    for (auto & p : ps){
                        allSampleNum++;
                        auto sample = view.camera.screenProjection(p);
                        double sampleScore = 0.0;
                        int nearbyEdgeNum = 0;
                        rtree.search(Box2(sample, sample).expand(5.0), 
                            [&scores, &sampleScore, &nearbyEdgeNum, &edges, &sample](int eid) -> bool {
                            double dist = Distance(sample, edges[eid]);
                            if (dist <= 3.0){
                                sampleScore += scores[eid] > 0 ? 1.0 : (scores[eid] == 0.0 ? 0.0 : -1.0);
                                nearbyEdgeNum++;
                            }
                            return true;
                        });
                        sampleScore /= std::max(nearbyEdgeNum, 1);
                        allSampleScoreSum += sampleScore;
                    }
                }
                occlusionResponse[bd.topo.hd] = allSampleScoreSum / std::max(allSampleNum, 1);
            }
        }

        void MixedGraph::installOcclusionResponce(const std::vector<std::vector<int>> & edges,
            const std::vector<double> & scores){
            assert(edges.size() == scores.size());
            std::vector<std::vector<PixelLoc>> pixels(edges.size());
            for (int i = 0; i < edges.size(); i++){
                auto & ps = pixels[i];
                ps.reserve(edges[i].size());
                for (int id : edges[i]){
                    ps.push_back(DecodeIndexToSubscript(id, regions.size()));
                }
            }
            installOcclusionResponce(pixels, scores);
        }

        void MixedGraph::showOcclusionResponse() const{
            gui::Visualizer2D vis(view.image);
            for (auto & b : regionLineGraph.constraints<RegionBoundaryData>()){
                double score = occlusionResponse[b.topo.hd];
                if (score <= 0)
                    continue;
                vis.params.color = gui::Red;
                vis.params.thickness = 1;
                for (auto & edge : b.data.normalizedEdges){
                    for (int i = 0; i < edge.size() - 1; i++){
                        vis << Line2(view.camera.screenProjection(edge[i]), view.camera.screenProjection(edge[i + 1]));
                    }
                }
            }
            vis << gui::manip2d::Show();
        }

        void MixedGraph::showSegmentations() const{
            int regionsNum = MinMaxValOfImage(regions).second + 1;
            auto ctable = gui::CreateRandomColorTableWithSize(regionsNum);
            auto segs = ctable(regions);
            gui::Visualizer2D(segs) << gui::manip2d::Show();
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