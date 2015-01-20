#include "../src/core/feature.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/algorithms.hpp"
#include "../src/core/view.hpp"
#include "../src/vis/visualize2d.hpp"
#include "../src/vis/visualizers.hpp"

#include "config.hpp"

using namespace panoramix;
using namespace test;

TEST(View, SampleViews) {

    auto panorama = cv::imread(ProjectDataDirStrings::PanoramaIndoor + "/13.jpg");
    core::ResizeToMakeHeightUnder(panorama, 800);
    auto panoView = core::CreatePanoramicView(panorama);

    std::vector<core::PerspectiveCamera> cams = {
        core::PerspectiveCamera(700, 700, 300, { 0, 0, 0 }, { 1, 0, 0 }, { 0, 0, -1 }),
        core::PerspectiveCamera(700, 700, 300, { 0, 0, 0 }, { 0, 1, 0 }, { 0, 0, -1 }),
        core::PerspectiveCamera(700, 700, 300, { 0, 0, 0 }, { -1, 0, 0 }, { 0, 0, -1 }),
        core::PerspectiveCamera(700, 700, 300, { 0, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 }),
        core::PerspectiveCamera(700, 700, 300, { 0, 0, 0 }, { 0, 0, 1 }, { 1, 0, 0 }),
        core::PerspectiveCamera(700, 700, 300, { 0, 0, 0 }, { 0, 0, -1 }, { 1, 0, 0 })
    };
    
    std::vector<core::View<core::PerspectiveCamera>> views;
    core::SampleViews(panoView, cams.begin(), cams.end(), std::back_inserter(views));

    core::SaveToDisk("./cache/test_view.View.SampleViews.panorama", panorama);
    core::SaveToDisk("./cache/test_view.View.SampleViews.views", views);

}

TEST(View, LinesGraph) {

    std::vector<core::View<core::PerspectiveCamera>> views;
    core::LoadFromDisk("./cache/test_view.View.SampleViews.views", views);

    std::vector<core::Vec3> vanishingPoints;
    std::vector<core::LinesGraph> linesGraphs;
    core::EstimateVanishingPointsAndBuildLinesGraphs(views, vanishingPoints, linesGraphs, 
        core::LineSegmentExtractor(), 10, 100, 5, false, true);

    for (int i = 0; i < views.size(); i++){
        vis::Visualizer2D viz(views[i].image);
        viz << vis::manip2d::SetColorTable(vis::ColorTableDescriptor::RGB)
            << vis::manip2d::SetThickness(2);
        for (auto & l : linesGraphs[i].elements<0>()){
            viz << l.data.line;
        }
        viz << vis::manip2d::SetThickness(1);
        for (auto & r : linesGraphs[i].elements<1>()){
            auto & line1 = linesGraphs[i].data(r.topo.lowers.front()).line;
            auto & line2 = linesGraphs[i].data(r.topo.lowers.back()).line;
            auto p1 = core::DistanceFromPointToLine(r.data.relationCenter, line1.component).second.position;
            auto p2 = core::DistanceFromPointToLine(r.data.relationCenter, line2.component).second.position;
            if (r.data.type == core::LineRelationData::Incidence){
                viz << vis::manip2d::SetColor(vis::ColorTag::Black);
                EXPECT_LT(core::Distance(p1, r.data.relationCenter), 100);
                EXPECT_LT(core::Distance(p2, r.data.relationCenter), 100);
            }
            else {
                viz << vis::manip2d::SetColor(vis::ColorTag::Cyan);
                EXPECT_LT(core::Distance(p1, r.data.relationCenter), 100);
                EXPECT_LT(core::Distance(p2, r.data.relationCenter), 100);
            }            
            viz << core::Line2(p1, r.data.relationCenter) << core::Line2(p2, r.data.relationCenter);
        }
        viz << vis::manip2d::Show();
    }

    core::SaveToDisk("./cache/test_view.View.LinesGraph.linesGraphs", linesGraphs);
    core::SaveToDisk("./cache/test_view.View.LinesGraph.vanishingPoints", vanishingPoints);  


}

TEST(View, RegionsGraph) {

    std::vector<core::View<core::PerspectiveCamera>> views;
    std::vector<core::LinesGraph> linesGraphs;
    core::LoadFromDisk("./cache/test_view.View.SampleViews.views", views);
    core::LoadFromDisk("./cache/test_view.View.LinesGraph.linesGraphs", linesGraphs);

    core::SegmentationExtractor seg;

    std::vector<core::Imagei> segmentedRegionsArray;
    std::vector<core::RegionsGraph> regionsGraphs;

    for (int i = 0; i < views.size(); i++){
        auto & v = views[i];
        std::vector<core::Line2> lines;
        for (auto & ld : linesGraphs[i].elements<0>()){
            auto line = ld.data.line.component;
            line.first -= core::normalize(line.direction()) * 5.0;
            line.second += core::normalize(line.direction()) * 5.0;
            lines.push_back(line);
        }
        auto segmentedRegions = seg(v.image, lines).first;
        int samplePointsOnBoundariesSum = 0;
        auto regions = core::CreateRegionsGraph(segmentedRegions, 10, 3);
        for (auto & r : regions.elements<1>()){
            for (auto & ps : r.data.sampledPoints){
                samplePointsOnBoundariesSum += ps.size();
            }
        }
        regionsGraphs.push_back(std::move(regions));
        segmentedRegionsArray.push_back(segmentedRegions);
    }

    core::SaveToDisk("./cache/test_view.View.RegionsGraph.regionsGraphs", regionsGraphs);
    core::SaveToDisk("./cache/test_view.View.RegionsGraph.segmentedRegionsArray", segmentedRegionsArray);
}


TEST(View, RegionLineConnections){

    std::vector<core::View<core::PerspectiveCamera>> views;
    std::vector<core::Imagei> segmentedRegionsArray;
    std::vector<core::RegionsGraph> regionsGraphs;
    std::vector<core::Vec3> vanishingPoints;
    std::vector<core::LinesGraph> linesGraphs;

    core::LoadFromDisk("./cache/test_view.View.SampleViews.views", views);
    core::LoadFromDisk("./cache/test_view.View.RegionsGraph.segmentedRegionsArray", segmentedRegionsArray);
    core::LoadFromDisk("./cache/test_view.View.RegionsGraph.regionsGraphs", regionsGraphs);
    core::LoadFromDisk("./cache/test_view.View.LinesGraph.linesGraphs", linesGraphs);
    core::LoadFromDisk("./cache/test_view.View.LinesGraph.vanishingPoints", vanishingPoints);

    std::vector<std::map<std::pair<core::RegionHandle, core::LineHandle>, std::vector<core::Point2>>> 
        regionLineConnectionsArray(views.size());
    for (int i = 0; i < views.size(); i++){
        
        vis::Visualizer2D vis(views[i].image);

        // visualize connections between regions and lines
        std::unordered_map<int, vis::Visualizer2D> vizs;

        //vis::Visualizer2D viz(vd.data.regionNet->image);
        int height = views[i].image.rows;
        int width = views[i].image.cols;

        core::ImageWithType<core::Vec3b> coloredOutput(segmentedRegionsArray[i].size());
        vis::ColorTable colors = vis::CreateRandomColorTableWithSize(regionsGraphs[i].internalElements<0>().size());
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                auto & c = colors[segmentedRegionsArray[i](cv::Point(x, y))];
                coloredOutput(cv::Point(x, y)) = core::Vec3b(c.red(), c.green(), c.blue());
            }
        }
        vizs[i].setImage(views[i].image);
        vizs[i].params.alphaForNewImage = 0.5;
        vizs[i] << coloredOutput;
        vizs[i] << vis::manip2d::SetColorTable(vis::ColorTableDescriptor::RGB);

        regionLineConnectionsArray[i] = core::RecognizeRegionLineConnections(segmentedRegionsArray[i], linesGraphs[i], 5);

        for (auto & l : regionLineConnectionsArray[i]) {
            auto & rh = l.first.first;
            auto & lh = l.first.second;
            auto & cline2 = linesGraphs[i].data(lh).line;
            auto & cam = views[i].camera;
            auto & viz = vizs[i];

            viz << vis::manip2d::SetColorTable(vis::ColorTableDescriptor::RGB) 
                << vis::manip2d::SetThickness(1) 
                << cline2;
            viz << vis::manip2d::SetColor(vis::ColorTag::Black)
                << vis::manip2d::SetThickness(1);
            auto & regionCenter = regionsGraphs[i].data(rh).center;
            for (auto & p : l.second) {
                viz << core::Line2(regionCenter, p);
            }
        }

        for (auto & viz : vizs) {
            viz.second << vis::manip2d::Show();
        }
    }

    core::SaveToDisk("./cache/test_view.View.RegionLineConnections.regionLineConnectionsArray", regionLineConnectionsArray);
}


TEST(View, ConstraintsAcrossViews){

    std::vector<core::View<core::PerspectiveCamera>> views;
    std::vector<core::Imagei> segmentedRegionsArray;
    std::vector<core::RegionsGraph> regionsGraphs;
    std::vector<core::Vec3> vanishingPoints;
    std::vector<core::LinesGraph> linesGraphs;

    core::LoadFromDisk("./cache/test_view.View.SampleViews.views", views);
    core::LoadFromDisk("./cache/test_view.View.RegionsGraph.segmentedRegionsArray", segmentedRegionsArray);
    core::LoadFromDisk("./cache/test_view.View.RegionsGraph.regionsGraphs", regionsGraphs);
    core::LoadFromDisk("./cache/test_view.View.LinesGraph.linesGraphs", linesGraphs);
    core::LoadFromDisk("./cache/test_view.View.LinesGraph.vanishingPoints", vanishingPoints);

    auto regionOverlappingsAcrossViews =
        core::RecognizeRegionOverlappingsAcrossViews(views, regionsGraphs);
    auto lineIncidencesAcrossViews = 
        core::RecognizeLineIncidencesAcrossViews(views, linesGraphs, M_PI / 4, 20.0 / 300.0);

    core::SaveToDisk("./cache/test_view.View.ConstraintsAcrossViews.regionOverlappingsAcrossViews", regionOverlappingsAcrossViews);
    core::SaveToDisk("./cache/test_view.View.ConstraintsAcrossViews.lineIncidencesAcrossViews", lineIncidencesAcrossViews);

}

