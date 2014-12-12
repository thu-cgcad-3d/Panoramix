#include "../src/core/feature.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/algorithms.hpp"
#include "../src/core/view.hpp"
#include "../src/vis/visualize2d.hpp"
#include "../src/vis/visualize3d.hpp"

#include "test_config.hpp"

using namespace panoramix;
using namespace test;

TEST(View, SampleViews) {

    std::cout << "cuda enabled device count: " << cv::gpu::getCudaEnabledDeviceCount() << std::endl;

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

    for (auto & v : views){
        cv::imshow("perspective", v.image);
        cv::waitKey();
    }

    core::SaveToDisk("./cache/test_view.SampleViews.panorama", panorama);
    core::SaveToDisk("./cache/test_view.SampleViews.views", views);

}

TEST(View, LinesGraph) {

    std::vector<core::View<core::PerspectiveCamera>> views;
    core::LoadFromDisk("./cache/test_view.SampleViews.views", views);

    std::vector<core::Vec3> vanishingPoints;
    std::vector<core::LinesGraph> linesGraphs;
    core::EstimateVanishingPointsAndBuildLinesGraphs(views, vanishingPoints, linesGraphs, 10, 50, 4);

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

    core::SaveToDisk("./cache/test_view.LinesGraph.linesGraphs", linesGraphs);
    core::SaveToDisk("./cache/test_view.LinesGraph.vanishingPoints", vanishingPoints);

}

TEST(View, RegionsGraph) {

    std::vector<core::View<core::PerspectiveCamera>> views;
    std::vector<core::LinesGraph> linesGraphs;
    core::LoadFromDisk("./cache/test_view.SampleViews.views", views);
    core::LoadFromDisk("./cache/test_view.LinesGraph.linesGraphs", linesGraphs);

    core::SegmentationExtractor seg;

    std::vector<core::Imagei> segmentedRegionsArray;
    std::vector<core::RegionsGraph> regionsGraphs;

    for (int i = 0; i < views.size(); i++){
        auto & v = views[i];
        std::vector<core::Line2> lines;
        for (auto & ld : linesGraphs[i].elements<0>()){
            auto line = ld.data.line.component;
            line.first -= core::normalize(line.direction()) * 100.0;
            line.second += core::normalize(line.direction()) * 100.0;
            lines.push_back(line);
        }
        auto segmentedRegions = seg(v.image, lines).first;
        int samplePointsOnBoundariesSum = 0;
        auto regions = core::CreateRegionsGraph(segmentedRegions, 15, 3);
        for (auto & r : regions.elements<1>()){
            for (auto & ps : r.data.sampledPoints){
                samplePointsOnBoundariesSum += ps.size();
            }
        }
        int samplePointsOnBoundariesSum2 = 0;
        auto regions2 = core::CreateRegionsGraph(segmentedRegions, 20, 3);
        for (auto & r : regions2.elements<1>()){
            for (auto & ps : r.data.sampledPoints){
                samplePointsOnBoundariesSum2 += ps.size();
            }
        }
        EXPECT_LE(samplePointsOnBoundariesSum2, samplePointsOnBoundariesSum);
        regionsGraphs.push_back(std::move(regions));
        segmentedRegionsArray.push_back(segmentedRegions);
    }

    core::SaveToDisk("./cache/test_view.RegionsGraph.regionsGraphs", regionsGraphs);
    core::SaveToDisk("./cache/test_view.RegionsGraph.segmentedRegionsArray", segmentedRegionsArray);
}


TEST(View, RegionLineConnections){

    std::vector<core::View<core::PerspectiveCamera>> views;
    std::vector<core::Imagei> segmentedRegionsArray;
    std::vector<core::RegionsGraph> regionsGraphs;
    std::vector<core::Vec3> vanishingPoints;
    std::vector<core::LinesGraph> linesGraphs;

    core::LoadFromDisk("./cache/test_view.SampleViews.views", views);
    core::LoadFromDisk("./cache/test_view.RegionsGraph.segmentedRegionsArray", segmentedRegionsArray);
    core::LoadFromDisk("./cache/test_view.RegionsGraph.regionsGraphs", regionsGraphs);
    core::LoadFromDisk("./cache/test_view.LinesGraph.linesGraphs", linesGraphs);
    core::LoadFromDisk("./cache/test_view.LinesGraph.vanishingPoints", vanishingPoints);

    std::vector<std::map<std::pair<core::RegionHandle, core::LineHandle>, std::vector<core::Point2>>> 
        regionLineConnectionsArray(views.size());
    for (int i = 0; i < views.size(); i++){
        
        vis::Visualizer2D vis(views[i].image);

        // visualize connections between regions and lines
        std::unordered_map<int, vis::Visualizer2D, core::HandleHasher<core::AtLevel<0>>> vizs;

        //vis::Visualizer2D viz(vd.data.regionNet->image);
        int height = views[i].image.rows;
        int width = views[i].image.cols;

        core::ImageWithType<core::Vec3b> coloredOutput(segmentedRegionsArray[i].size());
        vis::ColorTable colors = vis::CreateRandomColorTableWithSize(regionsGraphs[i].internalElements<0>().size());
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                coloredOutput(cv::Point(x, y)) =
                    vis::ToVec3b(colors[segmentedRegionsArray[i](cv::Point(x, y))]);
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

    core::SaveToDisk("./cache/test_view.RegionLineConnections.regionLineConnectionsArray", regionLineConnectionsArray);
}


TEST(View, ConstraintsAcrossViews){

    std::vector<core::View<core::PerspectiveCamera>> views;
    std::vector<core::Imagei> segmentedRegionsArray;
    std::vector<core::RegionsGraph> regionsGraphs;
    std::vector<core::Vec3> vanishingPoints;
    std::vector<core::LinesGraph> linesGraphs;

    core::LoadFromDisk("./cache/test_view.SampleViews.views", views);
    core::LoadFromDisk("./cache/test_view.RegionsGraph.segmentedRegionsArray", segmentedRegionsArray);
    core::LoadFromDisk("./cache/test_view.RegionsGraph.regionsGraphs", regionsGraphs);
    core::LoadFromDisk("./cache/test_view.LinesGraph.linesGraphs", linesGraphs);
    core::LoadFromDisk("./cache/test_view.LinesGraph.vanishingPoints", vanishingPoints);

    auto regionOverlappingsAcrossViews =
        core::RecognizeRegionOverlappingsAcrossViews(views, regionsGraphs);
    auto lineIncidencesAcrossViews = 
        core::RecognizeLineIncidencesAcrossViews(views, linesGraphs, M_PI / 4, 1e-3);

    core::SaveToDisk("./cache/test_view.ConstraintsAcrossViews.regionOverlappingsAcrossViews", regionOverlappingsAcrossViews);
    core::SaveToDisk("./cache/test_view.ConstraintsAcrossViews.lineIncidencesAcrossViews", lineIncidencesAcrossViews);

}


void VisualizeMixedGraph(const core::Image & panorama, core::MixedGraph & mg, const std::vector<core::Vec3> & vps) {

    vis::Visualizer3D viz;
    viz << vis::manip3d::SetDefaultColorTable(vis::ColorTableDescriptor::RGB);
    std::vector<vis::SpatialProjectedPolygon> spps;
    std::vector<core::Classified<core::Line3>> lines;
    for (auto & v : mg.elements<0>()){
        if (v.data.is<core::MGUnaryRegion>()){
            auto & region = v.data.ref<core::MGUnaryRegion>();
            vis::SpatialProjectedPolygon spp;
            // filter corners
            core::ForeachCompatibleWithLastElement(region.normalizedCorners.begin(), region.normalizedCorners.end(), std::back_inserter(spp.corners), 
                [](const core::Vec3 & a, const core::Vec3 & b) -> bool {
                return core::AngleBetweenDirections(a, b) > M_PI / 100.0;
            });
            if (spp.corners.size() < 3)
                continue;
            spp.projectionCenter = core::Point3(0, 0, 0);
            spp.plane = core::PlaneOfMGUnary(region, vps);
            spps.push_back(std::move(spp));
        }
        else {
            auto & line = v.data.ref<core::MGUnaryLine>();
            lines.push_back(core::ClassifyAs(core::LineOfMGUnary(line, vps), line.claz));
        }
    }
    viz << vis::manip3d::Begin(spps) << vis::manip3d::SetTexture(panorama) << vis::manip3d::End;
    viz << lines;
    viz << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Black);

    std::vector<core::Line3> connectionsRR, connectionsRL;
    for (auto & c : mg.elements<1>()){
        if (c.data.is<core::MGBinaryRegionRegionBoundary>()){
            std::vector<core::Vec3> filteredSamples;
            auto & samples = c.data.ref<core::MGBinaryRegionRegionBoundary>().samples;
            core::ForeachCompatibleWithLastElement(samples.begin(), samples.end(), std::back_inserter(filteredSamples),
                [](const core::Vec3 & a, const core::Vec3 & b) -> bool {
                //return core::AngleBetweenDirections(a, b) > M_PI / 80.0;
                return true;
            });
            for (auto & dir : filteredSamples){
                core::Line3 connection = {
                    core::LocationOnMGUnary(dir, mg.data(c.topo.lowers.front()), vps),
                    core::LocationOnMGUnary(dir, mg.data(c.topo.lowers.back()), vps)
                };
                connectionsRR.push_back(connection);
            }
        }
        else if (c.data.is<core::MGBinaryRegionRegionOverlapping>()){

        }
        else if (c.data.is<core::MGBinaryRegionLineConnection>()){
            std::vector<core::Vec3> filteredSamples;
            auto & samples = c.data.ref<core::MGBinaryRegionLineConnection>().samples;
            core::ForeachCompatibleWithLastElement(samples.begin(), samples.end(), std::back_inserter(filteredSamples),
                [](const core::Vec3 & a, const core::Vec3 & b) -> bool {
                //return core::AngleBetweenDirections(a, b) > M_PI / 80.0;
                return true;
            });
            for (auto & dir : filteredSamples){
                core::Line3 connection = {
                    core::LocationOnMGUnary(dir, mg.data(c.topo.lowers.front()), vps),
                    core::LocationOnMGUnary(dir, mg.data(c.topo.lowers.back()), vps)
                };
                connectionsRL.push_back(connection);
            }
        }
        else if (c.data.is<core::MGBinaryLineLineIntersection>()){

        }
        else if (c.data.is<core::MGBinaryLineLineIncidence>()){

        }
    }
    viz << vis::manip3d::SetDefaultLineWidth(1.0);
    viz << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::DarkGray) << connectionsRR;
    viz << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Black) << connectionsRL;
    auto bb = core::BoundingBoxOfContainer(lines);
    viz << vis::manip3d::Show();

}


TEST(View, MixedGraph) {

    std::vector<core::View<core::PerspectiveCamera>> views;
    std::vector<core::RegionsGraph> regionsGraphs;
    std::vector<core::LinesGraph> linesGraphs;
    std::vector<core::Vec3> vanishingPoints;
    core::Image panorama;

    core::ComponentIndexHashMap<std::pair<core::RegionIndex, core::RegionIndex>, double> regionOverlappingsAcrossViews;
    core::ComponentIndexHashMap<std::pair<core::LineIndex, core::LineIndex>, core::Vec3> lineIncidencesAcrossViews;
    std::vector<std::map<std::pair<core::RegionHandle, core::LineHandle>, std::vector<core::Point2>>>
        regionLineConnectionsArray;

    core::LoadFromDisk("./cache/test_view.SampleViews.panorama", panorama);
    core::LoadFromDisk("./cache/test_view.SampleViews.views", views);
    core::LoadFromDisk("./cache/test_view.RegionsGraph.regionsGraphs", regionsGraphs);
    core::LoadFromDisk("./cache/test_view.LinesGraph.linesGraphs", linesGraphs);
    core::LoadFromDisk("./cache/test_view.LinesGraph.vanishingPoints", vanishingPoints);

    core::LoadFromDisk("./cache/test_view.ConstraintsAcrossViews.regionOverlappingsAcrossViews", regionOverlappingsAcrossViews);
    core::LoadFromDisk("./cache/test_view.ConstraintsAcrossViews.lineIncidencesAcrossViews", lineIncidencesAcrossViews);
    core::LoadFromDisk("./cache/test_view.RegionLineConnections.regionLineConnectionsArray", regionLineConnectionsArray);

    auto mg = core::BuildMixedGraph(views, regionsGraphs, linesGraphs, 
        regionOverlappingsAcrossViews, lineIncidencesAcrossViews, regionLineConnectionsArray);

    core::InitializeMixedGraph(mg, vanishingPoints);
    VisualizeMixedGraph(panorama, mg, vanishingPoints);

    std::unordered_map<core::MixedGraphUnaryHandle, int> ccids;
    int ccNum = core::MarkConnectedComponentIds(mg, ccids);
    std::cout << "ccNum = " << ccNum << std::endl;
    
    core::SolveDepthsInMixedGraphMOSEK(mg, vanishingPoints, ccids);
    VisualizeMixedGraph(panorama, mg, vanishingPoints);

}







int main(int argc, char * argv[], char * envp[]) {
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    testing::GTEST_FLAG(filter) = "View.MixedGraph";
    return RUN_ALL_TESTS();
}
