#include "../src/core/feature.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/algorithms.hpp"
#include "../src/core/mixed_graph.hpp"
#include "../src/vis/visualize2d.hpp"
#include "../src/vis/visualizers.hpp"

#include "config.hpp"

using namespace panoramix;
using namespace test;

void VisualizeMixedGraph(const core::Image & panorama,
    const core::MixedGraph & mg,
    const std::vector<core::MGPatch> & patches,
    const std::vector<core::Vec3> & vps,
    bool doModal){

    vis::Visualizer viz("mixed graph");
    viz.installingOptions.discretizeOptions.colorTable = vis::ColorTableDescriptor::RGB;

    struct UnaryID {
        int patchId;
        core::MGUnaryHandle uh;
    };

    std::vector<std::pair<UnaryID, vis::SpatialProjectedPolygon>> spps;
    std::vector<core::Classified<core::Line3>> lines;

    for (int i = 0; i < patches.size(); i++){
        auto & patch = patches[i];
        for (auto & uhv : patch.uhs){
            auto uh = uhv.first;
            auto & v = mg.data(uh);
            if (v.type == core::MGUnary::Region){
                auto & region = v;
                vis::SpatialProjectedPolygon spp;
                // filter corners
                core::ForeachCompatibleWithLastElement(region.normalizedCorners.begin(), region.normalizedCorners.end(),
                    std::back_inserter(spp.corners),
                    [](const core::Vec3 & a, const core::Vec3 & b) -> bool {
                    return core::AngleBetweenDirections(a, b) > M_PI / 500.0;
                });
                if (spp.corners.size() < 3)
                    continue;

                spp.projectionCenter = core::Point3(0, 0, 0);                
                spp.plane = uhv.second.interpretAsPlane();
                spps.emplace_back(UnaryID{ i, uh }, std::move(spp));
            }
            else if (v.type == core::MGUnary::Line){
                auto & line = v;
                lines.push_back(core::ClassifyAs(uhv.second.interpretAsLine(mg.data(uh), vps), mg.data(uh).lineClaz));
            }
        }
    }

    vis::ResourceStore::set("texture", panorama);

    auto sppCallbackFun = [&patches, &vps, &mg](vis::InteractionID iid, const std::pair<UnaryID, vis::SpatialProjectedPolygon> & spp) {
        std::cout << "uh: " << spp.first.uh.id << "in patch " << spp.first.patchId << std::endl;
    };

    viz.begin(spps, sppCallbackFun).shaderSource(vis::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
    viz.installingOptions.lineWidth = 4.0;
    viz.add(lines);

   /* std::vector<core::Line3> connectionLines;
    for (auto & patch : patches){
        for (auto & bhv : patch.bhs){
            auto bh = bhv.first;
            auto & v = bhv.second;
            auto & samples = mg.data(bh).normalizedAnchors;
            for (int i = 0; i < samples.size(); i++){
                connectionLines.emplace_back(normalize(samples[i]) * v.sampleDepthsOnRelatedUnaries.front()[i],
                    normalize(samples[i]) * v.sampleDepthsOnRelatedUnaries.back()[i]);
            }
        }
    }*/

    viz.installingOptions.discretizeOptions.color = vis::ColorTag::DarkGray;
    viz.installingOptions.lineWidth = 2.0;
    viz.renderOptions.renderMode = vis::RenderModeFlag::Triangles | vis::RenderModeFlag::Lines;
    //viz.add(connectionLines);
    viz.camera(core::PerspectiveCamera(800, 800, 500, { 1.0, 1.0, -1.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, -1.0 }));
    viz.show(doModal);

}

template <class UhColorizerFunT = core::ConstantFunctor<vis::Color>>
void ManuallyOptimizeMixedGraph(const core::Image & panorama,
    const core::MixedGraph & mg,
    core::MGPatch & patch,
    const std::vector<core::Vec3> & vps,
    UhColorizerFunT uhColorizer = UhColorizerFunT(vis::ColorTag::White),
    bool optimizeInEachIteration = true) {

    bool modified = true;
    auto sppCallbackFun = [&patch, &vps, &mg, &modified](vis::InteractionID iid,
        const std::pair<core::MGUnaryHandle, vis::Colored<vis::SpatialProjectedPolygon>> & spp) {
        std::cout << "uh: " << spp.first.id << std::endl;
        if (iid == vis::InteractionID::PressSpace){
            std::cout << "space pressed!" << std::endl;
            int & claz = patch.uhs[spp.first].claz;
            claz = (claz + 1) % vps.size();
            std::cout << "current orientation is : " << vps[claz] << std::endl;
            modified = true;
        }
    };

    while (modified){

        vis::ResourceStore::set("texture", panorama);

        modified = false;
        if (optimizeInEachIteration){
            core::MGPatchDepthsOptimizer(mg, patch, vps, false, core::MGPatchDepthsOptimizer::EigenSparseQR)
                .optimize();
        }
        patch /= core::AverageDepthOfPatch(patch);

        vis::Visualizer viz("mixed graph optimizable");
        viz.renderOptions.bwColor = 1.0;
        viz.renderOptions.bwTexColor = 0.0;
        viz.installingOptions.discretizeOptions.colorTable = vis::ColorTableDescriptor::RGB;
        std::vector<std::pair<core::MGUnaryHandle, vis::Colored<vis::SpatialProjectedPolygon>>> spps;
        std::vector<vis::Colored<core::Line3>> lines;

        for (auto & uhv : patch.uhs){
            auto uh = uhv.first;
            auto & v = mg.data(uh);
            if (v.type == core::MGUnary::Region){
                auto & region = v;
                vis::SpatialProjectedPolygon spp;
                // filter corners
                core::ForeachCompatibleWithLastElement(region.normalizedCorners.begin(), region.normalizedCorners.end(),
                    std::back_inserter(spp.corners),
                    [](const core::Vec3 & a, const core::Vec3 & b) -> bool {
                    return core::AngleBetweenDirections(a, b) > M_PI / 100.0;
                });
                if (spp.corners.size() < 3)
                    continue;

                spp.projectionCenter = core::Point3(0, 0, 0);
                spp.plane = core::PlaneOfMGUnary(region, vps, uhv.second);
                spps.emplace_back(uh, std::move(vis::ColorAs(spp, uhColorizer(uh))));
            }
            else if (v.type == core::MGUnary::Line){
                auto & line = v;
                lines.push_back(vis::ColorAs(core::LineOfMGUnary(line, vps, uhv.second), uhColorizer(uh)));
            }
        }

        viz.begin(spps, sppCallbackFun).shaderSource(vis::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
        viz.installingOptions.lineWidth = 4.0;
        viz.add(lines);

        std::vector<core::Line3> connectionLines;
        for (auto & bhv : patch.bhs){
            auto bh = bhv.first;
            auto & v = bhv.second;
            auto & samples = mg.data(bh).normalizedAnchors;
            for (int i = 0; i < samples.size(); i++){
                connectionLines.emplace_back(normalize(samples[i]) * v.sampleDepthsOnRelatedUnaries.front()[i],
                    normalize(samples[i]) * v.sampleDepthsOnRelatedUnaries.back()[i]);
            }
        }

        viz.installingOptions.discretizeOptions.color = vis::ColorTag::DarkGray;
        viz.installingOptions.lineWidth = 2.0;
        viz.renderOptions.renderMode = vis::RenderModeFlag::Triangles | vis::RenderModeFlag::Lines;
        viz.add(connectionLines);
        viz.camera(core::PerspectiveCamera(800, 800, 500, { 1.0, 1.0, -1.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, -1.0 }));
        viz.show(true, false);

        vis::ResourceStore::clear();
    }

}


TEST(MixedGraph, Build) {

    std::vector<core::View<core::PerspectiveCamera>> views;
    std::vector<core::RegionsGraph> regionsGraphs;
    std::vector<core::LinesGraph> linesGraphs;
    std::vector<core::Vec3> vanishingPoints;
    core::Image panorama;

    core::ComponentIndexHashMap<std::pair<core::RegionIndex, core::RegionIndex>, double> regionOverlappingsAcrossViews;
    core::ComponentIndexHashMap<std::pair<core::LineIndex, core::LineIndex>, core::Vec3> lineIncidencesAcrossViews;
    std::vector<std::map<std::pair<core::RegionHandle, core::LineHandle>, std::vector<core::Point2>>>
        regionLineConnectionsArray;

    core::LoadFromDisk("./cache/test_view.View.SampleViews.panorama", panorama);
    core::LoadFromDisk("./cache/test_view.View.SampleViews.views", views);
    core::LoadFromDisk("./cache/test_view.View.RegionsGraph.regionsGraphs", regionsGraphs);
    core::LoadFromDisk("./cache/test_view.View.LinesGraph.linesGraphs", linesGraphs);
    core::LoadFromDisk("./cache/test_view.View.LinesGraph.vanishingPoints", vanishingPoints);

    core::LoadFromDisk("./cache/test_view.View.ConstraintsAcrossViews.regionOverlappingsAcrossViews", regionOverlappingsAcrossViews);
    core::LoadFromDisk("./cache/test_view.View.ConstraintsAcrossViews.lineIncidencesAcrossViews", lineIncidencesAcrossViews);
    core::LoadFromDisk("./cache/test_view.View.RegionLineConnections.regionLineConnectionsArray", regionLineConnectionsArray);

    core::MGUnaryVarTable unaryVars;
    core::MGBinaryVarTable binaryVars;

    core::MixedGraph mg = core::BuildMixedGraph(views, regionsGraphs, linesGraphs,
        regionOverlappingsAcrossViews, lineIncidencesAcrossViews, regionLineConnectionsArray, vanishingPoints,
        unaryVars, binaryVars, 1.0);

    core::SaveToDisk("./cache/test_view.MixedGraph.Build.mg", mg);
    core::SaveToDisk("./cache/test_view.MixedGraph.Build.unaryVars", unaryVars);
    core::SaveToDisk("./cache/test_view.MixedGraph.Build.binaryVars", binaryVars);

}


TEST(MixedGraph, RebuildOneView){

    std::vector<core::View<core::PerspectiveCamera>> views;
    core::Image panorama;

    core::LoadFromDisk("./cache/test_view.View.SampleViews.panorama", panorama);
    core::LoadFromDisk("./cache/test_view.View.SampleViews.views", views);    
  
    for (int i = 0; i < views.size(); i++){
        core::MGUnaryVarTable unaryVars;
        core::MGBinaryVarTable binaryVars;
        std::vector<core::Vec3> vanishingPoints;

        core::MixedGraph mg = core::BuildMixedGraph({ views[i] }, vanishingPoints, unaryVars, binaryVars);
        core::SaveToDisk("./cache/test_view.MixedGraph.RebuildOneView.mg[" + std::to_string(i) + "]", mg);
        std::vector<core::MGPatch> naivePatches =
            core::SplitMixedGraphIntoPatches(mg, unaryVars, binaryVars);

        core::MGPatch & largestPatch = *std::max_element(naivePatches.begin(), naivePatches.end(),
            [](const core::MGPatch & a, const core::MGPatch & b){return a.uhs.size() < b.uhs.size(); });

        core::MGPatchDepthsOptimizer pdo(mg, largestPatch, vanishingPoints, false,
            core::MGPatchDepthsOptimizer::EigenSparseQR);
        pdo.optimize();

        VisualizeMixedGraph(panorama, mg, { largestPatch }, vanishingPoints, true);
    }


}



TEST(MixedGraph, BasicOptimization) {

    core::Image panorama;
    std::vector<core::Vec3> vanishingPoints;

    core::MixedGraph mg;
    core::MGUnaryVarTable unaryVars;
    core::MGBinaryVarTable binaryVars;

    core::LoadFromDisk("./cache/test_view.View.SampleViews.panorama", panorama);
    core::LoadFromDisk("./cache/test_view.View.LinesGraph.vanishingPoints", vanishingPoints);

    core::LoadFromDisk("./cache/test_view.MixedGraph.Build.mg", mg);
    core::LoadFromDisk("./cache/test_view.MixedGraph.Build.unaryVars", unaryVars);
    core::LoadFromDisk("./cache/test_view.MixedGraph.Build.binaryVars", binaryVars);

    {
        core::MGBinaryHandle bh(74);
        core::MGPatch patch = core::MakePatchOnBinary(mg, bh, unaryVars, binaryVars);

        core::MGPatchDepthsOptimizer pdo(mg, patch, vanishingPoints, false);

        double errBefore = core::AverageBinaryDistanceOfPatch(mg, patch, vanishingPoints);// / core::AverageDepthOfPatch(starPatch);
        bool feasible = pdo.optimize();
        double errAfter = core::AverageBinaryDistanceOfPatch(mg, patch, vanishingPoints);// / core::AverageDepthOfPatch(starPatch);

        EXPECT_TRUE(feasible);
        EXPECT_GE(errBefore, errAfter);
        if (errBefore < errAfter){
            auto t = mg.data(bh).type;
        }
        if (!feasible){
            auto t = mg.data(bh).type;
        }
    }

}


TEST(MixedGraph, NaiveHolisticOptimization) {

    core::Image panorama;
    std::vector<core::Vec3> vanishingPoints;

    core::MixedGraph mg;
    core::MGUnaryVarTable unaryVars;
    core::MGBinaryVarTable binaryVars;

    core::LoadFromDisk("./cache/test_view.View.SampleViews.panorama", panorama);
    core::LoadFromDisk("./cache/test_view.View.LinesGraph.vanishingPoints", vanishingPoints);

    core::LoadFromDisk("./cache/test_view.MixedGraph.Build.mg", mg);
    core::LoadFromDisk("./cache/test_view.MixedGraph.Build.unaryVars", unaryVars);
    core::LoadFromDisk("./cache/test_view.MixedGraph.Build.binaryVars", binaryVars);

    std::vector<core::MGPatch> naivePatches =
        core::SplitMixedGraphIntoPatches(mg, unaryVars, binaryVars);

    for (auto & patch : naivePatches){
        //VisualizeMixedGraph(panorama, mg, { patch }, vanishingPoints, true);
        core::MGPatchDepthsOptimizer pdo(mg, patch, vanishingPoints, false,
            core::MGPatchDepthsOptimizer::EigenSparseQRSimplified);

        {
            core::Clock clock("holistic optimization");
            double distBefore = core::AverageBinaryDistanceOfPatch(mg, patch, vanishingPoints)
                / core::AverageUnaryCenterDepthOfPatch(mg, patch, vanishingPoints);
            pdo.optimize();
            double distAfter = core::AverageBinaryDistanceOfPatch(mg, patch, vanishingPoints)
                / core::AverageUnaryCenterDepthOfPatch(mg, patch, vanishingPoints);

            EXPECT_GE(distBefore, distAfter);
        }

        VisualizeMixedGraph(panorama, mg, { patch }, vanishingPoints, true);
    }

    core::SaveToDisk("./cache/test_view.MixedGraph.NaiveHolisticOptimization.naivePatches", naivePatches);

}



TEST(MixedGraph, LinesOptimization) {

    core::Image panorama;
    std::vector<core::Vec3> vanishingPoints;

    core::MixedGraph mg;
    core::MGUnaryVarTable unaryVars;
    core::MGBinaryVarTable binaryVars;

    core::LoadFromDisk("./cache/test_view.View.SampleViews.panorama", panorama);
    core::LoadFromDisk("./cache/test_view.View.LinesGraph.vanishingPoints", vanishingPoints);

    core::LoadFromDisk("./cache/test_view.MixedGraph.Build.mg", mg);
    core::LoadFromDisk("./cache/test_view.MixedGraph.Build.unaryVars", unaryVars);
    core::LoadFromDisk("./cache/test_view.MixedGraph.Build.binaryVars", binaryVars);

    std::vector<core::MGPatch> naivePatches;
    core::LoadFromDisk("./cache/test_view.MixedGraph.NaiveHolisticOptimization.naivePatches", naivePatches);

    std::vector<core::MGPatch> linePatches;
    std::vector<core::MGPatch> lineMSTPatches;

    for (auto & patch : naivePatches){
        std::vector<core::MGPatch> subPatches = core::SplitPatch(mg, patch, [&mg](core::MGBinaryHandle bh){
            return mg.data(bh).type == core::MGBinary::LineLineIncidence ||
                mg.data(bh).type == core::MGBinary::LineLineIntersection;
        });
        for (auto & subp : subPatches){
            if (mg.data(subp.uhs.begin()->first).type == core::MGUnary::Line){
                linePatches.push_back(std::move(subp));
            }
        }
    }

    for (auto & patch : linePatches){
        lineMSTPatches.push_back(core::MinimumSpanningTreePatch(mg, patch, vanishingPoints));
    }

    for (auto & patch : linePatches){
        core::MGPatchDepthsOptimizer pdo(mg, patch, vanishingPoints, false, core::MGPatchDepthsOptimizer::EigenSparseQR);
        pdo.optimize();
    }

    for (auto & patch : lineMSTPatches){
        core::MGPatchDepthsOptimizer pdo(mg, patch, vanishingPoints, false, core::MGPatchDepthsOptimizer::EigenSparseQR);
        pdo.optimize();
    }

    VisualizeMixedGraph(panorama, mg, linePatches, vanishingPoints, false);
    VisualizeMixedGraph(panorama, mg, lineMSTPatches, vanishingPoints, true);

    core::SaveToDisk("./cache/test_view.MixedGraph.LinesOptimization.linePatches", linePatches);
    core::SaveToDisk("./cache/test_view.MixedGraph.LinesOptimization.lineMSTPatches", lineMSTPatches);

}