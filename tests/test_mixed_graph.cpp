#include "../src/core/feature.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/algorithms.hpp"
#include "../src/core/mixed_graph.hpp"
#include "../src/vis/visualize2d.hpp"
#include "../src/vis/visualizers.hpp"

#include "test_config.hpp"

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
                spp.plane = core::PlaneOfMGUnary(region, vps, uhv.second);
                spps.emplace_back(UnaryID{ i, uh }, std::move(spp));
            }
            else if (v.type == core::MGUnary::Line){
                auto & line = v;
                lines.push_back(core::ClassifyAs(core::LineOfMGUnary(line, vps, uhv.second), uhv.second.claz));
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

    std::vector<core::Line3> connectionLines;
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
    }

    viz.installingOptions.discretizeOptions.color = vis::ColorTag::DarkGray;
    viz.installingOptions.lineWidth = 2.0;
    viz.renderOptions.renderMode = vis::RenderModeFlag::Triangles | vis::RenderModeFlag::Lines;
    viz.add(connectionLines);
    viz.camera(core::PerspectiveCamera(800, 800, 500, { 1.0, 1.0, -1.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, -1.0 }));
    viz.show(doModal);

}

struct DefaultColorizer {
    vis::Color operator()(core::MGUnaryHandle uh) const {
        return vis::ColorTag::White;
    }
};

template <class UhColorizerFunT = DefaultColorizer>
void ManuallyOptimizeMixedGraph(const core::Image & panorama,
    const core::MixedGraph & mg,
    core::MGPatch & patch,
    const std::vector<core::Vec3> & vps,
    UhColorizerFunT uhColorizer = UhColorizerFunT(),
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

        viz.installingOptions.discretizeOptions.color = vis::ColorFromTag(vis::ColorTag::DarkGray);
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

    core::InitializeUnaryVarDepths(unaryVars, 1.0);
    {
        core::MGBinaryHandle bh(74);
        core::MGPatch patch = core::MakePatchOnBinary(mg, bh, unaryVars, binaryVars);

        core::MGPatchDepthsOptimizer pdo(mg, patch, vanishingPoints, false);

        double errBefore = core::AverageBinaryDistanceOfPatch(patch);// / core::AverageDepthOfPatch(starPatch);
        bool feasible = pdo.optimize();
        double errAfter = core::AverageBinaryDistanceOfPatch(patch);// / core::AverageDepthOfPatch(starPatch);

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


TEST(MixedGraph, OptimizateBinaryPatch) {

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

    // test small patches
    // mosek
    for (auto & b : mg.elements<1>()){
        core::InitializeUnaryVarDepths(unaryVars, 1.0);
        core::MGPatch patch = core::MakePatchOnBinary(mg, b.topo.hd, unaryVars, binaryVars);
        auto oldPatch = patch;
        core::MGPatchDepthsOptimizer pdo(mg, patch, vanishingPoints, false);

        auto & uh1 = patch.uhs.begin()->first;
        auto & uh2 = std::next(patch.uhs.begin())->first;
        auto & ud1 = mg.data(uh1);
        auto & ud2 = mg.data(uh2);
        auto & uvar1 = unaryVars[uh1];
        auto & uvar2 = unaryVars[uh2];

        auto & bd = mg.data(patch.bhs.begin()->first);

        double errBefore = core::AverageBinaryDistanceOfPatch(patch);// / core::AverageDepthOfPatch(patch);
        bool feasible = pdo.optimize();
        double errAfter = core::AverageBinaryDistanceOfPatch(patch);// / core::AverageDepthOfPatch(patch);

        ASSERT_TRUE(feasible);
        if (core::IsGoodBinary(mg, b.topo.hd, unaryVars, vanishingPoints)){
            EXPECT_TRUE(errAfter - errBefore < 1e-2);
        }
    }

    // eigen
    for (auto & b : mg.elements<1>()){
        core::InitializeUnaryVarDepths(unaryVars, 1.0);
        core::MGPatch patch = core::MakePatchOnBinary(mg, b.topo.hd, unaryVars, binaryVars);
        auto oldPatch = patch;
        core::MGPatchDepthsOptimizer pdo(mg, patch, vanishingPoints, false, core::MGPatchDepthsOptimizer::EigenSparseQR);

        auto & uh1 = patch.uhs.begin()->first;
        auto & uh2 = std::next(patch.uhs.begin())->first;
        auto & ud1 = mg.data(uh1);
        auto & ud2 = mg.data(uh2);
        auto & uvar1 = unaryVars[uh1];
        auto & uvar2 = unaryVars[uh2];

        auto & bd = mg.data(patch.bhs.begin()->first);

        double errBefore = core::AverageBinaryDistanceOfPatch(patch, 2);
        bool feasible = pdo.optimize();
        double errAfter = core::AverageBinaryDistanceOfPatch(patch, 2);

        ASSERT_TRUE(feasible);
        if (core::IsGoodBinary(mg, b.topo.hd, unaryVars, vanishingPoints)){
            EXPECT_TRUE(errAfter - errBefore < 1e-2);
        }
    }


}



TEST(MixedGraph, OptimizateStarPatch){

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

    for (auto & u : mg.elements<0>()){
        core::InitializeUnaryVarDepths(unaryVars, 1.0);
        core::MGPatch starPatch = core::MakeStarPatchAroundUnary(mg, u.topo.hd, unaryVars, binaryVars);
        core::MGPatchDepthsOptimizer pdo(mg, starPatch, vanishingPoints, false, core::MGPatchDepthsOptimizer::MosekLinearProgramming);

        double errBefore = core::AverageBinaryDistanceOfPatch(starPatch);
        bool feasible = pdo.optimize();
        double errAfter = core::AverageBinaryDistanceOfPatch(starPatch);

        ASSERT_TRUE(feasible);
        if (core::IsGoodPatch(mg, starPatch, vanishingPoints)){
            EXPECT_TRUE(errAfter - errBefore < 1e-2);
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

    core::InitializeUnaryVarDepths(unaryVars, 1.0);

    std::vector<core::MGPatch> naivePatches =
        core::SplitMixedGraphIntoPatches(mg, unaryVars, binaryVars);

    for (auto & patch : naivePatches){
        VisualizeMixedGraph(panorama, mg, { patch }, vanishingPoints, false);
        core::MGPatchDepthsOptimizer pdo(mg, patch, vanishingPoints, false,
            core::MGPatchDepthsOptimizer::EigenSparseQR);

        double distBefore = core::AverageBinaryDistanceOfPatch(patch) / core::AverageDepthOfPatch(patch);
        pdo.optimize();
        double distAfter = core::AverageBinaryDistanceOfPatch(patch) / core::AverageDepthOfPatch(patch);

        EXPECT_GE(distBefore, distAfter);
        VisualizeMixedGraph(panorama, mg, { patch }, vanishingPoints, true);
    }

    core::SaveToDisk("./cache/test_view.MixedGraph.NaiveHolisticOptimization.naivePatches", naivePatches);

}



TEST(MixedGraph, GoodPatchOptimization){

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

    std::vector<core::MGPatch> goodPatches;
    for (auto & p : naivePatches){
        auto gps = core::SplitIntoGoodPatches(mg, p, vanishingPoints);
        for (auto & gp : gps){
            goodPatches.push_back(core::MinimumSpanningTreePatch(mg, gp));
        }
    }

    for (int i = 0; i < goodPatches.size(); i++){
        ASSERT_TRUE(core::IsTreePatch(mg, goodPatches[i]));
        VisualizeMixedGraph(panorama, mg, { goodPatches[i] }, vanishingPoints, i == goodPatches.size() - 1);
    }


    for (auto & patch : goodPatches){
        core::MGPatchDepthsOptimizer pdo(mg, patch, vanishingPoints, false,
            core::MGPatchDepthsOptimizer::MosekLinearProgramming);
        double distBefore = core::AverageBinaryDistanceOfPatch(patch) / core::AverageDepthOfPatch(patch);
        pdo.optimize();
        double distAfter = core::AverageBinaryDistanceOfPatch(patch) / core::AverageDepthOfPatch(patch);
    }

    VisualizeMixedGraph(panorama, mg, naivePatches, vanishingPoints, false);
    VisualizeMixedGraph(panorama, mg, goodPatches, vanishingPoints, true);

    core::SaveToDisk("./cache/test_view.MixedGraph.GoodPatchOptimization.goodPatches", goodPatches);

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
        lineMSTPatches.push_back(core::MinimumSpanningTreePatch(mg, patch));
    }

    for (auto & patch : linePatches){
        core::MGPatchDepthsOptimizer pdo(mg, patch, vanishingPoints, false, core::MGPatchDepthsOptimizer::MosekLinearProgramming);
        pdo.optimize();
    }

    for (auto & patch : lineMSTPatches){
        core::MGPatchDepthsOptimizer pdo(mg, patch, vanishingPoints, false, core::MGPatchDepthsOptimizer::MosekLinearProgramming);
        pdo.optimize();
    }

    VisualizeMixedGraph(panorama, mg, linePatches, vanishingPoints, false);
    VisualizeMixedGraph(panorama, mg, lineMSTPatches, vanishingPoints, true);

    core::SaveToDisk("./cache/test_view.MixedGraph.LinesOptimization.linePatches", linePatches);
    core::SaveToDisk("./cache/test_view.MixedGraph.LinesOptimization.lineMSTPatches", lineMSTPatches);

}


TEST(MixedGraph, Reconstruct){

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

    using namespace core;
    auto & vps = vanishingPoints;

    std::vector<core::MGPatch> patches =
        core::SplitMixedGraphIntoPatches(mg, unaryVars, binaryVars);

    for (auto & pp : patches) {
        for (auto uhv : pp.uhs) {
            if (mg.data(uhv.first).type == MGUnary::Region)
                continue;

            core::MGPatch patchSeed = { { uhv }, {} };

            core::HandledTable<MGUnaryHandle, double> uhScores(mg.internalElements<0>().size(), -1.0);
            core::HandledTable<MGBinaryHandle, double> bhScores(mg.internalElements<1>().size(), -1.0);

            std::vector<MGUnaryHandle> orderedUhs; // the installation order
            orderedUhs.reserve(mg.internalElements<0>().size());
            orderedUhs.push_back(uhv.first);

            {
                core::Clock clock = "extend patch";
                ExtandPatch(mg, patchSeed, unaryVars, binaryVars, vps, uhScores, bhScores, orderedUhs, 1e-3,
                    core::StaticConstantFunctor<bool, true>());
                patchSeed.updateBinaryVars(mg, vps);
            }

            if (patchSeed.uhs.size() == 1)
                continue;


            core::SaveToDisk("./cache/tmp", patchSeed, uhScores, bhScores, orderedUhs);
            return;

            core::HandledTable<MGBinaryHandle, int> containedInPatchSeed(mg.internalElements<1>().size(), false);
            for (auto bhv : patchSeed.bhs){
                containedInPatchSeed[bhv.first] = true;
            }

            core::HandledTable<MGUnaryHandle, int> uhOrders(mg.internalElements<0>().size(), -1);
            for (int o = 0; o < orderedUhs.size(); o++){
                uhOrders[orderedUhs[o]] = o;
            }
            std::vector<std::vector<MGBinaryHandle>> binariesDeterminedByOrderedUhs(orderedUhs.size());

            // for each uh-i, collect all binaries connecting two uhs: 
            //      one uh installed before uh-i, 
            //      another installed after uh-i (or IS uh-i)
            // sum the violations of these binaries for uh-i
            // we can iterate each bh to install bh into its determiner uhs
            // # FIXME!! TOO SLOW!!
            // record ancestries
            //{
            //    core::Clock clock = "...";
            //    for (auto & b : mg.elements<1>()){
            //        auto bh = b.topo.hd;
            //        auto uh1 = mg.topo(bh).lowers[0];
            //        auto uh2 = mg.topo(bh).lowers[1];
            //        if (!Contains(patchSeed, uh1) || !Contains(patchSeed, uh2))
            //            continue;

            //        int uh1InstallOrd = uhOrders.at(uh1);
            //        int uh2InstallOrd = uhOrders.at(uh2);

            //        if (uh1InstallOrd > uh2InstallOrd){
            //            std::swap(uh1, uh2);
            //            std::swap(uh1InstallOrd, uh2InstallOrd);
            //        }

            //        std::unordered_set<MGUnaryHandle> mayBeRelatedUhs;
            //        for (int ord = uh1InstallOrd + 1; ord <= uh2InstallOrd; ord++){
            //            mayBeRelatedUhs.insert(orderedUhs[ord]);
            //        }

            //        auto getUhParentsInPatch = [&mg, &patchSeed, &uhOrders, &containedInPatchSeed,
            //            uh1InstallOrd, uh2InstallOrd](MGUnaryHandle uh){
            //            std::vector<MGUnaryHandle> parents;
            //            int uhOrd = uhOrders.at(uh);
            //            for (auto bh : mg.topo(uh).uppers){
            //                if (!containedInPatchSeed[bh]) // only conisder bhs in patch
            //                    continue;
            //                auto anotherUh = mg.topo(bh).lowers[0];
            //                if (anotherUh == uh)
            //                    anotherUh = mg.topo(bh).lowers[1];
            //                int anotherUhOrd = uhOrders.at(anotherUh);
            //                if (anotherUhOrd > uhOrd)
            //                    continue;
            //                if (anotherUhOrd <= uh1InstallOrd || anotherUhOrd > uh2InstallOrd)
            //                    continue;
            //                parents.push_back(anotherUh);
            //            }
            //            return parents;
            //        };
            //        DepthFirstSearch(mayBeRelatedUhs.begin(), mayBeRelatedUhs.end(), getUhParentsInPatch, 
            //            [&binariesDeterminedByOrderedUhs, &uhOrders, bh](MGUnaryHandle uh){
            //            binariesDeterminedByOrderedUhs[uhOrders.at(uh)].push_back(bh);
            //            return true;
            //        });
            //    }
            //}

            { // not accurate ? do we need the absolute accuracy?
                core::Clock clock = "...";
                for (auto & b : mg.elements<1>()){
                    auto bh = b.topo.hd;
                    auto uh1 = mg.topo(bh).lowers[0];
                    auto uh2 = mg.topo(bh).lowers[1];
                    if (!Contains(patchSeed, uh1) || !Contains(patchSeed, uh2))
                        continue;

                    int uh1InstallOrd = uhOrders.at(uh1);
                    int uh2InstallOrd = uhOrders.at(uh2);

                    if (uh1InstallOrd > uh2InstallOrd){
                        std::swap(uh1InstallOrd, uh2InstallOrd);
                    }
                    for (int ord = uh1InstallOrd + 1; ord <= uh2InstallOrd; ord++){
                        binariesDeterminedByOrderedUhs[ord].push_back(bh);
                    }
                }
            }



            // compute the confidence of each uh installation
            std::vector<double> installationConfidenceOfOrderedUhs(orderedUhs.size(), 0.0);
            double avgDepthOfPatch = AverageDepthOfPatch(patchSeed);
            for (int i = 0; i < orderedUhs.size(); i++){
                auto uh = orderedUhs[i];
                //double uhScore = uhScores.at(uh);
                auto & affectedBhs = binariesDeterminedByOrderedUhs[i];
                if (affectedBhs.size() == 0){
                    assert(uh == uhv.first);
                    installationConfidenceOfOrderedUhs[i] = 1e5;
                    continue;
                }

                double distanceSum = 0.0;
                for (auto bh : affectedBhs){
                    distanceSum += AnchorDistanceSumOnBinary(mg, bh, patchSeed, vanishingPoints);
                }
                distanceSum /= affectedBhs.size();
                distanceSum /= avgDepthOfPatch;
                if (distanceSum == 0.0){
                    distanceSum = 1e-5;
                }
                double uhConfidence = 1.0 - distanceSum;

                installationConfidenceOfOrderedUhs[i] = uhConfidence;
            }

            // check the most unconfident uhs
            std::vector<MGUnaryHandle> uhsOrderedByUnConfidencies = orderedUhs;
            std::sort(uhsOrderedByUnConfidencies.begin(), uhsOrderedByUnConfidencies.end(),
                [&uhOrders, &installationConfidenceOfOrderedUhs](MGUnaryHandle uh1, MGUnaryHandle uh2){
                return installationConfidenceOfOrderedUhs[uhOrders.at(uh1)] <
                    installationConfidenceOfOrderedUhs[uhOrders.at(uh2)];
            });

            ManuallyOptimizeMixedGraph(panorama, mg, patchSeed, vps,
                [&uhsOrderedByUnConfidencies, &uhOrders](MGUnaryHandle curUh){
                /*  auto pos = std::find(uhsOrderedByUnConfidencies.begin(), uhsOrderedByUnConfidencies.end(), curUh) -
                uhsOrderedByUnConfidencies.begin();
                double r = double(pos) / uhsOrderedByUnConfidencies.size();
                return vis::Color(r, r, r, 1.0);*/
                return curUh == uhsOrderedByUnConfidencies.front() ?
                    vis::ColorTag::Red : vis::ColorTag::Blue;
                /*double r = double(uhOrders[curUh]) / uhsOrderedByUnConfidencies.size();
                return vis::ColorFromHSV(r, 0.8, 0.8) / 255.0;*/
            });

        }
    }

}

TEST(MixedGraph, Test){
    using namespace core;
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

    core::MGPatch patchSeed;

    core::HandledTable<MGUnaryHandle, double> uhScores;
    core::HandledTable<MGBinaryHandle, double> bhScores;

    std::vector<MGUnaryHandle> orderedUhs;
    core::LoadFromDisk("./cache/tmp", patchSeed, uhScores, bhScores, orderedUhs);

    core::HandledTable<MGBinaryHandle, int> containedInPatchSeed(mg.internalElements<1>().size(), false);
    for (auto bhv : patchSeed.bhs){
        containedInPatchSeed[bhv.first] = true;
    }

    core::HandledTable<MGUnaryHandle, int> uhOrders(mg.internalElements<0>().size(), -1);
    for (int o = 0; o < orderedUhs.size(); o++){
        uhOrders[orderedUhs[o]] = o;
    }
    std::vector<std::vector<MGBinaryHandle>> binariesDeterminedByOrderedUhs(orderedUhs.size());

    { // not accurate ? do we need the absolute accuracy?
        core::Clock clock = "...";
        for (auto & b : mg.elements<1>()){
            auto bh = b.topo.hd;
            auto uh1 = mg.topo(bh).lowers[0];
            auto uh2 = mg.topo(bh).lowers[1];
            if (!Contains(patchSeed, uh1) || !Contains(patchSeed, uh2))
                continue;

            int uh1InstallOrd = uhOrders.at(uh1);
            int uh2InstallOrd = uhOrders.at(uh2);

            if (uh1InstallOrd > uh2InstallOrd){
                std::swap(uh1InstallOrd, uh2InstallOrd);
            }
            for (int ord = uh1InstallOrd + 1; ord <= uh2InstallOrd; ord++){
                if (mg.data(orderedUhs[ord]).type == MGUnary::Region){
                    binariesDeterminedByOrderedUhs[ord].push_back(bh);
                }
            }
        }
    }



    // compute the confidence of each uh installation
    std::vector<double> installationConfidenceOfOrderedUhs(orderedUhs.size(), 0.0);
    double avgDepthOfPatch = AverageDepthOfPatch(patchSeed);
    for (int i = 0; i < orderedUhs.size(); i++){
        auto uh = orderedUhs[i];
        //double uhScore = uhScores.at(uh);
        auto & affectedBhs = binariesDeterminedByOrderedUhs[i];
        if (affectedBhs.size() == 0){
            installationConfidenceOfOrderedUhs[i] = 1e5;
            continue;
        }

        double distanceSum = 0.0;
        for (auto bh : affectedBhs){
            distanceSum += AnchorDistanceSumOnBinary(mg, bh, patchSeed, vanishingPoints);
        }
        distanceSum /= affectedBhs.size();
        distanceSum /= avgDepthOfPatch;
        if (distanceSum == 0.0){
            distanceSum = 1e-5;
        }
        double uhConfidence = 1.0 - distanceSum;

        installationConfidenceOfOrderedUhs[i] = uhConfidence;
    }

    // check the most unconfident uhs
    std::vector<MGUnaryHandle> uhsOrderedByUnConfidencies = orderedUhs;
    std::sort(uhsOrderedByUnConfidencies.begin(), uhsOrderedByUnConfidencies.end(),
        [&uhOrders, &installationConfidenceOfOrderedUhs](MGUnaryHandle uh1, MGUnaryHandle uh2){
        return installationConfidenceOfOrderedUhs[uhOrders.at(uh1)] <
            installationConfidenceOfOrderedUhs[uhOrders.at(uh2)];
    });

    ManuallyOptimizeMixedGraph(panorama, mg, patchSeed, vanishingPoints,
        [&uhsOrderedByUnConfidencies, &uhOrders](MGUnaryHandle curUh) -> vis::Color {
        /*  auto pos = std::find(uhsOrderedByUnConfidencies.begin(), uhsOrderedByUnConfidencies.end(), curUh) -
        uhsOrderedByUnConfidencies.begin();
        double r = double(pos) / uhsOrderedByUnConfidencies.size();
        return vis::Color(r, r, r, 1.0);*/
        /* return vis::ColorFromTag(curUh == uhsOrderedByUnConfidencies.front() ?
        vis::ColorTag::Red : vis::ColorTag::Blue);*/
        /* double r = double(uhOrders[curUh]) / uhsOrderedByUnConfidencies.size();
        return vis::ColorFromHSV(r, 0.8, 0.8) / 255.0;*/

        bool v = false;
        for (int i = 0; i < 10; i++){
            if (curUh == uhsOrderedByUnConfidencies[i]){
                return vis::Color(i / 10.0, 0.0, 1.0 - i / 10.0);
            }
        }
        return vis::ColorTag::Red;

    }, false);
}

TEST(MixedGraph, Test2){

    using namespace core;
    std::vector<MGUnaryHandle> uhsOrderedByUnConfidencies;

    core::LoadFromDisk("./cache/uhsOrderedByUnConfidencies", uhsOrderedByUnConfidencies);


}