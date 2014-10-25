#define TEST_VIEWS_NET_INCREMENTED

#include <iostream>
#include <string>
#include <random>

#include <GCoptimization.h>

#include "../src/core/mesh_maker.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/debug.hpp"
#include "../src/rec/reconstruction_engine.hpp"
#include "../src/rec/regions_net_visualize.hpp"
#include "test_config.hpp"


using namespace panoramix;

void ShowPanoramaVPs(const rec::ReconstructionEngine & engine) {
    auto vps = engine.globalData().vanishingPoints;
    for (auto & vp : vps)
        vp /= core::norm(vp);
    double ortho = core::norm(core::Vec3(vps[0].dot(vps[1]), vps[1].dot(vps[2]), vps[2].dot(vps[0])));
    assert(ortho < 1e-1);

    auto antivps = vps;
    for (auto & p : antivps)
        p = -p;

    std::vector<core::Vec3> allvps(vps.begin(), vps.end());
    allvps.insert(allvps.end(), antivps.begin(), antivps.end());

    std::vector<core::Point2> vp2s(allvps.size());
    std::transform(allvps.begin(), allvps.end(), vp2s.begin(),
        [&engine](const core::Vec3 & p3){
        return engine.params().camera.screenProjection(p3);
    });

    vis::Visualizer2D(engine.globalData().panorama)
        << vis::manip2d::SetThickness(2)
        << vis::manip2d::SetColor(vis::ColorTag::Red) << vp2s[0]
        << vis::manip2d::SetColor(vis::ColorTag::Green) << vp2s[1]
        << vis::manip2d::SetColor(vis::ColorTag::Blue) << vp2s[2]
        << vis::manip2d::Show();
}

inline std::string ComposeCacheFileName(const std::string & inputFile, const std::string & tag) {
    std::string name = inputFile;
    for (auto & c : name){
        if (c == '/' || c == '\\')
            c = '.';
        if (c == ':')
            c = '@';
    }
    std::string d = "demo10.10";
    return "./cache/" + name + '.' + tag + d + ".state";
}

struct All {
    rec::ReconstructionEngine engine;
    core::Image panorama;
    core::PanoramicCamera originCam;
    std::vector<core::PerspectiveCamera> cams;
    std::vector<rec::ReconstructionEngine::ViewHandle> viewHandles;
    template <class Archiver>
    inline void serialize(Archiver & ar) {
        ar(engine, panorama, originCam, cams, viewHandles);
    }
};

int main(int argc, char * argv[], char * envp[]) {

    std::string originalFile = test::ProjectDataDirStrings::PanoramaIndoor + "/13.jpg";
    std::string cacheFileBeforeComputingFeatures = ComposeCacheFileName(originalFile, "1_before_fea");
    std::string cacheFileAfterComputingFeatures = ComposeCacheFileName(originalFile, "2_after_fea");
    std::string cacheFileAfterEstimatingVPs = ComposeCacheFileName(originalFile, "3_after_vp");
    std::string cacheFileAfterEstimatingLineDepths = ComposeCacheFileName(originalFile, "4_after_linedepths");

    All all;

    core::UpdateIfFileIsTooOld(originalFile, cacheFileBeforeComputingFeatures, 
        [&](const std::string & in, const std::string & out) {
        all.panorama = cv::imread(in);
        cv::resize(all.panorama, all.panorama, cv::Size(2000, 1000));
        all.engine.setPanorama(all.panorama);
        all.originCam = core::PanoramicCamera(all.panorama.cols / M_PI / 2.0);
        all.cams = {
            core::PerspectiveCamera(700, 700, all.originCam.focal(), { 0, 0, 0 }, { 1, 0, 0 }, { 0, 0, -1 }),
            core::PerspectiveCamera(700, 700, all.originCam.focal(), { 0, 0, 0 }, { 0, 1, 0 }, { 0, 0, -1 }),
            core::PerspectiveCamera(700, 700, all.originCam.focal(), { 0, 0, 0 }, { -1, 0, 0 }, { 0, 0, -1 }),
            core::PerspectiveCamera(700, 700, all.originCam.focal(), { 0, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 }),
            core::PerspectiveCamera(700, 700, all.originCam.focal(), { 0, 0, 0 }, { 0, 0, 1 }, { 1, 0, 0 }),
            core::PerspectiveCamera(700, 700, all.originCam.focal(), { 0, 0, 0 }, { 0, 0, -1 }, { 1, 0, 0 })
        };
        all.viewHandles.clear();
        all.viewHandles.reserve(all.cams.size());
        for (auto & camera : all.cams) {
            auto viewHandle = all.engine.insertPhoto(
                core::CameraSampler<core::PerspectiveCamera, core::PanoramicCamera>(camera, all.originCam)(all.panorama),
                camera);
            all.viewHandles.push_back(viewHandle);
            all.engine.updateConnections(viewHandle);
        }
        core::SaveToDisk(out, all);
    });

    core::UpdateIfFileIsTooOld(cacheFileBeforeComputingFeatures, cacheFileAfterComputingFeatures, 
        [&](const std::string & in, const std::string & out){
        core::LoadFromDisk(in, all);
        for (auto & viewHandle : all.viewHandles) {
            all.engine.computeFeatures(viewHandle);
        }        
        core::SaveToDisk(out, all);
    });
    
    core::UpdateIfFileIsTooOld(cacheFileAfterComputingFeatures, cacheFileAfterEstimatingVPs, 
        [&](const std::string & in, const std::string & out) {
        core::LoadFromDisk(in, all);
        all.engine.estimateVanishingPointsAndClassifyLines();
        core::SaveToDisk(out, all);
    });

    core::UpdateIfFileIsTooOld(cacheFileAfterEstimatingVPs, cacheFileAfterEstimatingLineDepths,
        [&](const std::string & in, const std::string & out) {
        core::LoadFromDisk(in, all);
        all.engine.recognizeRegionLineRelations();
        all.engine.estimateSpatialLineDepths();
        core::SaveToDisk(out, all);
    });

    core::UpdateIfFileIsTooOld(cacheFileAfterEstimatingLineDepths, "xxx",
        [&](const std::string & in, const std::string &) {
        core::LoadFromDisk(in, all);
        all.engine.estimateRegionPlanes();
    });


    //try {
    //    all.engine.initializeRegionOrientations();
    //} catch (GCException e) {
    //    e.Report();
    //}


    return 0;
}
