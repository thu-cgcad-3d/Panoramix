#define TEST_VIEWS_NET_INCREMENTED

#include "../src/core/mesh_maker.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/debug.hpp"
#include "../src/rec/reconstruction_engine.hpp"
#include "../src/rec/reconstruction_engine_visualize.hpp"
#include "../src/rec/regions_net_visualize.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

#include <Windows.h>


using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;
static const std::string ProjectTestDataDirStr_Normal = ProjectTestDataDirStr + "/normal";
static const std::string ProjectTestDataDirStr_PanoramaIndoor = ProjectTestDataDirStr + "/panorama/indoor";
static const std::string ProjectTestDataDirStr_PanoramaOutdoor = ProjectTestDataDirStr + "/panorama/outdoor";


int FilterException(int code, PEXCEPTION_POINTERS ex) {
    std::cout << "Filtering " << std::hex << code << std::endl;
    return EXCEPTION_EXECUTE_HANDLER;
};

void ShowPanoramaVPs(const rec::ReconstructionEngine & engine) {
    auto vps = engine.globalData().vanishingPoints;
    for (auto & vp : vps)
        vp /= core::norm(vp);
    double ortho = core::norm(core::Vec3(vps[0].dot(vps[1]), vps[1].dot(vps[2]), vps[2].dot(vps[0])));
    EXPECT_LT(ortho, 1e-1);

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

//TEST(ViewsNet, FixedCamera) {
void run(){

    cv::Mat panorama = cv::imread(ProjectTestDataDirStr_PanoramaIndoor + "/13.jpg");
    cv::resize(panorama, panorama, cv::Size(2000, 1000));
    core::PanoramicCamera originCam(panorama.cols / M_PI / 2.0);

    std::vector<core::PerspectiveCamera> cams = {
        core::PerspectiveCamera(700, 700, originCam.focal(), { 0, 0, 0 }, { 1, 0, 0 }, { 0, 0, -1 }),
        core::PerspectiveCamera(700, 700, originCam.focal(), { 0, 0, 0 }, { 0, 1, 0 }, { 0, 0, -1 }),
        core::PerspectiveCamera(700, 700, originCam.focal(), { 0, 0, 0 }, { -1, 0, 0 }, { 0, 0, -1 }),
        core::PerspectiveCamera(700, 700, originCam.focal(), { 0, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 })
    };

    rec::ReconstructionEngine engine;

    for (auto & camera : cams){
        auto viewHandle = engine.insertPhoto(
            core::CameraSampler<core::PerspectiveCamera, core::PanoramicCamera>(camera, originCam)(panorama), 
            camera);
        engine.computeFeatures(viewHandle);
        engine.updateConnections(viewHandle);
    }

    // estimate vanishing points and classify lines
    engine.estimateVanishingPointsAndClassifyLines();
    engine.reconstructLinesAndFaces();

}



int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    run();
    return 0;
}
