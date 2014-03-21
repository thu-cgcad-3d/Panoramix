#include "../src/core/views_net.hpp"
#include "../src/core/feature_visualize.hpp"
#include "../src/core/views_net_visualize.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;

TEST(ViewsNet, ViewsNet) {
    /// get perspective pieces

    cv::Mat panorama = cv::imread(ProjectTestDataDirStr + "/google0.png");

    cv::resize(panorama, panorama, cv::Size(1000, 500));
    core::PanoramicCamera originCam(panorama.cols / M_PI / 2.0);

    std::vector<core::PerspectiveCamera> cams(30);
    std::generate(cams.begin(), cams.end(), [](){
        return core::PerspectiveCamera(500, 500, 200, 
            core::PerspectiveCamera::Vec3(0, 0, 0), 
            core::PerspectiveCamera::Vec3((rand() % 10000)-5000, (rand() % 10000)-5000, (rand() % 10000-5000)),
            core::PerspectiveCamera::Vec3(0, 0, -1));
    });
    
    std::vector<core::Image> ims(cams.size());
    std::transform(cams.begin(), cams.end(), ims.begin(),
        [&panorama, &originCam](const core::PerspectiveCamera & pcam){
        core::Image im;
        core::CameraSampler<core::PerspectiveCamera, core::PanoramicCamera>(pcam, originCam)(panorama, im);
        return im;
    });

    /// insert into views net

    core::ViewsNet net;
    for (int i = 0; i < cams.size(); i++){
        auto & camera = cams[i];
        auto & im = ims[i];
        auto viewHandle = net.insertPhoto(im, camera);
        
        net.computeFeatures(viewHandle);
        core::ImageFeatureVisualizer(im) << 
            core::manip::SetColor(core::Color(0, 0, 255)) <<
            core::manip::SetThickness(2) <<
            net.views().data(viewHandle).lineSegments <<
            core::manip::SetColor(core::Color(255, 0, 0)) << 
            core::manip::SetThickness(1) <<
            net.views().data(viewHandle).lineSegmentIntersections;

        net.updateConnections(viewHandle);
        net.computeTransformationOnConnections(viewHandle);
        net.calibrateCamera(viewHandle);
        net.calibrateAllCameras();

        if (net.isTooCloseToAnyExistingView(viewHandle).isValid()){
            continue;
        }

        net.estimateVanishingPointsAndClassifyLines();
        auto vps = net.globalData().vanishingPoints;
        for (auto & vp : vps)
            vp /= cv::norm(vp);
        double ortho = cv::norm(core::Vec3(vps[0].dot(vps[1]), vps[1].dot(vps[2]), vps[2].dot(vps[0])));
        EXPECT_LT(ortho, 1e-1);

        auto antivps = vps;
        for (auto & p : antivps)
            p = -p;

        std::vector<core::Vec3> allvps(vps.begin(), vps.end());
        allvps.insert(allvps.end(), antivps.begin(), antivps.end());

        std::vector<core::Point2> vp2s(allvps.size());
        std::transform(allvps.begin(), allvps.end(), vp2s.begin(),
            [&originCam](const core::Vec3 & p3){
            return originCam.screenProjection(p3);
        });

        core::ImageFeatureVisualizer(panorama) <<
            core::manip::SetWindowName("Vanishing points") <<
            core::manip::SetColor(core::Color(0, 0, 255)) <<
            core::manip::SetThickness(3) <<
            vp2s;

        core::ImageFeatureVisualizer() << net.views().data(viewHandle);

    }

    
}

int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
