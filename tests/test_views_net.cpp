#include "../src/core/views_net.hpp"
#include "../src/core/feature_visualize.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;

TEST(ViewsNet, ViewsNet) {
    /// get perspective pieces

    cv::Mat panorama = cv::imread(ProjectTestDataDirStr + "/panofactory.jpg");
    ASSERT_EQ(4000, panorama.cols);
    ASSERT_EQ(2000, panorama.rows);

    cv::resize(panorama, panorama, cv::Size(1000, 500));
    core::PanoramicCamera originCam(panorama.cols / M_PI / 2.0);

    std::vector<core::PerspectiveCamera> cams(20);
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
            net.views().data(viewHandle).featureLineSegment <<
            core::manip::SetColor(core::Color(255, 0, 0)) << 
            core::manip::SetThickness(1) <<
            net.views().data(viewHandle).featureLineIntersections;

        net.updateConnections(viewHandle);
        net.updateConnections(viewHandle);
        net.calibrateCamera(viewHandle);
        net.calibrateAllCameras();

        net.computeGlobalFeatures();

    }
}

int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
