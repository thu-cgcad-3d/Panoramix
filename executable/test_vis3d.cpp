#include "../src/vis/visualize3d.hpp"
//#include "../src/vis/ogre_resources.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;
static const std::string ProjectTestDataDirStr_Normal = ProjectTestDataDirStr + "/normal";
static const std::string ProjectTestDataDirStr_PanoramaIndoor = ProjectTestDataDirStr + "/panorama/indoor";
static const std::string ProjectTestDataDirStr_PanoramaOutdoor = ProjectTestDataDirStr + "/panorama/outdoor";

TEST(Visualizer3D, 3D) { 

    vis::Visualizer3D()
        << vis::manip3d::SetDefaultColor(vis::ColorTag::Red)
        << core::Line3(core::Vec3(0, 0, 1), core::Vec3(0, 0, -1))
        << vis::manip3d::SetDefaultColor(vis::ColorTag::Yellow)
        << core::Line3(core::Vec3(-1, 0, 0), core::Vec3(1, 0, 0))
        << vis::manip3d::Show();

}

TEST(Visualizer3D, Texture) {
    std::vector<std::vector<std::pair<core::Point3, core::Point2>>> polys = {
            {
                { { -1.0, -1.0, -1.0 }, { 0.0, 0.0 } },
                { { -1.0, 1.0, -1.0 }, { 0.0, 1.0 } },
                { { 1.0, 1.0, -1.0 }, { 1.0, 1.0 } },
                { { 1.0, -1.0, -1.0 }, { 1.0, 0.0 } },
            },
            {
                { { -1.0, -1.0, -1.0 }, { 0.0, 0.0 } },
                { { -1.0, -1.0, 1.0 }, { 0.0, 1.0 } },
                { { 1.0, -1.0, 1.0 }, { 1.0, 1.0 } },
                { { 1.0, -1.0, -1.0 }, { 1.0, 0.0 } },
            },
            {
                { { -1.0, -1.0, 1.0 }, { 0.0, 0.0 } },
                { { -1.0, 1.0, 1.0 }, { 0.0, 1.0 } },
                { { 1.0, 1.0, 1.0 }, { 1.0, 1.0 } },
                { { 1.0, -1.0, 1.0 }, { 1.0, 0.0 } },
            },
            {
                { { -1.0, 1.0, -1.0 }, { 0.0, 0.0 } },
                { { -1.0, 1.0, 1.0 }, { 0.0, 1.0 } },
                { { 1.0, 1.0, 1.0 }, { 1.0, 1.0 } },
                { { 1.0, 1.0, -1.0 }, { 1.0, 0.0 } }
            }
    };

    cv::Mat texture = cv::imread(ProjectTestDataDirStr_PanoramaIndoor + "/14.jpg");
    cv::resize(texture, texture, cv::Size(1024, 512));
    vis::Visualizer3D() 
        << texture
        << polys 
        << vis::manip3d::AutoSetCamera 
        << vis::manip3d::Show();
//
}


TEST(Visualizer3D, Background) {

    for (auto c : vis::AllColorTags()) {
        std::stringstream ss;
        ss << c;
        vis::Visualizer3D()
            << vis::manip3d::SetBackgroundColor(c)
            << vis::manip3d::SetWindowName(ss.str())
            << vis::manip3d::Show(false);
    }

}

int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    //testing::FLAGS_gtest_filter = "Texture";
    return RUN_ALL_TESTS();
    //run();
    //return 0;
}
