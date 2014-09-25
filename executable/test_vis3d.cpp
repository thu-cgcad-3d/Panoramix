#include "../src/vis/renderable_object_tree.hpp"
#include "../src/vis/visualize3d.hpp"
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

static_assert(vis::CanMakeRenderable<core::Line3>::value, "Line3 is not renderable!");
static_assert(vis::CanMakeRenderable<core::Point3>::value, "Point3 is not renderable!");

//TEST(Visualizer3D, RenderableObject) { 
void run() {
    std::vector<core::Line3> lines = {
        { {-1, 1}, {5, 5} },
        { { -2, -2 }, {8, 6} }
    };
    
    vis::Visualizer3D()
        << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Black)
        << lines
        << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
        << vis::manip3d::Show();
}

int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    //testing::FLAGS_gtest_filter = "Texture";
    //return RUN_ALL_TESTS();
    run();
    //return 0;
}
