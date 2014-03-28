#include "../src/vis/visualize3d.hpp"
#include "../src/vis/qt_resources.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;

TEST(Visualizer3D, Visualizer3D) {

    core::Line3 line1, line2, line3;
    line1.first = core::Point3(1, 0, 0);
    line1.second = core::Point3(-1, 0, 0);
    line2.first = { 0, 1, 0 };
    line2.second = { 0, -1, 0 };
    line3.first = { 0, 0, 1 };
    line3.second = { 0, 0, -1 };
    vis::Visualizer3D()        
        << line1
        << core::Point3(0, 0, 1)
        << core::Point3(0, 0, -1)
        //<< line1
        << vis::manip3d::SetDefaultColor(core::Color(255,0,0))
        << line2
        << line3
        << vis::manip3d::AutoSetCamera
        << vis::manip3d::SetRenderMode(vis::All)
        << vis::manip3d::Show();
    
}

int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    vis::InitGui(argc, argv);
    return RUN_ALL_TESTS();
}
