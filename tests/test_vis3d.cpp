#include "../src/vis/visualize3d.hpp"
#include "../src/vis/ogre_resources.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;

TEST(Visualizer3D, Background){

    for (auto c : vis::AllColorTags()){
        std::stringstream ss;
        ss << c;
        vis::Visualizer3D()
            << vis::manip3d::SetBackgroundColor(c)
            << vis::manip3d::SetWindowName(ss.str())
            << vis::manip3d::Show(false);
    }

}

TEST(Visualizer3D, 3D) { 

    vis::Visualizer3D()
        << vis::manip3d::SetLineWidth(10)
        << vis::manip3d::SetDefaultColor(vis::ColorTag::Red)
        << core::Line3(core::Vec3(0, 0, 1), core::Vec3(0, 0, -1))
        << vis::manip3d::SetDefaultColor(vis::ColorTag::Yellow)
        << core::Line3(core::Vec3(-1, 0, 0), core::Vec3(1, 0, 0))
        << vis::manip3d::Show();

}

int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
